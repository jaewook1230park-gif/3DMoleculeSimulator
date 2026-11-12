import streamlit as st
import streamlit.components.v1 as components
import requests
import pandas as pd
import re
import json
from urllib.parse import quote
from rdkit import Chem
from rdkit.Chem import Descriptors, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

st.set_page_config(
    page_title="분자 탐험실",
    page_icon="🧪",
    layout="wide"
)

st.markdown("""
<style>
@import url('https://fonts.googleapis.com/css2?family=Nanum+Gothic:wght@400;700;800&display=swap');
html, body, [class*="css"] {
    font-family: 'Nanum Gothic', sans-serif;
}
.stApp {
    background: linear-gradient(135deg,#fffaf0,#f5f8ff);
}
.hero {
    padding: 1.4rem 1.6rem;
    border-radius: 22px;
    background: linear-gradient(110deg,#fff0c2,#e9f3ff);
    border: 1px solid #f1dfad;
}
.hero h1 {
    margin: 0;
    color: #24324a;
    font-size: 2.2rem;
}
.hero p {
    color: #53627a;
}
.smallnote {
    color: #64748b;
    font-size: .85rem;
}
</style>

<div class="hero">
  <h1>🧪 분자 탐험실</h1>
  <p>분자의 모양과 물성을 탐험해 보세요.</p>
</div>
""", unsafe_allow_html=True)

PUBCHEM = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
HEADERS = {
    "User-Agent": "MoleculeExplorerStreamlit/1.0"
}


@st.cache_data(ttl=3600, show_spinner=False)
def pubchem_json(url):
    response = requests.get(
        url, timeout=20, headers=HEADERS
    )
    response.raise_for_status()
    return response.json()


def get_cid(query):
    url = (
        f"{PUBCHEM}/compound/name/"
        f"{quote(query, safe='')}/cids/JSON"
    )

    try:
        data = pubchem_json(url)
        return int(data["IdentifierList"]["CID"][0])
    except Exception:
        url = (
            f"{PUBCHEM}/compound/smiles/"
            f"{quote(query, safe='')}/cids/JSON"
        )
        data = pubchem_json(url)
        return int(data["IdentifierList"]["CID"][0])


def get_properties(cid):
    props = (
        "Title,IUPACName,MolecularFormula,"
        "CanonicalSMILES,IsomericSMILES,InChI"
    )
    url = (
        f"{PUBCHEM}/compound/cid/{cid}/property/"
        f"{props}/JSON"
    )
    return pubchem_json(url)["PropertyTable"]["Properties"][0]


def recursive_sections(obj):
    if isinstance(obj, dict):
        heading = obj.get("TOCHeading", "")

        for info in obj.get("Information", []):
            value = info.get("Value", {})
            strings = []

            if isinstance(value, dict):
                for item in value.get("StringWithMarkup", []):
                    if item.get("String"):
                        strings.append(item["String"])

                for item in value.get("Number", []):
                    if "Number" in item:
                        number = str(item["Number"])
                        unit = item.get("Unit", "")
                        strings.append(
                            number + (" " + unit if unit else "")
                        )

            if strings:
                yield heading, " / ".join(strings)

        for child in obj.get("Section", []):
            yield from recursive_sections(child)

    elif isinstance(obj, list):
        for item in obj:
            yield from recursive_sections(item)


@st.cache_data(ttl=3600, show_spinner=False)
def get_pubchem_record(cid):
    url = (
        "https://pubchem.ncbi.nlm.nih.gov/rest/"
        f"pug_view/data/compound/{cid}/JSON"
    )
    try:
        data = pubchem_json(url)
        return list(
            recursive_sections(
                data.get("Record", {}).get("Section", [])
            )
        )
    except Exception:
        return []


def find_property(sections, patterns):
    for heading, value in sections:
        if any(p in heading.lower() for p in patterns):
            if value.strip():
                return value.strip()
    return "자료 없음"


def mol_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError("분자 구조를 해석하지 못했습니다.")
    return mol
def svg_for_molecule(mol, legend=""):
    molecule = Chem.Mol(mol)
    Chem.rdDepictor.Compute2DCoords(molecule)

    drawer = rdMolDraw2D.MolDraw2DSVG(650, 400)
    options = drawer.drawOptions()
    options.addStereoAnnotation = True
    options.bondLineWidth = 2.2
    options.padding = 0.08

    drawer.DrawMolecule(molecule, legend=legend)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def get_3d_sdf(cid):
    url = (
        f"{PUBCHEM}/compound/cid/{cid}/SDF"
        "?record_type=3d"
    )
    response = requests.get(
        url, timeout=25, headers=HEADERS
    )
    response.raise_for_status()

    if "V2000" not in response.text and "V3000" not in response.text:
        raise ValueError("3D 좌표가 제공되지 않습니다.")

    return response.text


def show_3d(sdf_text, mode="ballstick"):
    safe_sdf = json.dumps(sdf_text)

    if mode == "spacefill":
        style = "{sphere:{scale:0.34}}"
    else:
        style = "{stick:{radius:0.16},sphere:{scale:0.25}}"

    html = f"""
    <div id="mol3d"
         style="width:100%;height:440px;
                border-radius:14px;background:#fbfdff;">
    </div>

    <script src="https://3Dmol.org/build/3Dmol-min.js"></script>
    <script>
    const sdf = {safe_sdf};
    const viewer = $3Dmol.createViewer(
        document.getElementById('mol3d'),
        {{backgroundColor:'#fbfdff'}}
    );
    viewer.addModel(sdf, 'sdf');
    viewer.setStyle({{}}, {style});
    viewer.zoomTo();
    viewer.render();
    </script>
    """

    components.html(html, height=455, scrolling=False)


def render_structure(mol, legend=""):
    components.html(
        svg_for_molecule(mol, legend),
        height=420,
        scrolling=False
    )


st.markdown("### 🔎 어떤 분자를 알아볼까요?")

with st.form("molecule_form"):
    query = st.text_input(
        "분자 이름 또는 SMILES",
        value="물",
        placeholder="물, 에탄올, 포도당, benzene, CCO"
    )
    submitted = st.form_submit_button(
        "분자 구조 찾기",
        type="primary",
        use_container_width=True
    )

if "molecule_query" not in st.session_state:
    st.session_state.molecule_query = "물"

if submitted:
    st.session_state.molecule_query = query.strip()

query = st.session_state.molecule_query
if query:
    try:
        with st.spinner("분자 정보를 가져오는 중입니다..."):
            cid = get_cid(query)
            p = get_properties(cid)
            sections = get_pubchem_record(cid)

        smiles = (
            p.get("IsomericSMILES")
            or p.get("CanonicalSMILES")
        )
        mol = mol_from_smiles(smiles)

        formula = p.get("MolecularFormula", "자료 없음")
        title = p.get("Title", query)
        iupac = p.get("IUPACName", "자료 없음")
        mass = Descriptors.MolWt(mol)

        # 실험식 계산
        from collections import Counter
        from functools import reduce
        from math import gcd

        mol_h = Chem.AddHs(mol)
        counts = Counter(
            atom.GetSymbol() for atom in mol_h.GetAtoms()
        )
        divisor = reduce(gcd, counts.values()) if counts else 1

        order = ["C", "H"] + sorted(
            e for e in counts if e not in ("C", "H")
        )

        empirical_parts = []
        for element in order:
            n = counts.get(element, 0)
            if n:
                reduced = n // divisor
                empirical_parts.append(
                    element + (str(reduced) if reduced != 1 else "")
                )

        empirical_formula = "".join(empirical_parts)

        melting = find_property(
            sections, ["melting point"]
        )
        boiling = find_property(
            sections, ["boiling point"]
        )
        density = find_property(
            sections, ["density"]
        )
        flash = find_property(
            sections, ["flash point"]
        )

        if flash == "자료 없음":
            flash = "없음"

        st.markdown("---")
        st.markdown(f"## {title}")
        st.caption(f"PubChem CID: {cid}")

        c1, c2, c3 = st.columns(3)
        c1.metric("분자식", formula)
        c2.metric("화학식량", f"{mass:.3f} g/mol")
        c3.metric("원자 수", str(mol_h.GetNumAtoms()))

        st.markdown("### 🧬 분자 구조")

        tabs = st.tabs([
            "루이스 구조식",
            "공-막대 모형",
            "공간 채움 모형",
            "골격 구조식",
            "회전 가능한 3D"
        ])

        with tabs[0]:
            st.info(
                "현재는 2D 결합 구조를 표시합니다. "
                "비공유 전자쌍을 정확히 그리는 전용 루이스 "
                "전자점 표기는 추가 구현이 필요합니다."
            )
            render_structure(mol, "2D 결합 구조")

        with tabs[1]:
            try:
                show_3d(get_3d_sdf(cid), "ballstick")
            except Exception:
                st.warning("3D 좌표를 불러오지 못했습니다.")
                render_structure(mol)

        with tabs[2]:
            try:
                show_3d(get_3d_sdf(cid), "spacefill")
            except Exception:
                st.warning("공간 채움 모형을 표시할 수 없습니다.")

        with tabs[3]:
            if any(a.GetSymbol() == "C" for a in mol.GetAtoms()):
                st.caption(
                    "탄소를 포함하는 물질의 2D 구조입니다. "
                    "정식 골격 구조식 표기와는 차이가 있을 수 있습니다."
                )
                render_structure(mol, "골격 구조")
            else:
                st.info(
                    "탄소를 포함하지 않아 골격 구조식을 "
                    "표시하지 않습니다."
                )

        with tabs[4]:
            st.caption("마우스로 드래그하여 회전할 수 있습니다.")
            try:
                show_3d(get_3d_sdf(cid), "ballstick")
            except Exception:
                st.warning("회전 가능한 3D 좌표가 없습니다.")

        # 일부 대표 물질에 대한 시성식 예시
        condensed = {
            962: "H₂O",
            702: "CH₃CH₂OH",
            180: "CH₃COOH",
            176: "CH₃COCH₃",
            241: "C₆H₆ (벤젠 고리)"
        }.get(cid, "자료 없음")

        st.markdown("---")
        st.markdown("### 📋 분자 정보")

        info_rows = [
            ("분자식", formula),
            ("시성식", condensed),
            ("실험식", empirical_formula),
            ("화학식량", f"{mass:.3f} g/mol"),
            ("녹는점", melting),
            ("끓는점", boiling),
            ("밀도 (25°C)", density),
            ("인화점", flash),
        ]

        st.dataframe(
            pd.DataFrame(info_rows, columns=["항목", "값"]),
            hide_index=True,
            use_container_width=True
        )

        st.markdown("### 📝 참고")
        st.markdown(f"""
- **분자명:** {title}
- **IUPAC 명칭:** {iupac}
- **SMILES:** `{smiles}`
- **데이터 출처:** [PubChem](https://pubchem.ncbi.nlm.nih.gov/compound/{cid})
        """)

        st.caption(
            "교육용 참고 도구입니다. 물성값은 물질의 상태와 "
            "측정 조건에 따라 달라질 수 있습니다."
        )

    except requests.HTTPError:
        st.error(
            "분자를 찾지 못했거나 PubChem 응답에 문제가 있습니다. "
            "영문 분자명 또는 SMILES를 바꿔 입력해 보세요."
        )
    except requests.RequestException:
        st.error(
            "인터넷 또는 PubChem 서버에 연결할 수 없습니다."
        )
    except Exception as e:
        st.error(f"분자 정보를 처리하지 못했습니다: {e}")
        st.info("예시: 물, 에탄올, 아세트산, benzene, CCO")

st.markdown("---")
st.caption(
    "분자 구조는 RDKit 및 PubChem 데이터를 사용합니다. "
    "3D 모형은 3Dmol.js를 사용하므로 인터넷 연결이 필요합니다."
)
