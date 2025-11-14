import streamlit as st
import requests
import pandas as pd
import base64
from Bio import Entrez
from Bio.Blast import NCBIWWW, NCBIXML
from html import escape
from io import StringIO, BytesIO

# New imports for coordinate analysis
from Bio.PDB import PDBParser, MMCIFParser
from Bio.PDB.Polypeptide import is_aa
import numpy as np
import matplotlib.pyplot as plt
from scipy.spatial.distance import pdist, squareform

# ---------- App config ----------
st.set_page_config(page_title="PDB Explorer", layout="wide")
st.title("GeneScope — PDB Explorer")
st.write("Structure-first pipeline using the RCSB Data API. Enter a PDB ID or upload a PDB/mmCIF file.")

# ---------- Simple CSS ----------
st.markdown(
    """
    <style>
      body { background:#f6f8fb; }
      .card { background:white; padding:14px; border-radius:8px; box-shadow:0 1px 6px rgba(16,24,40,0.06); margin-bottom:16px; }
      .muted { color:#6b7280; font-size:13px; }
    </style>
    """,
    unsafe_allow_html=True,
)

# ---------- Sidebar: menu + settings ----------
st.sidebar.header("Controls")
menu = st.sidebar.radio(
    "Page",
    options=[
        "Overview",
        "Visualization",
        "Chains",
        "Ligands",
        "Experimental",
        "Validation",
        "BLAST",
        "Downloads",
        "Coordinate Analysis"   # new page
    ],
)

st.sidebar.markdown("---")
email = st.sidebar.text_input("NCBI email (required for BLAST)", value="")
if email:
    Entrez.email = email
st.sidebar.caption("Enter a valid email for NCBI qblast. No email -> BLAST blocked.")

# ---------- Input: PDB ID or upload ----------
with st.container():
    col1, col2 = st.columns([3, 1])
    with col1:
        pdb_input = st.text_input("PDB ID (e.g. 4O4S)", value="").strip().upper()
    with col2:
        uploaded = st.file_uploader("Or upload PDB/mmCIF", type=["pdb", "cif", "ent"])
    fetch_btn = st.button("Fetch / Load")

# session store for loaded data
if "entry_json" not in st.session_state:
    st.session_state.entry_json = None
if "struct_summary" not in st.session_state:
    st.session_state.struct_summary = None
if "poly_seqs" not in st.session_state:
    st.session_state.poly_seqs = {}  # entity_id -> sequence
if "uploaded_b64" not in st.session_state:
    st.session_state.uploaded_b64 = None
if "pdb_source" not in st.session_state:
    st.session_state.pdb_source = None  # "rcsb" or "upload"
if "selected_entity" not in st.session_state:
    st.session_state.selected_entity = None
if "parsed_structure" not in st.session_state:
    st.session_state.parsed_structure = None  # Bio.PDB Structure object

# ---------- Helper: RCSB endpoints ----------
def rcsb_get_json(url):
    try:
        r = requests.get(url, timeout=20)
        r.raise_for_status()
        return r.json()
    except Exception:
        return None

def fetch_entry(pdb_id):
    return rcsb_get_json(f"https://data.rcsb.org/rest/v1/core/entry/{pdb_id}")

def fetch_struct_summary(pdb_id):
    return rcsb_get_json(f"https://data.rcsb.org/rest/v1/core/structure_summary/{pdb_id}")

def fetch_polymer_entity(pdb_id, entity_id):
    return rcsb_get_json(f"https://data.rcsb.org/rest/v1/core/polymer_entity/{pdb_id}/{entity_id}")

def fetch_chem_comp(chem_id):
    return rcsb_get_json(f"https://data.rcsb.org/rest/v1/core/chem_comp/{chem_id}")

def fetch_validation(pdb_id):
    return rcsb_get_json(f"https://data.rcsb.org/rest/v1/core/validation_summary/{pdb_id}")

def fetch_experiment(pdb_id):
    return rcsb_get_json(f"https://data.rcsb.org/rest/v1/core/experiment/{pdb_id}")

# ---------- Load logic ----------
if fetch_btn:
    # clear previous
    st.session_state.entry_json = None
    st.session_state.struct_summary = None
    st.session_state.poly_seqs = {}
    st.session_state.uploaded_b64 = None
    st.session_state.pdb_source = None
    st.session_state.selected_entity = None
    st.session_state.parsed_structure = None

    if uploaded is not None:
        # handle uploaded file: read bytes and create base64 data URI for NGL
        content = uploaded.getvalue()
        b64 = base64.b64encode(content).decode("utf-8")
        st.session_state.uploaded_b64 = b64
        st.session_state.pdb_source = "upload"
        st.success("Loaded uploaded structure file (viewer will use the uploaded file).")
        # attempt to parse uploaded file into a Bio.PDB structure
        try:
            text = content.decode("utf-8", errors="ignore")
            # choose parser by extension
            if uploaded.name.lower().endswith(".cif"):
                parser = MMCIFParser()
                structure = parser.get_structure("uploaded", StringIO(text))
            else:
                parser = PDBParser(QUIET=True)
                structure = parser.get_structure("uploaded", StringIO(text))
            st.session_state.parsed_structure = structure
            # extract SEQRES-based chains if present (simple)
            chains = {}
            for line in text.splitlines():
                if line.startswith("SEQRES"):
                    parts = line.split()
                    if len(parts) >= 5:
                        chain = parts[2]
                        seqs = parts[4:]
                        chains.setdefault(chain, []).extend(seqs)
            for c, s in chains.items():
                st.session_state.poly_seqs[c] = " ".join(s)
        except Exception:
            # parsing failed — still allow visualization via base64
            st.warning("Uploaded file parsed partially. Visualization should still work.")
    elif pdb_input:
        # fetch RCSB JSON and polymer entities/sequences
        entry = fetch_entry(pdb_input)
        if not entry:
            st.error("PDB ID not found or RCSB API error.")
        else:
            st.session_state.entry_json = entry
            st.session_state.pdb_source = "rcsb"
            st.success(f"Fetched RCSB entry for {pdb_input}.")
            # fetch structure summary and try to collect polymer entities
            st.session_state.struct_summary = fetch_struct_summary(pdb_input)
            # collect polymer entity ids
            cids = entry.get("rcsb_entry_container_identifiers", {}).get("polymer_entity_ids", [])
            if not cids:
                # fallback: try polymer_entities list
                cids = [str(p.get("entity_id")) for p in entry.get("polymer_entities", []) if p.get("entity_id")]
            # fetch sequences for entities
            for eid in cids:
                try:
                    poly = fetch_polymer_entity(pdb_input, eid)
                    seq = poly.get("entity_poly", {}).get("pdbx_seq_one_letter_code_can") or poly.get("entity_poly", {}).get("pdbx_seq_one_letter_code", "")
                    if seq:
                        seq = "".join(seq.split())
                        st.session_state.poly_seqs[str(eid)] = seq
                except Exception:
                    continue
            # attempt to download coordinate file and parse for analysis (optional)
            try:
                pdb_url = f"https://files.rcsb.org/download/{pdb_input}.pdb"
                r = requests.get(pdb_url, timeout=20)
                if r.status_code == 200:
                    text = r.text
                    parser = PDBParser(QUIET=True)
                    structure = parser.get_structure(pdb_input, StringIO(text))
                    st.session_state.parsed_structure = structure
            except Exception:
                # parsing the coordinate file failed; keep going
                st.session_state.parsed_structure = None

# ---------- small helper to ensure dataframe arrow compatibility ----------
def safe_table(df):
    # convert all columns to strings to avoid pyarrow issues
    try:
        df2 = df.copy()
        for col in df2.columns:
            df2[col] = df2[col].astype(str)
        st.dataframe(df2, use_container_width=True)
    except Exception:
        st.table(df)

# ---------- Utility: build NGL HTML (works with RCSB URL or uploaded base64) ----------
def ngl_viewer_html_from_url(url, show_hetero=True, highlight_resname=None):
    # url is a direct fetchable PDB file url or data URL
    sel_high = f'var highlightSel = "resname {highlight_resname}";' if highlight_resname else "var highlightSel = null;"
    hetero_repr = 'o.addRepresentation("ball+stick", {sele: "hetero", color: "element"});' if show_hetero else ""
    html = f"""
    <div id="viewport" style="width:100%; height:600px; border-radius:8px; overflow:hidden;"></div>
    <script src="https://unpkg.com/ngl@0.10.4/dist/ngl.js"></script>
    <script>
        {sel_high}
        var stage = new NGL.Stage("viewport");
        // try loading file
        stage.loadFile("{url}", {{defaultRepresentation: true}}).then(function(o) {{
            o.addRepresentation("cartoon", {{ sele: "protein", color: "chainname" }});
            {hetero_repr}
            if (highlightSel) {{
                o.addRepresentation("ball+stick", {{ sele: highlightSel, color: "element", scale: 2 }});
            }}
            o.autoView();
        }}).catch(function(e) {{
            var el = document.getElementById("viewport");
            el.innerHTML = "<div style='padding:10px;color:#b91c1c'>Viewer failed to load. Open file URL in new tab.</div>";
        }});
        window.addEventListener("resize", function() {{ stage.handleResize(); }}, false);
    </script>
    """
    return html

# ---------- Coordinate analysis module ----------
def analyze_structure_and_show(structure, max_ca_for_matrix=500):
    """
    Analyze a Bio.PDB Structure object and render tables/plots in Streamlit.
    Returns a pandas DataFrame of atom-level records.
    """
    atoms = []
    lig_atoms = []
    # iterate atoms
    for model in structure:
        for chain in model:
            chain_id = chain.id
            for residue in chain:
                # residue.id is a tuple: (hetfield, resseq, icode)
                resname = residue.get_resname()
                res_id = residue.id[1]
                is_lig = not is_aa(residue, standard=True)
                for atom in residue.get_atoms():
                    coord = atom.get_coord()
                    b = atom.get_bfactor()
                    atoms.append({
                        "model": model.id,
                        "chain": chain_id,
                        "residue": resname,
                        "res_id": res_id,
                        "atom": atom.get_name(),
                        "x": float(coord[0]),
                        "y": float(coord[1]),
                        "z": float(coord[2]),
                        "b_factor": float(b),
                        "is_ligand": bool(is_lig)
                    })
                    if is_lig:
                        lig_atoms.append({
                            "chain": chain_id,
                            "residue": resname,
                            "res_id": res_id,
                            "atom": atom.get_name(),
                            "x": float(coord[0]),
                            "y": float(coord[1]),
                            "z": float(coord[2]),
                            "b_factor": float(b)
                        })
    if not atoms:
        st.info("No atoms parsed from the structure.")
        return None

    df_atoms = pd.DataFrame(atoms)

    # Show basic atom/residue table (first N rows)
    st.subheader("Atom table (first 500 rows)")
    safe_table(df_atoms.head(500))

    # Provide CSV download of full atom table
    csv_buf = df_atoms.to_csv(index=False)
    st.download_button("Download full atom table (CSV)", data=csv_buf, file_name=f"{pdb_input or 'uploaded'}_atoms.csv", mime="text/csv")

    # B-factor distribution
    st.subheader("B-factor distribution (all atoms)")
    fig, ax = plt.subplots()
    ax.hist(df_atoms["b_factor"].dropna(), bins=40)
    ax.set_xlabel("B-factor")
    ax.set_ylabel("Frequency")
    st.pyplot(fig)

    # Ligand-specific analysis
    if len(lig_atoms) > 0:
        st.subheader("Ligand atoms (summary)")
        df_lig = pd.DataFrame(lig_atoms)
        safe_table(df_lig.head(200))
        st.write("Ligand B-factor distribution")
        fig2, ax2 = plt.subplots()
        ax2.hist(df_lig["b_factor"].dropna(), bins=30)
        ax2.set_xlabel("B-factor (ligand atoms)")
        ax2.set_ylabel("Frequency")
        st.pyplot(fig2)
    else:
        st.info("No ligand atoms detected in the parsed structure.")

    # Cα distance matrix (limit size to avoid heavy plotting)
    ca_atoms = df_atoms[df_atoms["atom"] == "CA"]
    n_ca = len(ca_atoms)
    if n_ca > 1:
        st.subheader(f"Cα pairwise distances (n={n_ca})")
        if n_ca > max_ca_for_matrix:
            st.warning(f"Cα count {n_ca} exceeds {max_ca_for_matrix}. Skipping full distance heatmap to avoid browser lag.")
        else:
            coords = ca_atoms[["x","y","z"]].to_numpy()
            D = squareform(pdist(coords))
            fig3, ax3 = plt.subplots(figsize=(6,5))
            im = ax3.imshow(D, cmap="viridis")
            ax3.set_title("Cα pairwise distances (Å)")
            fig3.colorbar(im, ax=ax3)
            st.pyplot(fig3)

    # Chain-level statistics: residue count vs avg B-factor
    st.subheader("Chain-level stats")
    try:
        chain_stats = df_atoms.groupby("chain").agg({"res_id": "nunique", "b_factor": "mean"}).rename(columns={"res_id":"Residue Count","b_factor":"Avg B-factor"})
        safe_table(chain_stats.reset_index().rename(columns={"index":"Chain"}))
        fig4, ax4 = plt.subplots()
        ax4.scatter(chain_stats["Residue Count"], chain_stats["Avg B-factor"])
        ax4.set_xlabel("Residue count")
        ax4.set_ylabel("Average B-factor")
        ax4.set_title("Chain size vs Avg B-factor")
        st.pyplot(fig4)
    except Exception:
        st.info("Failed to compute chain-level stats.")

    return df_atoms

# ---------- Page: Overview ----------
if menu == "Overview":
    st.header("Overview")
    st.markdown("Main structure-level fields from RCSB (summary).")
    if st.session_state.pdb_source == "rcsb" and st.session_state.entry_json:
        entry = st.session_state.entry_json
        # try to gather straightforward summary fields
        title = entry.get("struct", {}).get("title", "")
        methods = ", ".join([e.get("method", "") for e in entry.get("exptl", [])]) if entry.get("exptl") else ""
        release = entry.get("rcsb_accession_info", {}).get("initial_release_date", "")
        deposition = entry.get("rcsb_accession_info", {}).get("deposit_date", "")
        rev = entry.get("rcsb_accession_info", {}).get("revision_date", "")
        pdb_info = [
            ["PDB ID", pdb_input],
            ["Title", title],
            ["Experimental Method(s)", methods],
            ["Release date", release],
            ["Deposition date", deposition],
            ["Revision date", rev],
            ["Source", "RCSB Data API"]
        ]
        df = pd.DataFrame(pdb_info, columns=["Field", "Value"]).astype(str)
        safe_table(df)
        # quick summary from structure_summary if available
        if st.session_state.struct_summary:
            ss = st.session_state.struct_summary
            try:
                quick = [
                    ["Authors", ", ".join(ss.get("rcsb_authors", []) )],
                    ["Polymer entity count", str(ss.get("polymer_entity_count",""))],
                    ["Non-polymer entity count", str(ss.get("nonpolymer_entity_count",""))],
                    ["Molecular weight (Da)", str(ss.get("molecular_weight",""))],
                ]
                safe_table(pd.DataFrame(quick, columns=["Metric","Value"]))
            except Exception:
                pass
    elif st.session_state.pdb_source == "upload":
        st.info("Overview limited for uploaded files. Viewer and basic parsing available.")
        if st.session_state.uploaded_b64:
            st.write("Uploaded structure is available to visualize and to extract limited chain info.")
    else:
        st.info("Provide a PDB ID (fetch) or upload a PDB/mmCIF file, then click Fetch/Load.")

# ---------- Page: Visualization ----------
elif menu == "Visualization":
    st.header("Visualization")
    st.markdown("Interactive 3D view. You can highlight ligands (if any).")
    # determine URL or data URI
    highlight = None
    if st.session_state.pdb_source == "rcsb" and st.session_state.entry_json:
        url = f"https://files.rcsb.org/download/{pdb_input}.pdb"
        # collect ligands to allow highlight select
        nonpoly = st.session_state.entry_json.get("nonpolymer_entities", {})
        ligand_ids = []
        if isinstance(nonpoly, dict):
            for v in nonpoly.values():
                chem = v.get("nonpolymer_comp", {}).get("chem_comp", {})
                if chem.get("id"):
                    ligand_ids.append(chem.get("id"))
        elif isinstance(nonpoly, list):
            for v in nonpoly:
                chem = v.get("nonpolymer_comp", {}).get("chem_comp", {})
                if chem.get("id"):
                    ligand_ids.append(chem.get("id"))
        ligand_ids = sorted(set(ligand_ids))
        if ligand_ids:
            highlight = st.selectbox("Highlight ligand (resname)", options=["None"] + ligand_ids)
            if highlight == "None":
                highlight = None
        ngl_html = ngl_viewer_html_from_url(url, show_hetero=True, highlight_resname=highlight)
        st.components.v1.html(ngl_html, height=640, scrolling=True)
    elif st.session_state.pdb_source == "upload" and st.session_state.uploaded_b64:
        # create data URI for uploaded file; assume PDB text
        b64 = st.session_state.uploaded_b64
        data_url = f"data:text/plain;base64,{b64}"
        highlight = st.text_input("Highlight ligand resname (e.g. HEM) — optional", value="")
        highlight = highlight.strip().upper() if highlight.strip() else None
        ngl_html = ngl_viewer_html_from_url(data_url, show_hetero=True, highlight_resname=highlight)
        st.components.v1.html(ngl_html, height=640, scrolling=True)
    else:
        st.info("Fetch a PDB ID or upload a structure file to visualize.")

# ---------- Page: Chains ----------
elif menu == "Chains":
    st.header("Macromolecules / Chains")
    st.markdown("Table of polymer entities and their sequences (if available). You can download FASTA of a selected chain.")
    if st.session_state.poly_seqs:
        # build table
        rows = []
        for eid, seq in st.session_state.poly_seqs.items():
            rows.append([str(eid), str(len(seq)), seq[:60] + ("..." if len(seq) > 60 else "")])
        df = pd.DataFrame(rows, columns=["Entity ID", "Length", "Sequence (preview)"])
        safe_table(df)
        sel = st.selectbox("Select entity to view/download", options=list(st.session_state.poly_seqs.keys()))
        seq_full = st.session_state.poly_seqs.get(sel, "")
        st.code(seq_full[:5000] + ("..." if len(seq_full) > 5000 else ""))
        fasta = f">{pdb_input}_entity_{sel}\n" + "\n".join([seq_full[i:i+70] for i in range(0, len(seq_full), 70)])
        st.download_button("Download chain as FASTA", data=fasta, file_name=f"{pdb_input}_entity_{sel}.fasta", mime="text/plain")
        st.session_state.selected_entity = sel
    else:
        st.info("No chain sequences available. Fetch a PDB entry from RCSB to extract polymer sequences.")

# ---------- Page: Ligands ----------
elif menu == "Ligands":
    st.header("Small molecules / Ligands")
    st.markdown("List of ligands from the entry. Click a ligand to view details and highlight it in the 3D viewer.")
    if st.session_state.pdb_source == "rcsb" and st.session_state.entry_json:
        nonpoly = st.session_state.entry_json.get("nonpolymer_entities", {})
        items = []
        if isinstance(nonpoly, dict):
            items = list(nonpoly.values())
        elif isinstance(nonpoly, list):
            items = nonpoly
        ligand_rows = []
        ligand_ids = []
        for lig in items:
            chem = lig.get("nonpolymer_comp", {}).get("chem_comp", {})
            chem_id = chem.get("id", "")
            name = chem.get("name", "")
            formula = chem.get("formula", "")
            ligand_rows.append([chem_id, name, formula])
            if chem_id:
                ligand_ids.append(chem_id)
        if ligand_rows:
            df_lig = pd.DataFrame(ligand_rows, columns=["Ligand ID", "Name", "Formula"])
            safe_table(df_lig)
            chosen = st.selectbox("Select ligand to fetch details", options=["None"] + sorted(set(ligand_ids)))
            if chosen and chosen != "None":
                comp = fetch_chem_comp(chosen)
                if comp:
                    # show main fields in table
                    fields = [
                        ["Chem ID", comp.get("id","")],
                        ["Name", comp.get("name","")],
                        ["Formula", comp.get("formula","")],
                        ["SMILES", comp.get("chemical_fragments", {}).get("smiles","") if comp.get("chemical_fragments") else ""],
                        ["InChI", comp.get("chemical_fragments", {}).get("inchi","") if comp.get("chemical_fragments") else ""]
                    ]
                    safe_table(pd.DataFrame(fields, columns=["Field","Value"]))
                    # provide highlight option: open visualization page to highlight resname
                    if st.button("Highlight this ligand in viewer"):
                        url = f"https://files.rcsb.org/download/{pdb_input}.pdb" if st.session_state.pdb_source == "rcsb" else f"data:text/plain;base64,{st.session_state.uploaded_b64}"
                        ngl_html = ngl_viewer_html_from_url(url, show_hetero=True, highlight_resname=chosen)
                        st.components.v1.html(ngl_html, height=640, scrolling=True)
                else:
                    st.info("No chem_comp details from RCSB for this ligand.")
        else:
            st.info("No ligands listed in the RCSB entry JSON.")
    else:
        st.info("Ligand details require a fetched RCSB entry. Uploads have limited ligand parsing support.")

# ---------- Page: Experimental ----------
elif menu == "Experimental":
    st.header("Experimental details")
    st.markdown("Show basic experiment and method metadata from RCSB.")
    if st.session_state.pdb_source == "rcsb" and st.session_state.entry_json:
        expt = fetch_experiment(pdb_input)
        if expt:
            # expt can be complex; pick a few standard fields
            try:
                method = ", ".join([e.get("method","") for e in expt.get("exptl", [])]) if expt.get("exptl") else st.session_state.entry_json.get("exptl", [{}])[0].get("method","")
                sample_prep = expt.get("sample_preparation", "N/A")
                instrument = expt.get("instrument", "N/A")
                ex_rows = [
                    ["Method", method],
                    ["Sample preparation", sample_prep],
                    ["Instrument", instrument]
                ]
                safe_table(pd.DataFrame(ex_rows, columns=["Field","Value"]))
            except Exception:
                st.write("Experiment info present but could not be parsed fully.")
        else:
            # fallback: show entry exptl if available
            entry = st.session_state.entry_json
            if entry and entry.get("exptl"):
                arr = [[i, e.get("method","")] for i, e in enumerate(entry.get("exptl"))]
                df = pd.DataFrame(arr, columns=["Index","Method"]).astype(str)
                safe_table(df)
            else:
                st.info("No experimental details in the entry JSON.")
    else:
        st.info("Experimental details need a RCSB entry (fetch by PDB ID).")

# ---------- Page: Validation ----------
elif menu == "Validation":
    st.header("Validation metrics")
    st.markdown("Global validation metrics from RCSB (if available).")
    if st.session_state.pdb_source == "rcsb":
        val = fetch_validation(pdb_input)
        if val:
            g = val.get("global_validation_metrics", [])
            if g:
                rows = [[m.get("metric_name",""), str(m.get("value",""))] for m in g]
                safe_table(pd.DataFrame(rows, columns=["Metric","Value"]))
            else:
                st.info("No global validation metrics in the validation summary.")
        else:
            st.info("No validation summary available for this entry.")
    else:
        st.info("Validation summary requires a fetched RCSB entry.")

# ---------- Page: BLAST ----------
elif menu == "BLAST":
    st.header("BLASTp (protein vs PDB)")
    st.markdown("Select a chain sequence (Chains page) or paste a protein sequence and run BLASTp against the PDB database.")
    # sequence input area
    seq_area = st.text_area("Protein sequence (FASTA or raw). If left empty, selected chain sequence will be used.", value="")
    # allow user to pick selected chain from session
    if st.session_state.poly_seqs:
        default_choice = st.session_state.selected_entity or list(st.session_state.poly_seqs.keys())[0]
        chosen_chain = st.selectbox("Or choose chain sequence", options=["(use text area)"] + list(st.session_state.poly_seqs.keys()))
    else:
        chosen_chain = "(use text area)"
    blast_btn = st.button("Run BLASTp (top 5)")

    if blast_btn:
        # determine sequence
        seq = ""
        if seq_area.strip():
            # remove FASTA header if present
            lines = seq_area.strip().splitlines()
            if lines and lines[0].startswith(">"):
                seq = "".join([l.strip() for l in lines if not l.startswith(">")])
            else:
                seq = seq_area.strip()
        elif chosen_chain != "(use text area)":
            seq = st.session_state.poly_seqs.get(chosen_chain, "")
        else:
            seq = ""

        if not seq:
            st.error("No sequence provided. Paste a sequence or select a chain that has a sequence.")
        elif not email:
            st.error("Set your NCBI email in the sidebar before running BLASTp.")
        else:
            with st.spinner("Running remote BLASTp against PDB (NCBI qblast)..."):
                try:
                    handle = NCBIWWW.qblast("blastp", "pdb", seq, hitlist_size=5)
                    # parse results robustly
                    blast_iter = NCBIXML.parse(handle)
                    hits_all = []
                    for brec in blast_iter:
                        for aln in brec.alignments:
                            hsp = aln.hsps[0] if aln.hsps else None
                            hits_all.append({
                                "Title": str(aln.title),
                                "Length": str(aln.length),
                                "Score": str(hsp.score) if hsp else "",
                                "E-value": str(hsp.expect) if hsp else ""
                            })
                        break
                    try:
                        handle.close()
                    except Exception:
                        pass
                    if hits_all:
                        df_hits = pd.DataFrame(hits_all).astype(str)
                        safe_table(df_hits)
                    else:
                        st.info("No BLASTp hits returned.")
                except Exception as e:
                    st.error(f"BLASTp failed: {e}")

# ---------- Page: Downloads ----------
elif menu == "Downloads":
    st.header("Downloads")
    st.markdown("Download legacy PDB or mmCIF files directly from RCSB (if using a fetched PDB ID). You can also download selected chain FASTA.")
    if st.session_state.pdb_source == "rcsb" and pdb_input:
        pdb_url = f"https://files.rcsb.org/download/{pdb_input}.pdb"
        cif_url = f"https://files.rcsb.org/download/{pdb_input}.cif"
        st.markdown(f"- [Download PDB (.pdb)]({pdb_url})")
        st.markdown(f"- [Download mmCIF (.cif)]({cif_url})")
        # chain FASTA downloads handled on Chains page (download button)
    elif st.session_state.pdb_source == "upload":
        st.info("You uploaded a file; use Visualization menu to view. To download the same file, re-download from your local machine.")
    else:
        st.info("No PDB loaded (fetch or upload first).")

# ---------- Page: Coordinate Analysis ----------
elif menu == "Coordinate Analysis":
    st.header("Coordinate file analysis")
    st.markdown("Geometric and biophysical analysis of the loaded coordinate file (PDB/mmCIF). Requires a parsed structure (auto-parsed after fetch/upload if possible).")

    structure = st.session_state.get("parsed_structure", None)
    if structure is None:
        st.info("No parsed coordinate file available. Fetch a PDB entry (so the code can download coordinates) or upload a PDB/mmCIF file and click Fetch/Load.")
    else:
        df_atoms = analyze_structure_and_show(structure)
        if df_atoms is not None:
            st.success("Coordinate analysis complete.")

# ---------- Footer ----------
st.markdown("---")
st.caption("Data from RCSB PDB (data.rcsb.org). BLASTp uses NCBI remote BLAST (qblast). NGL viewer is used for 3D visualization.")