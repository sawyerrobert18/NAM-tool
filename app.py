"""
NAM Animal Impact Calculator  v1.1
Run with:  streamlit run app.py
License:   MIT
"""

import re
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd
import streamlit as st
from fpdf import FPDF

# ─── Page config ──────────────────────────────────────────────────────────────
st.set_page_config(
    page_title="NAM Animal Impact Calculator",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ─── Data directory ────────────────────────────────────────────────────────────
_BASE = Path(__file__).parent
DATA_DIR = next(
    (p for name in ("data", "Database") if (p := _BASE / name).is_dir()), None
)
if DATA_DIR is None:
    st.error("Cannot find data directory (data/ or Database/) next to app.py.")
    st.stop()

# ─── Vessel volumes ────────────────────────────────────────────────────────────
VESSEL_FILL_ML = {
    "T25 flask": 6.0,
    "T75 flask": 17.5,
    "T175 flask": 45.0,
    "6-well plate (per well)": 2.0,
    "12-well plate (per well)": 1.0,
    "24-well plate (per well)": 0.5,
    "96-well plate (per well)": 0.2,
}
TRYPSIN_ML = {
    "T25 flask": 2.0,
    "T75 flask": 5.0,
    "T175 flask": 10.0,
    "6-well plate (per well)": 0.5,
    "12-well plate (per well)": 0.3,
    "24-well plate (per well)": 0.2,
    "96-well plate (per well)": 0.1,
}

# ─── Experiment types ─────────────────────────────────────────────────────────
EXPERIMENT_TYPES = {
    "Cell line maintenance": {
        "uses_cell_line": True,
        "uses_vessels": True,
        "extra_reagents": {},
        "desc": "Routine subculture and passaging of adherent or suspension cell lines.",
    },
    "Cell-based assay (cytotoxicity / proliferation / viability)": {
        "uses_cell_line": True,
        "uses_vessels": True,
        "extra_reagents": {},
        "desc": "Functional assay run on cultured cells (MTT, LDH, BrdU, etc.). Uses standard culture reagents.",
    },
    "ELISA": {
        "uses_cell_line": False,
        "uses_vessels": False,
        "extra_reagents": {"R007": ("BSA used for blocking buffer (g)", 0.30)},
        "desc": "Enzyme-linked immunosorbent assay. Main animal-derived reagent is BSA in blocking buffer.",
    },
    "Western Blot": {
        "uses_cell_line": False,
        "uses_vessels": False,
        "extra_reagents": {"R007": ("BSA used for blocking (g)", 0.50)},
        "desc": "Protein detection by SDS-PAGE and antibody probing. BSA used for membrane blocking.",
    },
    "Immunofluorescence / IHC": {
        "uses_cell_line": False,
        "uses_vessels": False,
        "extra_reagents": {"R007": ("BSA used for blocking buffer (g)", 0.10)},
        "desc": "Antibody-based imaging. BSA used in blocking and antibody dilution buffer.",
    },
    "Flow Cytometry": {
        "uses_cell_line": False,
        "uses_vessels": False,
        "extra_reagents": {"R007": ("BSA used in staining buffer (g)", 0.20)},
        "desc": "Cell surface/intracellular marker detection. BSA used in FACS buffer.",
    },
    "Organoid / 3D culture": {
        "uses_cell_line": True,
        "uses_vessels": True,
        "extra_reagents": {},
        "desc": "3D organoid culture requiring Matrigel or equivalent matrix. Matrigel impact is NOT CALCULABLE.",
    },
}

NOT_CALCULABLE_RIDS = {"R005", "R010"}
NO_ANIMAL_RIDS     = {"R012", "R013", "R014", "R015", "R016", "R017"}

VESSEL_CITATION = (
    "ATCC Complete Guide to Cell Culture; Corning Cell Culture Guide (standard fill volumes)."
)
TRYPSIN_VOL_CITATION = "Corning Cell Culture Guide -- standard trypsin dissociation volumes."

# ─── Data loading ─────────────────────────────────────────────────────────────
@st.cache_data
def load_data():
    cell_lines = pd.read_csv(DATA_DIR / "cell_lines_FINAL.csv", dtype=str).fillna("")
    reagents   = pd.read_csv(DATA_DIR / "reagents_FINAL.csv",   dtype=str).fillna("")
    alts       = pd.read_csv(DATA_DIR / "alternatives_FINAL.csv", dtype=str).fillna("")
    vendors    = pd.read_csv(DATA_DIR / "vendors_FINAL.csv",    dtype=str).fillna("")
    return cell_lines, reagents, alts, vendors

# ─── Helpers ──────────────────────────────────────────────────────────────────
def sanitize_pdf(s: str) -> str:
    """Replace unicode chars that latin-1 cannot encode."""
    table = {
        "—": "--", "–": "-", "’": "'", "‘": "'",
        "“": '"',  "”": '"', "†": "(dagger)", "°": "deg",
        "µ": "u",  "×": "x", "≥": ">=", "≤": "<=",
        "®": "(R)", "™": "(TM)", "±": "+/-",
    }
    for ch, rep in table.items():
        s = s.replace(ch, rep)
    return s.encode("latin-1", errors="replace").decode("latin-1")


def parse_components(comp_str: str) -> list:
    s = comp_str.lower()
    ids = []
    def _add(r):
        if r not in ids: ids.append(r)
    if "fbs" in s or "fetal bovine serum" in s:            _add("R001")
    if ("calf bovine serum" in s or "calf serum" in s) and "fetal" not in s: _add("R002")
    if "horse serum" in s or "equine serum" in s:          _add("R003")
    if "bovine trypsin" in s:                              _add("R004")
    if "porcine trypsin" in s:                             _add("R018")
    if "matrigel" in s:                                    _add("R005")
    if "laminin" in s and "recombinant" not in s:          _add("R010")
    if ("bovine collagen" in s or "purecol" in s) and "rat" not in s: _add("R006")
    if ("bsa" in s or "bovine serum albumin" in s) and "recombinant" not in s: _add("R007")
    if ("bovine insulin" in s or "insulin bovine" in s) and "recombinant" not in s: _add("R008")
    if "murine" in s and ("antibod" in s or "mab" in s):  _add("R009")
    if "fibronectin" in s and "recombinant" not in s:      _add("R011")
    if "rat tail collagen" in s:                           _add("R019")
    if "collagenase" in s:                                 _add("R017")
    return ids


def parse_serum_pct(pct_str: str, cl_id: str) -> dict:
    s = str(pct_str).strip()
    if not s or s == "0": return {}
    if "HS" in s and "FBS" in s:
        result = {}
        m = re.search(r"(\d+(?:\.\d+)?)\s*HS", s)
        if m: result["R003"] = float(m.group(1))
        m = re.search(r"(\d+(?:\.\d+)?)\s*FBS", s)
        if m: result["R001"] = float(m.group(1))
        return result
    try:
        pct = float(re.search(r"\d+(?:\.\d+)?", s).group())
        return {"R002": pct} if cl_id == "CL009" else {"R001": pct}
    except: return {}


def parse_cost(raw: str):
    if not raw or any(k in str(raw).lower() for k in ("login", "require", "n/a")):
        return None
    m = re.search(r"\d[\d,]*(?:\.\d+)?", str(raw).replace(",", ""))
    return float(m.group().replace(",", "")) if m else None


def pack_to_L(pack_str: str):
    ps = pack_str.strip().lower()
    m = re.match(r"^(\d+(?:\.\d+)?)\s*ml$", ps)
    if m: return float(m.group(1)) / 1000.0
    m = re.match(r"^(\d+(?:\.\d+)?)\s*l$", ps)
    if m: return float(m.group(1))
    return None


def vendor_cpl(vrow):
    cost = parse_cost(vrow["cost_per_unit_usd"])
    if cost is None: return None
    vol  = pack_to_L(vrow["pack_size"])
    if not vol: return None
    return cost / vol


def alt_cpl(arow):
    cost = parse_cost(arow["cost_per_unit_usd"])
    if cost is None: return None
    if arow["unit"].strip().lower() == "liter": return cost
    vol = pack_to_L(arow["unit"])
    return cost / vol if (vol and vol > 0) else None


def fmt_n(n: float) -> str:
    if n < 0.001:  return f"{n:.5f}"
    if n < 0.01:   return f"{n:.4f}"
    if n < 1:      return f"{n:.3f}"
    if n < 10:     return f"{n:.2f}"
    return f"{n:.1f}"


def perf_label(p: str) -> str:
    pl = p.strip().lower()
    if pl == "drop-in":        return "Drop-in replacement"
    if "optimization" in pl:   return "Optimization required"
    return p


def cat_label(c: str) -> str:
    cu = c.upper()
    if "USED-NOT-KILLED" in cu: return "Used, not killed"
    if "BYPRODUCT"        in cu: return "Killed (byproduct)"
    if "KILLED"           in cu: return "Killed"
    if "NO-ANIMAL"        in cu: return "No animal"
    return c


# ─── Core calculation ─────────────────────────────────────────────────────────
def calculate_impact(
    cl_row, vessel_type, n_vessels, duration_days, change_freq_days,
    serum_override, vendor_sels, extra_quantities,
    reagents_df, vendors_df, alts_df, exp_type_key,
):
    results   = []
    warnings  = []
    citations = []

    exp_type = EXPERIMENT_TYPES[exp_type_key]
    uses_vessels   = exp_type["uses_vessels"]
    uses_cell_line = exp_type["uses_cell_line"]

    # ── Media volumes ─────────────────────────────────────────────────────────
    total_media_L = 0.0
    tryp_L        = 0.0
    serum_dict    = {}

    if uses_cell_line and uses_vessels:
        fill_mL       = VESSEL_FILL_ML[vessel_type]
        n_changes     = duration_days / change_freq_days
        total_media_L = fill_mL * n_vessels * n_changes / 1000.0
        citations.append({
            "what": f"Vessel fill volumes ({vessel_type}: {fill_mL} mL)",
            "source": VESSEL_CITATION,
            "url": "https://www.atcc.org/resources/culture-guides/animal-cell-culture-guide",
        })

        cl_id = cl_row["cell_line_id"]
        serum_dict = (
            {"R001": float(serum_override)} if serum_override is not None
            else parse_serum_pct(cl_row["serum_pct"], cl_id)
        )

        if cl_row["trypsin_required"].lower().startswith("yes"):
            n_passages = duration_days / change_freq_days
            tryp_mL    = TRYPSIN_ML[vessel_type] * n_vessels * n_passages
            tryp_L     = tryp_mL / 1000.0
            citations.append({
                "what": f"Trypsin dissociation ({vessel_type}: {TRYPSIN_ML[vessel_type]} mL per vessel/passage)",
                "source": TRYPSIN_VOL_CITATION,
                "url": "",
            })

    # ── Which reagent IDs to process ──────────────────────────────────────────
    if uses_cell_line:
        reagent_ids = parse_components(cl_row["animal_derived_components"])
    else:
        reagent_ids = []

    # Add IDs from extra_quantities (ELISA/WB/Flow user inputs)
    for rid in extra_quantities:
        if rid not in reagent_ids:
            reagent_ids.append(rid)

    # ── Per-reagent loop ───────────────────────────────────────────────────────
    for rid in reagent_ids:
        if rid in NO_ANIMAL_RIDS:
            continue

        rdf = reagents_df[reagents_df["reagent_id"] == rid]
        if rdf.empty: continue
        r = rdf.iloc[0]

        mid_raw   = r["animals_per_unit_mid"].strip()
        is_tier2  = r["estimate_flag"].strip() == "*"
        these_alts = alts_df[alts_df["replaces_reagent_id"].str.contains(rid, regex=False)]

        # NOT CALCULABLE
        if rid in NOT_CALCULABLE_RIDS or mid_raw in ("N/A", ""):
            warnings.append({"type": "not_calculable", "reagent_id": rid, "name": r["name"]})
            citations.append({
                "what": f"{r['name']} -- sourcing data (NOT CALCULABLE)",
                "source": r["source_citation"], "url": r["source_url"],
            })
            results.append({
                "reagent_id": rid, "name": r["name"],
                "not_calculable": True, "category": r["animal_use_category"],
                "animals_killed": r["animals_killed"], "is_tier2": is_tier2,
                "alternatives": these_alts,
                "source_citation": r["source_citation"], "source_url": r["source_url"],
                "amount_display": "N/A", "amount_unit": "",
            })
            continue

        mid = float(mid_raw)
        unit = r["unit"]

        # ── Volume / weight for this reagent ──────────────────────────────────
        if rid in extra_quantities:
            # Weight-based input (ELISA/WB/Flow): quantity in grams -> kg
            qty_g  = extra_quantities[rid]
            qty_kg = qty_g / 1000.0
            # BSA unit is kg; animals_per_unit_mid = 5.4 per kg
            vol_L         = qty_kg
            amount_display = f"{qty_g:.2f} g"
            amount_unit    = "kg"
        elif rid in ("R001", "R002", "R003"):
            pct   = serum_dict.get(rid, 0.0)
            vol_L = total_media_L * pct / 100.0
            amount_display = f"{vol_L * 1000:.1f} mL"
            amount_unit    = "liter"
        elif rid in ("R004", "R018"):
            if tryp_L <= 0: continue
            vol_L          = tryp_L
            amount_display = f"{tryp_L * 1000:.1f} mL"
            amount_unit    = "liter of 0.25% solution"
        else:
            warnings.append({"type": "needs_input", "reagent_id": rid, "name": r["name"]})
            continue

        if vol_L <= 0: continue

        # Normalise vol_L to the reagent's native unit
        # (serum: liter, trypsin: liter 0.25%, BSA: kg)
        animals_total = vol_L * mid

        # ── Vendor ────────────────────────────────────────────────────────────
        v_opts  = vendors_df[vendors_df["reagent_id"] == rid]
        sel_vid = vendor_sels.get(rid, "")
        if sel_vid and not v_opts[v_opts["vendor_id"] == sel_vid].empty:
            vrow = v_opts[v_opts["vendor_id"] == sel_vid].iloc[0]
        elif not v_opts.empty:
            vrow = v_opts.iloc[0]
        else:
            vrow = None

        cpl  = vendor_cpl(vrow) if vrow is not None else None
        # For weight-based: cpl is per-liter but we have kg. Need cost per kg.
        if rid == "R007" and vrow is not None:
            # BSA sold by gram/kg; pack_size might be "100g"
            pack = vrow["pack_size"].strip().lower()
            m_g  = re.match(r"^(\d+)\s*g$", pack)
            m_kg = re.match(r"^(\d+(?:\.\d+)?)\s*kg$", pack)
            cost_raw = parse_cost(vrow["cost_per_unit_usd"])
            if cost_raw is not None:
                if m_g:
                    cpl = cost_raw / (float(m_g.group(1)) / 1000.0)  # per kg
                elif m_kg:
                    cpl = cost_raw / float(m_kg.group(1))
                else:
                    cpl = None
        cost = cpl * vol_L if cpl is not None else None

        citations.append({
            "what": f"{r['name']} -- {mid} animals per {unit}",
            "source": r["source_citation"], "url": r["source_url"],
        })

        results.append({
            "reagent_id": rid, "name": r["name"],
            "not_calculable": False,
            "vol_L": vol_L,
            "amount_display": amount_display,
            "amount_unit": amount_unit,
            "animals_per_unit": mid,
            "unit": unit,
            "animals_total": animals_total,
            "category": r["animal_use_category"],
            "animals_killed": r["animals_killed"],
            "is_tier2": is_tier2,
            "vendor": vrow, "cost_per_L": cpl, "current_cost": cost,
            "alternatives": these_alts,
            "source_citation": r["source_citation"], "source_url": r["source_url"],
        })

    return results, warnings, citations


# ─── Matplotlib chart ─────────────────────────────────────────────────────────
def make_bar_chart(results):
    calc = [r for r in results if not r.get("not_calculable") and r.get("animals_total", 0) > 0]
    if not calc:
        return None

    names  = [r["name"] for r in calc]
    values = [r["animals_total"] for r in calc]
    colors = []
    for r in calc:
        ak = r["animals_killed"].upper()
        if "NOT-KILLED" in ak:
            colors.append("#f59e0b")    # amber
        else:
            colors.append("#ef4444")    # red

    fig, ax = plt.subplots(figsize=(7, max(1.5, len(calc) * 0.7)))
    fig.patch.set_facecolor("#0e1117")
    ax.set_facecolor("#0e1117")

    bars = ax.barh(names, values, color=colors, height=0.5)
    ax.set_xlabel("Animals", color="#9ca3af", fontsize=9)
    ax.tick_params(colors="#9ca3af", labelsize=9)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    for spine in ("bottom", "left"):
        ax.spines[spine].set_color("#374151")

    # Value labels
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_width() + max(values) * 0.01,
            bar.get_y() + bar.get_height() / 2,
            fmt_n(val),
            va="center", ha="left", color="#9ca3af", fontsize=8,
        )

    ax.set_xlim(0, max(values) * 1.25)
    plt.tight_layout(pad=0.5)
    return fig


# ─── PDF export ───────────────────────────────────────────────────────────────
def build_pdf(cl_name, exp_type_key, inputs, results, warnings, citations):
    pdf = FPDF()
    pdf.add_page()
    pdf.set_margins(20, 20, 20)
    pdf.set_auto_page_break(auto=True, margin=20)

    def ln(text, bold=False, sz=10, indent=0):
        pdf.set_font("Helvetica", "B" if bold else "", sz)
        if indent:
            pdf.set_x(20 + indent)
        pdf.multi_cell(0, 5, sanitize_pdf(text))

    ln("NAM Animal Impact Calculator -- Results Report", bold=True, sz=15)
    ln("Every number cites a source. Open source, MIT license.", sz=8)
    pdf.ln(3)

    ln("Experiment Inputs", bold=True, sz=12)
    for k, v in inputs.items():
        ln(f"  {k}: {v}", sz=9)
    pdf.ln(3)

    # NOT CALCULABLE warnings
    not_calc = [w for w in warnings if w["type"] == "not_calculable"]
    if not_calc:
        pdf.set_text_color(180, 0, 0)
        ln("WARNING -- NOT CALCULABLE", bold=True, sz=11)
        pdf.set_text_color(0, 0, 0)
        for w in not_calc:
            ln(
                f"  {w['name']}: Animal impact cannot be reliably calculated. "
                "The supplier does not publish sourcing yield data. "
                "Do not cite in a publication.", sz=8
            )
        pdf.ln(2)

    # Totals
    calc = [r for r in results if not r.get("not_calculable")]
    total_killed = sum(
        r["animals_total"] for r in calc
        if "NOT-KILLED" not in r["animals_killed"].upper()
    )
    total_unp = sum(
        r["animals_total"] for r in calc
        if "NOT-KILLED" in r["animals_killed"].upper()
    )
    ln("Animal Impact Results", bold=True, sz=12)
    ln(f"  Total animals killed: {fmt_n(total_killed)}", bold=True, sz=10)
    ln(f"  Total animals used (not killed): {fmt_n(total_unp)}", bold=True, sz=10)
    pdf.ln(2)

    # Per-reagent table
    ln("Reagent breakdown:", bold=True, sz=10)
    for r in results:
        if r["not_calculable"]:
            ln(f"  {r['name']}: NOT CALCULABLE -- see warning", sz=8)
        else:
            flag = " *" if r["is_tier2"] else ""
            unp  = "NOT-KILLED" in r["animals_killed"].upper()
            cost_str = f"  |  Cost: ${r['current_cost']:.2f}" if r.get("current_cost") else ""
            ln(
                f"  {r['name']}{flag}: {r['amount_display']} used  |  "
                + ("0 (used not killed)" if unp else f"{fmt_n(r['animals_total'])} animals killed")
                + cost_str,
                sz=8
            )
    pdf.ln(3)

    # Alternatives
    ln("Available Animal-Free Alternatives", bold=True, sz=12)
    for r in results:
        if r["alternatives"].empty: continue
        ln(f"  Alternatives for {r['name']}:", bold=True, sz=9)
        for _, alt in r["alternatives"].iterrows():
            ln(
                f"    - {alt['name']}  |  {perf_label(alt['performance_equivalence'])}  "
                f"|  ${alt['cost_per_unit_usd']}/{alt['unit']}",
                sz=8
            )
    pdf.ln(3)

    # Citations
    ln("Citations", bold=True, sz=12)
    seen = set()
    idx  = 1
    for c in citations:
        key = c["source"][:100]
        if key in seen: continue
        seen.add(key)
        ln(f"[{idx}] {c['what']}", sz=7)
        ln(f"    Source: {c['source']}", sz=7, indent=4)
        if c.get("url"):
            ln(f"    {c['url']}", sz=7, indent=4)
        pdf.ln(1)
        idx += 1

    pdf.set_font("Helvetica", "I", 7)
    pdf.ln(2)
    pdf.multi_cell(0, 4, sanitize_pdf(
        "* Tier 2 estimate -- derived from biological data; see back_calculation_methodology.md\n"
        "(dagger) Animals used but not killed (venipuncture); counted as 0 in total killed.\n"
        "Tool: NAM Animal Impact Calculator | License: MIT"
    ))
    return bytes(pdf.output())


# ─── Streamlit UI ─────────────────────────────────────────────────────────────
def main():
    cell_lines, reagents, alts, vendors = load_data()

    # ── Sidebar ───────────────────────────────────────────────────────────────
    st.sidebar.header("Experiment Design")

    # 1. Experiment type first
    exp_type_key = st.sidebar.selectbox(
        "Experiment type", list(EXPERIMENT_TYPES.keys())
    )
    exp_type = EXPERIMENT_TYPES[exp_type_key]
    st.sidebar.caption(exp_type["desc"])
    st.sidebar.markdown("---")

    # 2. Cell line (only if experiment uses cells)
    cl_row      = None
    component_ids = []
    vendor_sels   = {}

    if exp_type["uses_cell_line"]:
        cl_options = {
            f"{row['common_name']} ({row['atcc_catalog']})": row["cell_line_id"]
            for _, row in cell_lines.iterrows()
        }
        cl_label = st.sidebar.selectbox("Cell line", list(cl_options.keys()))
        cl_id    = cl_options[cl_label]
        cl_row   = cell_lines[cell_lines["cell_line_id"] == cl_id].iloc[0]
        if cl_id == "CL013":
            st.sidebar.warning(
                "U2OS: LOW-CONFIDENCE -- confirmed via literature, not ATCC product page directly."
            )

    # 3. Vessel inputs
    vessel_type, n_vessels, duration, freq = "T75 flask", 6, 14, 3.0
    serum_override = None

    if exp_type["uses_vessels"]:
        vessel_type = st.sidebar.selectbox(
            "Vessel type", list(VESSEL_FILL_ML.keys()), index=1
        )
        n_vessels = int(st.sidebar.number_input("Number of vessels", 1, value=6, step=1))
        duration  = int(st.sidebar.number_input("Experiment duration (days)", 1, value=14, step=1))
        freq      = st.sidebar.number_input(
            "Media change / passage every N days", 0.5, value=3.0, step=0.5
        )

        st.sidebar.markdown("---")
        default_pct = cl_row["serum_pct"] if cl_row is not None else "10"
        if st.sidebar.checkbox(f"Override serum % (default: {default_pct})"):
            serum_override = float(st.sidebar.slider("FBS %", 0, 20, 10))

    # 4. Extra reagent quantities (ELISA, WB, etc.)
    extra_quantities = {}
    if exp_type["extra_reagents"]:
        st.sidebar.markdown("---")
        st.sidebar.subheader("Reagent quantities")
        for rid, (label, default) in exp_type["extra_reagents"].items():
            qty = st.sidebar.number_input(label, min_value=0.0, value=default, step=0.05)
            if qty > 0:
                extra_quantities[rid] = qty

    # 5. Vendor dropdowns
    if exp_type["uses_cell_line"] and cl_row is not None:
        component_ids = parse_components(cl_row["animal_derived_components"])
    for rid in extra_quantities:
        if rid not in component_ids:
            component_ids.append(rid)

    if component_ids:
        st.sidebar.markdown("---")
        st.sidebar.subheader("Vendor selection")
        for rid in component_ids:
            if rid in NO_ANIMAL_RIDS or rid in NOT_CALCULABLE_RIDS:
                continue
            v_opts = vendors[vendors["reagent_id"] == rid]
            if v_opts.empty:
                continue
            r_name = reagents[reagents["reagent_id"] == rid]["name"].values
            label  = r_name[0] if len(r_name) else rid
            v_labels = {
                f"{vr['vendor_name']} {vr['catalog_number']} ({vr['pack_size']}, ${vr['cost_per_unit_usd']})": vr["vendor_id"]
                for _, vr in v_opts.iterrows()
            }
            chosen = st.sidebar.selectbox(f"Vendor -- {label}", list(v_labels.keys()))
            vendor_sels[rid] = v_labels[chosen]

    # ── Calculate ─────────────────────────────────────────────────────────────
    dummy_cl = cell_lines.iloc[0] if cl_row is None else cl_row
    results, warnings, citations = calculate_impact(
        dummy_cl, vessel_type, n_vessels, duration, float(freq),
        serum_override, vendor_sels, extra_quantities,
        reagents, vendors, alts, exp_type_key,
    )

    # ── Derived totals ────────────────────────────────────────────────────────
    calc_results = [r for r in results if not r.get("not_calculable")]
    total_killed = sum(
        r["animals_total"] for r in calc_results
        if "NOT-KILLED" not in r["animals_killed"].upper()
    )
    total_unk    = sum(
        r["animals_total"] for r in calc_results
        if "NOT-KILLED" in r["animals_killed"].upper()
    )
    all_costs    = [r["current_cost"] for r in calc_results if r["current_cost"] is not None]
    total_cost   = sum(all_costs) if all_costs else None

    # ── Page header ───────────────────────────────────────────────────────────
    st.title("NAM Animal Impact Calculator")
    st.caption(
        f"Experiment type: **{exp_type_key}** "
        + (f"| Cell line: **{cl_row['common_name']}**" if cl_row is not None else "")
        + "  |  Every number cites a source.  |  * = Tier 2 estimate  |  (dagger) = used, not killed"
    )

    # NOT CALCULABLE alerts
    for w in [w for w in warnings if w["type"] == "not_calculable"]:
        st.error(
            f"**{w['name']} -- NOT CALCULABLE.** "
            "The supplier does not publish per-animal yield data. "
            "Do not cite this figure in a publication or grant."
        )

    # ── Tabs ─────────────────────────────────────────────────────────────────
    tab1, tab2, tab3 = st.tabs(["Animal Impact", "Alternatives", "Citations"])

    # ════════ TAB 1 ════════════════════════════════════════════════════════════
    with tab1:
        # Top metrics
        m1, m2, m3 = st.columns(3)
        m1.metric("Animals killed", fmt_n(total_killed))
        m2.metric("Animals used, not killed (dagger)", fmt_n(total_unk))
        m3.metric(
            "Current reagent cost",
            f"${total_cost:.2f}" if total_cost is not None else "N/A",
        )

        if not calc_results and not [w for w in warnings if w["type"] == "not_calculable"]:
            st.info("No calculable animal-derived reagents for this experiment configuration.")
        else:
            st.markdown("---")
            # Chart
            fig = make_bar_chart(results)
            if fig:
                col_chart, col_space = st.columns([2, 1])
                with col_chart:
                    st.pyplot(fig, use_container_width=True)
                plt.close(fig)

            st.markdown("---")
            # Quantities table
            table_rows = []
            for r in results:
                unp  = "NOT-KILLED" in r["animals_killed"].upper()
                row = {
                    "Reagent": r["name"] + (" *" if r.get("is_tier2") else ""),
                    "Amount used": r.get("amount_display", "N/A"),
                    "Category": cat_label(r["category"]),
                    "Animals/unit": (
                        "N/A" if r.get("not_calculable")
                        else f"{r['animals_per_unit']} per {r['unit']}"
                    ),
                    "Animals killed": (
                        "NOT CALCULABLE" if r.get("not_calculable")
                        else ("0 (dagger)" if unp else fmt_n(r["animals_total"]))
                    ),
                    "Cost": (
                        f"${r['current_cost']:.2f}" if r.get("current_cost") else "N/A"
                    ),
                }
                table_rows.append(row)

            if table_rows:
                df_display = pd.DataFrame(table_rows)
                st.dataframe(
                    df_display,
                    use_container_width=True,
                    hide_index=True,
                )

            # Source footnotes
            st.markdown("---")
            for r in calc_results:
                if r.get("is_tier2"):
                    st.caption(
                        f"* **{r['name']}**: Tier 2 estimate -- derived from published biology. "
                        "See back_calculation_methodology.md."
                    )
                src = r["source_citation"]
                url = r["source_url"]
                st.caption(
                    f"**{r['name']}** source: [{src}]({url})" if url else f"**{r['name']}** source: {src}"
                )

        # "Needs input" notices
        ni = [w for w in warnings if w["type"] == "needs_input"]
        if ni:
            st.markdown("---")
            for w in ni:
                st.info(
                    f"**{w['name']}**: used in this protocol but requires a quantity "
                    "input (g or mL). Weight-based calculation is a planned feature."
                )

        st.markdown("---")
        st.caption(
            "* Tier 2: derived from organ mass x protein content x extraction yield.  "
            "(dagger) Animals used but not killed (venipuncture from live donor herd).  "
            "All figures use mid-point estimates. Min/max ranges in database."
        )

        # PDF export
        inputs_summary = {
            "Experiment type": exp_type_key,
        }
        if cl_row is not None:
            inputs_summary["Cell line"] = f"{cl_row['common_name']} ({cl_row['atcc_catalog']})"
        if exp_type["uses_vessels"]:
            inputs_summary.update({
                "Vessel": f"{vessel_type} x {n_vessels}",
                "Duration": f"{duration} days",
                "Change frequency": f"every {freq} days",
            })
        if extra_quantities:
            for rid, qty in extra_quantities.items():
                r_name = reagents[reagents["reagent_id"] == rid]["name"].values
                inputs_summary[r_name[0] if len(r_name) else rid] = f"{qty} g"

        try:
            pdf_bytes = build_pdf(
                cl_row["common_name"] if cl_row is not None else exp_type_key,
                exp_type_key, inputs_summary, results, warnings, citations,
            )
            cl_name_safe = (
                cl_row["common_name"].replace(" ", "_") if cl_row is not None else exp_type_key.split()[0]
            )
            st.download_button(
                "Export PDF (grant-ready)",
                data=pdf_bytes,
                file_name=f"NAM_impact_{cl_name_safe}.pdf",
                mime="application/pdf",
            )
        except Exception as e:
            st.warning(f"PDF export error: {e}")

    # ════════ TAB 2: ALTERNATIVES ══════════════════════════════════════════════
    with tab2:
        # Total savings metric at top
        # (same reagents, same numbers as Tab 1 -- because results list is shared)
        st.metric(
            "Total animals that could be saved by switching all reagents to alternatives",
            fmt_n(total_killed + total_unk),
        )
        st.caption(
            "This matches the sum shown in Animal Impact tab. "
            "Cost delta shown per-reagent below."
        )
        st.markdown("---")

        any_alts = False
        for r in results:
            if r["alternatives"].empty:
                continue
            any_alts = True
            unp = "NOT-KILLED" in r["animals_killed"].upper()

            with st.expander(f"{r['name']}", expanded=True):
                hdr1, hdr2 = st.columns(2)
                if not r.get("not_calculable"):
                    hdr1.metric(
                        "Current: animals killed",
                        "0 (dagger)" if unp else fmt_n(r["animals_total"]),
                    )
                    hdr2.metric(
                        "Current: cost",
                        f"${r['current_cost']:.2f}" if r.get("current_cost") else "N/A",
                    )
                else:
                    st.info("NOT CALCULABLE -- see warning at top of page.")

                alt_rows = []
                for _, alt in r["alternatives"].iterrows():
                    acpl     = alt_cpl(alt)
                    validated = cl_row is not None and cl_row["common_name"] in alt["validated_cell_lines"]
                    cost_str  = f"${alt['cost_per_unit_usd']}/{alt['unit']}"
                    if r.get("cost_per_L") and acpl:
                        delta     = acpl - r["cost_per_L"]
                        sign      = "+" if delta >= 0 else ""
                        delta_str = f"{sign}${delta:.0f}/L vs current"
                    else:
                        delta_str = "N/A"
                    saved_str = (
                        fmt_n(r["animals_total"])
                        if not r.get("not_calculable") and not unp else "N/A"
                    )
                    alt_rows.append({
                        "Alternative": alt["name"],
                        "Performance": perf_label(alt["performance_equivalence"]),
                        "Validated for this line": "Yes" if validated else "No",
                        "Price": cost_str,
                        "Cost delta": delta_str,
                        "Animals saved": saved_str,
                        "Animal impact": "0",
                    })

                if alt_rows:
                    st.dataframe(
                        pd.DataFrame(alt_rows),
                        use_container_width=True,
                        hide_index=True,
                    )
                if alt["notes"]:
                    st.caption(alt["notes"][:300])

        if not any_alts:
            st.info("No alternatives found in database for the reagents in this experiment.")

        if any(r["reagent_id"] in NOT_CALCULABLE_RIDS for r in results):
            st.warning(
                "Matrigel/Laminin note: GrowDex and VitroGel did NOT support vascular organoid "
                "survival (Scientific Reports 2025, DOI:10.1038/s41598-025-20091-w). "
                "PEG hydrogels work for iPSC-CM but not all organoid types. "
                "Verify per-application."
            )

    # ════════ TAB 3: CITATIONS ════════════════════════════════════════════════
    with tab3:
        st.markdown(
            "Every source used in this calculation. "
            "Copy directly into a grant Methods or supplementary table."
        )
        if cl_row is not None:
            st.markdown(f"**Cell line protocol:** {cl_row['source_citation']}")
            if cl_row["source_url"]:
                st.markdown(cl_row["source_url"])
            st.markdown("---")

        seen = set()
        idx  = 1
        for c in citations:
            key = c["source"][:120]
            if key in seen: continue
            seen.add(key)
            st.markdown(f"**[{idx}]** {c['what']}")
            st.markdown(f"&nbsp;&nbsp;&nbsp;&nbsp;*{c['source']}*")
            if c.get("url"):
                st.markdown(f"&nbsp;&nbsp;&nbsp;&nbsp;{c['url']}")
            st.markdown("")
            idx += 1

        alt_seen = set()
        for r in results:
            for _, alt in r["alternatives"].iterrows():
                key = alt["source_citation"][:80]
                if key in alt_seen or not alt["source_citation"]: continue
                alt_seen.add(key)
                st.markdown(f"**[{idx}]** Alternative: {alt['name']}")
                st.markdown(f"&nbsp;&nbsp;&nbsp;&nbsp;*{alt['source_citation']}*")
                if alt["source_url"]:
                    st.markdown(f"&nbsp;&nbsp;&nbsp;&nbsp;{alt['source_url']}")
                st.markdown("")
                idx += 1

        st.markdown("---")
        st.caption(
            "* Tier 2 derived estimates are documented in full in "
            "back_calculation_methodology.md (included in this repository)."
        )


if __name__ == "__main__":
    main()
