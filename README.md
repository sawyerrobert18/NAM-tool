# NAM Animal Impact Calculator

A web tool for wet-lab researchers to audit the **animal impact of their cell-culture reagent choices** — and to find certified animal-free alternatives. You pick an experiment and cell line, enter the amounts you actually use, and the tool back-calculates how many animals were required to produce those reagents, with **every figure traceable to a published source**.

Built in the spirit of the 3Rs and New Approach Methodologies (NAMs): make the hidden animal cost of routine reagents visible, quantifiable, and citable.

> **Author:** Muteeb Syed, Long island university, brooklyn.
> **Live app:** https://nam-tool-hmpkgigrqmsgagbjkrwavi.streamlit.app
> **Source:** https://github.com/sawyerrobert18/repo
> **Contact:** sawyerrobert18@gmail.com
> **License:** MIT

---

## What it does

Back-calculates animal impact of animal-derived reagents (fetal bovine serum, other sera, trypsin, BSA, collagens, EHS-derived matrices, etc.) from the amount you actually consume in an experiment.
Distinguishes killing from use. Reagents are categorised as *killed*, *killed (by-product)*, or *used not killed* (e.g. serum from live-donor herds via venipuncture), so the headline number is not inflated.
Cites every number. Each figure links to a primary source (product sheet, peer-reviewed paper, or vendor Certificate of Analysis). Estimates that are derived rather than directly published are flagged with an asterisk (`*`, Tier-2 — see below).
Surfaces animal-free alternatives for each reagent (e.g. human platelet lysate, recombinant laminins, synthetic matrices), with performance notes, price, and literature citations.
Exports a grant-ready PDF summarising inputs, per-reagent impact, costs, alternatives, and the full citation list — suitable for a Methods section or an ethics/funding application.

## How the numbers are derived

Every reagent carries a **data tier**:

Tier 1 — Direct. The animals-per-unit figure comes straight from a primary source (e.g. FBS: industry-average serum volume per fetal calf).
Tier 2 — Derived (`*`). No direct per-animal figure exists, so the value is reconstructed from *cited* anatomical/biochemical quantities multiplied by one clearly-stated efficiency term (e.g. bovine trypsin from pancreas mass × trypsin content × extraction yield). Tier-2 figures are always shown with an asterisk and a range.

The full derivation for every Tier-2 reagent — including sources and assumptions — is documented in [`back_calculation_methodology.md`](back_calculation_methodology.md). The guiding rule: *if a number can't be defended there, it doesn't appear in the tool.

## Running it locally

```bash
pip install -r requirements.txt
streamlit run app.py
```

Requires Python 3.9+. Dependencies: `streamlit`, `pandas`, `matplotlib`, `fpdf2`.

## Data

Reference data lives in `Database/` as CSVs:

| File | Contents |
|------|----------|
| `cell_lines_FINAL.csv` | Cell lines with base media, serum type/%, and animal-derived components (ATCC-verified where possible) |
| `reagents_FINAL.csv` | Animal-derived reagents with per-animal yields, tiers, and citations |
| `vendors_FINAL.csv` | Vendors with catalogue numbers and publicly-listed prices |
| `alternatives_FINAL.csv` | Animal-free alternatives with performance and citations |

## Limitations & roadmap

This is an honest research tool, not a certified accounting system. Known limitations:

Tier-2 estimates carry real uncertainty - the EHS/Matrigel figure spans roughly a 15× range; it is intended for order-of-magnitude comparison, not exact counts.
One protocol per run. Maintenance and passaging currently share a single event clock; separating them (media-change vs passage intervals) is planned.
Coverage is a curated subset of common cell lines, reagents, and experiment types - not exhaustive.
A small number of cell-line data points are literature-derived rather than direct-from-vendor and are flagged as lower-confidence in the app.

Planned: a checkbox-based reagent picker scoped per experiment type; separate maintenance/passaging protocols; expanded cell-line and reagent coverage; broader vendor price coverage.

## Feedback

Corrections to any cited figure are especially welcome - the tool is only as good as its sources. Please open an issue on the repository or email _[CONTACT]_.

---

*Built to make the animal cost of research reagents visible and, wherever possible, avoidable.*
