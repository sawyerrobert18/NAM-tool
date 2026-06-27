# Back-Calculation Methodology

**NAM Animal Impact Calculator** — how every animal-per-unit figure is derived.
Last updated: 2026-06-26.

This document is the audit trail behind each number the tool reports. It exists
so that a grant reviewer, ethics board, or skeptical colleague can trace any
figure back to a primary source or an explicitly-stated assumption. If a number
cannot be defended here, it should not appear in the tool.

---

## 1. Tiering system

Every reagent's `animals_per_unit` value carries a data tier:

**Tier 1 — Direct.** The animals-per-unit figure comes straight from a primary
source or product specification. No modelling. Example: fetal bovine serum, where
the industry-average serum volume per fetal calf is published.

**Tier 2 — Derived (marked with `*`).** No direct per-animal figure exists, so the
value is reconstructed from cited anatomical/biochemical quantities multiplied by
one or more efficiency terms. The cited quantities are real; the efficiency terms
are reasoned assumptions, stated explicitly. Every Tier-2 figure in the tool is
displayed with an asterisk and a "derived estimate" note.

The principle: a Tier-2 estimate is honest **only if its assumed terms are
disclosed**. A derived number presented as if it were measured is worse than no
number at all.

---

## 2. Tier-1 figures (direct)

**Fetal Bovine Serum (R001)** — 2.5 calves/L (range 2–3).
Industry average ~0.4 L serum per fetal calf. Pecora et al. 2020, *Translational
Animal Science*, DOI:10.1093/tas/txag044. FBS is **not** a slaughter by-product:
the dam is selected for pregnancy and the fetus is bled by cardiac puncture
post-slaughter. Counted as KILLED.

**Calf / Horse Serum (R002, R003)** — 0 animals killed.
Collected by venipuncture from living donor herds (ATCC 302030 spec; ATCC Animal
Cell Culture Guide). Animals are used but not killed. Reported as USED-NOT-KILLED†
with a welfare note, impact = 0.

**Base media, Pen-Strep, bacterial Collagenase II (R012–R017)** — 0 animals.
Fully synthetic or microbial/fermentation origin. No animal use.

---

## 3. Tier-2 figures (derived) — slaughterhouse by-products

These share a common structure:

```
mass of target protein per animal = organ/tissue mass
                                     x protein content fraction
                                     x extraction efficiency   [assumption]
animals per unit = (unit mass) / (mass per animal)
```

**Bovine Trypsin (R004)** — 0.78 animals/L of 0.25% solution (range 0.28–2.1).
Pancreas ~400 g × ~2% trypsin × ~80% protein fraction × ~50% extraction ≈
3.2 g/animal; 0.25% solution = 2.5 g/L → 2.5 / 3.2 = 0.78. Slaughterhouse
by-product (KILLED-BYPRODUCT). Worthington Biochemical documentation.

**Porcine Trypsin (R018)** — 0.6 animals/L (range 0.2–1.5).
Same method, porcine pancreas anatomy (~150–300 g). Thermo Fisher 25200056
(confirmed porcine).

**Bovine Serum Albumin (R007)** — 5.4 animals/kg (range 5.0–7.4).
30 L blood × 50% serum = 15 L; albumin 35 g/L × 15 L = 525 g; × 35% recovery ≈
184 g/animal → 1000 / 184 = 5.4. Slaughterhouse by-product. Sigma A7906/A1470.

**Bovine Collagen Type I (R006)** — 0.31 animals/kg (range 0.2–0.6).
Hide 27 kg × 30% collagen × 40% extraction ≈ 3240 g/animal → 1000 / 3240 = 0.31.
By-product. Advanced BioMatrix PureCol.

**Bovine Insulin (R008)** — 21 animals/g (range 11–44).
Pancreas 400 g × 200 µg/g insulin × 60% extraction ≈ 48 mg/animal →
1000 / 48 = 20.8. By-product. *Note: most culture insulin is now recombinant
human (0 animals) — verify before logging.*

**Bovine Fibronectin (R011)** — 0.44 animals/g (range 0.28–0.83).
15 L plasma × 300 µg/mL × 50% recovery ≈ 2.25 g/animal → 1 / 2.25 = 0.44.
By-product.

**Rat-Tail Collagen Type I (R019)** — ~10 animals/L (range 5–20).
Tail tendon ~0.5–1 g/rat × ~70% collagen × ~30–50% extraction; 3 mg/mL solution.
**Dedicated kill** — rats are killed specifically for tail tendon, *not* a
by-product, hence ~30–50× higher impact than bovine collagen. Corning 354249.

---

## 4. Tier-2 figures (derived) — EHS murine matrix (NEW MODEL)

Matrigel (R005) and murine Laminin-111 (R010) were previously marked
**NOT CALCULABLE** because no source publishes "mice per mL" directly. This was
*more* conservative than the rest of the tool: the by-product figures above all
rely on an undisclosed-in-source extraction efficiency too. By the tool's own
Tier-2 standard, EHS products are no less calculable than trypsin or BSA — they
simply stack one additional assumed term and therefore carry wider uncertainty.

The model below makes that explicit. It is implemented in `matrigel_model.py`
(runnable, every parameter adjustable).

### 4.1 Background

Matrigel is a basement-membrane extract (BME) from the Engelbreth-Holm-Swarm
(EHS) sarcoma, grown in C57BL/6 mice. Composition ≈ 60% laminin, 30% collagen IV,
8% nidogen/entactin, plus heparan-sulfate proteoglycans and bound growth factors.
Murine laminin-111 is the purified laminin fraction of the same tumour.

*According to PubMed:* canonical preparation protocols — Kleinman 2001,
*Current Protocols in Cell Biology*, [DOI:10.1002/0471143030.cb1002s00](https://doi.org/10.1002/0471143030.cb1002s00);
Kleinman et al. 1982, *Biochemistry*, [DOI:10.1021/bi00350a005](https://doi.org/10.1021/bi00350a005).
Composition figures: Kastana et al. 2019 (*Methods Mol Biol*); Aisenbrey et al.
2020, *Nature Reviews Materials* (synthetic alternatives review). Product protein
concentration: Corning Matrigel Certificate of Analysis (per-lot, 8–12 mg/mL).
Purified EHS laminin is sold at ~1–2 mg/mL (Sigma/Corning product specs).

### 4.2 The model

```
matrix protein per mouse (mg) = tumor mass (g) x 1000 x extractable fraction
mL product per mouse          = matrix protein per mouse / product conc (mg/mL)
                                [Laminin: x laminin fraction x purification recovery]
animals per mL                = 1 / (mL product per mouse)
```

| Parameter | Min | Mid | Max | Basis |
|-----------|-----|-----|-----|-------|
| EHS tumor wet mass per mouse (g) | 3 | 6 | 10 | Protocol descriptions (Kleinman) — **sourced range** |
| Extractable BM-protein fraction of wet mass | 2% | 4% | 6% | **Flagged assumption** — no open per-mouse mg figure |
| Matrigel protein conc (mg/mL) | 8 | 10 | 12 | Corning CoA — **sourced** |
| Laminin product conc (mg/mL) | 1 | 1 | 2 | Sigma/Corning laminin specs — **sourced** |
| Laminin fraction of matrix | — | 60% | — | Kastana 2019 — **sourced** |
| Laminin purification recovery | — | 50% | — | **Flagged assumption** |

`animals/mL` is **maximised** (worst case, most animals) by the *least* product
per mouse: small tumour × low extractable fraction × high product concentration.

### 4.3 Results

| Product | Min animals/mL | **Mid** | Max animals/mL |
|---------|---------------|---------|----------------|
| **Matrigel (R005)** | 0.013 | **0.042** | 0.20 |
| **Laminin-111 (R010)** | 0.0056 | **0.014** | 0.11 |

Mid-point trace (Matrigel): 6 g tumour × 4% = 240 mg matrix protein/mouse;
240 mg ÷ 10 mg/mL = 24 mL Matrigel/mouse; 1 ÷ 24 = 0.042 animals/mL.

Worked context (mid): coating one well with ~0.5 mL ≈ 0.02 mice; a 5 mL
experiment ≈ 0.21 mice.

### 4.4 Honesty statement

This is a **high-uncertainty Tier-2 estimate**. The min–max spread is ~15× for
Matrigel, driven mainly by the extractable-fraction assumption. It should be
reported with the asterisk, the range, and a note that the dominant uncertainty
is extraction efficiency. It is suitable for relative comparison and
order-of-magnitude awareness — **not** for a precise claim like "this experiment
killed exactly N mice." Where a project needs a hard number, the xeno-free
alternatives (Vitronectin XF, recombinant laminin-521, synthetic hydrogels)
sidestep the question entirely.

---

## 5. What is deliberately NOT modelled

- **EHS growth-factor content / batch variation** — out of scope; affects
  function, not animal count.
- **Mortality during tumour passaging** (mice that don't take the implant) —
  would raise the estimate; omitted for lack of a cited rate. The figures above
  are therefore conservative (lower bound on true impact).
- **Recombinant-equivalent products** — recombinant insulin, in-vitro mAbs,
  recombinant laminin-521 are 0 animals and must be logged as such, not via these
  bovine/murine figures.

---

*Attribution: literature figures retrieved via PubMed; DOIs linked inline above.
Product specifications from manufacturer Certificates of Analysis. Tool license: MIT.*
