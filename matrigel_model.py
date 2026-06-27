"""
EHS / Matrigel / Laminin animal-impact back-calculation model
==============================================================

Goal: estimate animals (mice) consumed per mL of product, so the NAM
calculator can show a Tier-2 (asterisked) number instead of "NOT CALCULABLE".

Method mirrors the app's existing Tier-2 derivations (trypsin, BSA, collagen):
  cited anatomical/product figures  x  one explicitly-flagged efficiency term.

The chain:
  matrix protein recoverable per mouse (mg)
      = tumor mass per mouse (g) x 1000 (mg/g)
        x extractable BM-protein fraction of wet tumor  [FLAGGED ASSUMPTION]
  mL product per mouse
      = matrix protein per mouse (mg) / product protein concentration (mg/mL)
        x product-specific recovery factor              [Laminin only]
  animals per mL
      = 1 / (mL product per mouse)

Every parameter is given as (min, mid, max). Worst case for animals/mL
(the MOST animals) = least product per mouse = small tumor x low fraction
x high concentration. Best case = the opposite.

Sources (see back_calculation_methodology.md for full citations):
  - EHS tumor grown in C57BL/6 mice; Matrigel = BME extract.
    Kleinman 2001 Curr Protoc Cell Biol  DOI:10.1002/0471143030.cb1002s00
    Kleinman et al. 1982 Biochemistry     DOI:10.1021/bi00350a005
  - Composition ~60% laminin / 30% collagen IV / 8% nidogen.
    Kastana et al. 2019 Methods Mol Biol; Aisenbrey et al. 2020 Nat Rev Mater.
  - Matrigel protein concentration 8-12 mg/mL: Corning Matrigel
    Certificate of Analysis (per-lot spec).
  - Purified EHS laminin sold at ~1 mg/mL: Sigma/Corning laminin product specs.
"""

from dataclasses import dataclass


@dataclass
class Tri:
    """A triangular (min, mid, max) parameter."""
    lo: float
    mid: float
    hi: float

    def __iter__(self):
        return iter((self.lo, self.mid, self.hi))


# ---------------------------------------------------------------- parameters

# Tumor wet mass harvested per mouse (grams).
# EHS tumors are grown 2-4 weeks then harvested at several grams each.
# Bounded by protocol descriptions (Kleinman); treated as a sourced range.
TUMOR_MASS_G = Tri(3.0, 6.0, 10.0)

# Fraction of wet tumor mass recoverable as extractable basement-membrane
# protein. THIS IS THE FLAGGED ASSUMPTION (no open per-mouse mg figure).
# EHS is an unusually matrix-rich tumor; soft-tissue protein is ~10-20% of
# wet mass and a large share of EHS protein is extractable BM. We take a
# deliberately conservative 2-6% of wet mass as *recovered* Matrigel protein.
EXTRACTABLE_FRACTION = Tri(0.02, 0.04, 0.06)

# Matrigel product protein concentration (mg/mL). Corning CoA range.
MATRIGEL_CONC_MG_ML = Tri(8.0, 10.0, 12.0)

# Purified laminin-111 product concentration (mg/mL). Sold dilute.
LAMININ_CONC_MG_ML = Tri(1.0, 1.0, 2.0)

# Laminin as fraction of matrix protein (~60%) x purification recovery (~50%).
LAMININ_MASS_FRACTION = 0.60
LAMININ_PURIFICATION_RECOVERY = 0.50


# ---------------------------------------------------------------- core model

def protein_per_mouse_mg(tumor_g, frac):
    return tumor_g * 1000.0 * frac


def matrigel_animals_per_ml(tumor_g, frac, conc):
    protein = protein_per_mouse_mg(tumor_g, frac)
    ml_per_mouse = protein / conc
    return 1.0 / ml_per_mouse, ml_per_mouse


def laminin_animals_per_ml(tumor_g, frac, conc):
    protein = protein_per_mouse_mg(tumor_g, frac)
    laminin_mg = protein * LAMININ_MASS_FRACTION * LAMININ_PURIFICATION_RECOVERY
    ml_per_mouse = laminin_mg / conc
    return 1.0 / ml_per_mouse, ml_per_mouse


def three_point(product):
    """Return (min, mid, max) animals/mL.
    animals/mL is maximized by the LEAST product per mouse:
      small tumor x low fraction x high concentration (for laminin: low conc->more mL->fewer animals,
      so laminin 'most animals' uses high conc too)."""
    if product == "matrigel":
        f = matrigel_animals_per_ml
        # most animals: small tumor, low frac, high conc
        hi_animals, _ = f(TUMOR_MASS_G.lo, EXTRACTABLE_FRACTION.lo, MATRIGEL_CONC_MG_ML.hi)
        mid_animals, mid_ml = f(TUMOR_MASS_G.mid, EXTRACTABLE_FRACTION.mid, MATRIGEL_CONC_MG_ML.mid)
        lo_animals, _ = f(TUMOR_MASS_G.hi, EXTRACTABLE_FRACTION.hi, MATRIGEL_CONC_MG_ML.lo)
    else:
        f = laminin_animals_per_ml
        hi_animals, _ = f(TUMOR_MASS_G.lo, EXTRACTABLE_FRACTION.lo, LAMININ_CONC_MG_ML.hi)
        mid_animals, mid_ml = f(TUMOR_MASS_G.mid, EXTRACTABLE_FRACTION.mid, LAMININ_CONC_MG_ML.mid)
        lo_animals, _ = f(TUMOR_MASS_G.hi, EXTRACTABLE_FRACTION.hi, LAMININ_CONC_MG_ML.lo)
    return lo_animals, mid_animals, hi_animals, mid_ml


def fmt(x):
    # Kept identical to app.py fmt_n so model and app displays agree.
    if x == 0:
        return "0"
    if x < 0.001:
        return f"{x:.5f}"
    if x < 0.01:
        return f"{x:.4f}"
    if x < 1:
        return f"{x:.3f}"
    if x < 10:
        return f"{x:.2f}"
    return f"{x:.1f}"


if __name__ == "__main__":
    print("=" * 64)
    print("EHS Matrigel / Laminin animal-impact model")
    print("=" * 64)

    for product, label, conc in (
        ("matrigel", "Matrigel (R005)", MATRIGEL_CONC_MG_ML),
        ("laminin", "Laminin-111 (R010)", LAMININ_CONC_MG_ML),
    ):
        lo, mid, hi, mid_ml = three_point(product)
        print(f"\n{label}")
        print(f"  product conc (mg/mL): {conc.lo}-{conc.hi} (mid {conc.mid})")
        print(f"  mid: ~{mid_ml:.0f} mL product per mouse")
        print(f"  animals per mL:  min {fmt(lo)}  |  MID {fmt(mid)}  |  max {fmt(hi)}")
        # worked context examples
        for use_ml, ctx in ((0.5, "0.5 mL (one coated well)"),
                            (5.0, "5 mL (typical experiment)")):
            print(f"    {ctx}: {fmt(mid * use_ml)} mice (mid)")

    print("\n" + "-" * 64)
    print("Mid-point parameter trace (Matrigel):")
    p = protein_per_mouse_mg(TUMOR_MASS_G.mid, EXTRACTABLE_FRACTION.mid)
    print(f"  tumor {TUMOR_MASS_G.mid} g x {EXTRACTABLE_FRACTION.mid:.0%} "
          f"= {p:.0f} mg matrix protein / mouse")
    print(f"  {p:.0f} mg / {MATRIGEL_CONC_MG_ML.mid} mg/mL "
          f"= {p/MATRIGEL_CONC_MG_ML.mid:.0f} mL Matrigel / mouse")
    print(f"  1 / {p/MATRIGEL_CONC_MG_ML.mid:.0f} mL "
          f"= {fmt(MATRIGEL_CONC_MG_ML.mid/p)} animals / mL")
