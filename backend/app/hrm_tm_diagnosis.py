import math
from app.engines.hrm_engine import _amplicon_tm

# SMOKING GUN 1: hold GC% fixed at 50%, vary LENGTH. Real amplicon Tm should keep
# rising with length (more base pairs = more stable). Two-state NN should asymptote.
unit = "ATGC"  # 50% GC repeating
print("Fixed 50% GC, varying length — does two-state NN Tm respond to length?")
print(f"{'length':>8}{'NN Tm':>10}")
for L in [20, 40, 60, 80, 120, 200, 400, 1000]:
    s = (unit * (L//4 + 1))[:L]
    print(f"{L:>8}{_amplicon_tm(s):>10.2f}")
print()
print("-> The two-state NN Tm flatlines: the concentration term R*ln(C/4) becomes")
print("   negligible vs the length-scaled dS, so Tm -> dH/dS ratio = f(GC) ONLY.")
print()

# SMOKING GUN 2: the GC-Tm slope compresses at the top of the human range.
print("Fixed length 90 bp, varying GC% — Tm sensitivity per 5% GC step:")
print(f"{'GC%':>6}{'NN Tm':>10}{'dTm/+5%GC':>12}")
prev=None
for gcp in [30,35,40,45,50,55,60,65,70]:
    n_gc = round(90*gcp/100); n_at = 90-n_gc
    # interleave to avoid runs
    s = ("GC"*45 + "AT"*45)
    s = "".join(("G" if i< n_gc else "A") for i in range(90))
    # better: deterministic GC fraction, alternating to keep NN sane
    seq=[]
    g=0
    for i in range(90):
        want_gc = (g+0.5)/(i+1) < gcp/100
        if want_gc and g< n_gc:
            seq.append("G" if i%2 else "C"); g+=1
        else:
            seq.append("A" if i%2 else "T")
    s="".join(seq)
    t=_amplicon_tm(s)
    slope = "" if prev is None else f"{t-prev:>12.2f}"
    print(f"{gcp:>6}{t:>10.2f}{slope}")
    prev=t
