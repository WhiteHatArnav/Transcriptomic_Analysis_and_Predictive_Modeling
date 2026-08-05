import re
import pandas as pd
from collections import Counter, OrderedDict

# =========================
# Configuration
# =========================
input_expression_csv = "/Users/arnavjoshi/Desktop/LUADTranscAnalysis/luad_expression_matrix_with_race_sampletype_and_subtypes.csv"

# =========================
# Helpers
# =========================
def _str(s):
    return "" if pd.isna(s) else str(s)

def _simplify(text):
    txt = _str(text).lower()
    txt = re.sub(r"[_\-\/]+", " ", txt)
    txt = re.sub(r"[^\w\s\(\):\.]", " ", txt)
    txt = re.sub(r"\s+", " ", txt).strip()
    return txt

def normalize_race(value):
    t = _simplify(value)
    if "white" in t:
        return "white"
    if "black" in t and "african" in t:
        return "black or african american"
    if t in {"black", "african american", "black or african-american"}:
        return "black or african american"
    return t  # preserves other categories if present

def normalize_sample_type(value):
    t = _simplify(value)
    return "normal" if "normal" in t else t

def parse_integrative(value):
    """
    Returns '1'..'6' when a clear iCluster digit is present.
    Accepts forms like '6', 'iCluster:6', 'iCluster 6'.
    Returns '' for missing/unclear.
    """
    t = _simplify(value)
    m = re.search(r"icluster[:\s\-]*([1-6])\b", t)
    if m:
        return m.group(1)
    # exact digit token
    m = re.search(r"\b([1-6])\b", t)
    if m:
        return m.group(1)
    return ""

def parse_immune_code(value):
    """
    Returns 'C1'..'C6' when a clear immune subtype code is present.
    Accepts forms like '... (Immune C3)', 'C2', or names mapped below.
    """
    t = _simplify(value)
    m = re.search(r"\bc([1-6])\b", t)
    if m:
        return f"C{m.group(1)}"
    # name-based fallbacks
    if "wound healing" in t:
        return "C1"
    if "ifn" in t or "gamma" in t:
        return "C2"
    if "inflammatory" in t:
        return "C3"
    if "lymphocyte depleted" in t:
        return "C4"
    if "immunologically quiet" in t:
        return "C5"
    if "tgf" in t:
        return "C6"
    return ""

# =========================
# Load matrix (6-level column header)
# Levels: [Patient ID, Sample ID, Race, Sample Type, Integrative Subtype, Immune Subtype]
# =========================
expr_df = pd.read_csv(input_expression_csv, header=[0,1,2,3,4,5], index_col=0)

# Build a tidy view of column annotations
cols = expr_df.columns
meta = pd.DataFrame({
    "patient_id": [c[0] for c in cols],
    "sample_id": [c[1] for c in cols],
    "race_raw":   [c[2] for c in cols],
    "stype_raw":  [c[3] for c in cols],
    "integr_raw": [c[4] for c in cols],
    "immune_raw": [c[5] for c in cols],
})

# Normalize annotations
meta["race"] = meta["race_raw"].map(normalize_race)
meta["sample_type"] = meta["stype_raw"].map(normalize_sample_type)
meta["integr"] = meta["integr_raw"].map(parse_integrative)       # '1'..'6' or ''
meta["immune"] = meta["immune_raw"].map(parse_immune_code)       # 'C1'..'C6' or ''

# =========================
# Counting logic
# =========================
race_white = meta["race"] == "white"
race_black = meta["race"] == "black or african american"

# 1) Race = white & Integrative subtype == "6"
white_integr6 = meta[race_white & (meta["integr"] == "6")].shape[0]

# 2) Race = black or african american & Integrative subtype == "6"
black_integr6 = meta[race_black & (meta["integr"] == "6")].shape[0]

# 3) Race = white & each immune subtype separately
immune_labels = [f"C{i}" for i in range(1,7)]
white_immune_counts = OrderedDict(
    (lab, meta[race_white & (meta["immune"] == lab)].shape[0]) for lab in immune_labels
)

# 4) Race = black or african american & each immune subtype separately
black_immune_counts = OrderedDict(
    (lab, meta[race_black & (meta["immune"] == lab)].shape[0]) for lab in immune_labels
)

# 5) Race = white & "Normal" samples (sample_type contains 'normal')
white_normal = meta[race_white & (meta["sample_type"] == "normal")].shape[0]

# 6) Race = black or african american & "Normal" samples
black_normal = meta[race_black & (meta["sample_type"] == "normal")].shape[0]

# =========================
# Report
# =========================
print("Counts for LUAD (case-insensitive, normalized matching)\n")

print("Integrative Subtype == 6")
print(f"  White: {white_integr6}")
print(f"  Black or African American: {black_integr6}\n")

print("Immune Subtype by Race (C1..C6)")
print("  White:")
for lab, cnt in white_immune_counts.items():
    print(f"    {lab}: {cnt}")
print("  Black or African American:")
for lab, cnt in black_immune_counts.items():
    print(f"    {lab}: {cnt}")
print()

print('Sample Type == "Normal"')
print(f"  White: {white_normal}")
print(f"  Black or African American: {black_normal}")
