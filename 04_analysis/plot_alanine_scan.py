"""
FoldX alanine scanning — CC12.1 / RBD interface
Generates 5 figures in 04_analysis/plots/

Run from the project root:
    python3 04_analysis/plot_alanine_scan.py
"""

import os
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import seaborn as sns
import numpy as np

BASE   = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA   = os.path.join(BASE, "03_data")
PLOTS  = os.path.join(BASE, "04_analysis", "plots")
os.makedirs(PLOTS, exist_ok=True)

# ── load data ──────────────────────────────────────────────────────────────────
ddg    = pd.read_csv(os.path.join(DATA, "binding_ddg.csv"))
bm     = pd.read_csv(os.path.join(DATA, "foldx_summary.csv"))
master = pd.read_csv(os.path.join(DATA, "analysecomplex_master.csv"))

# residue position for sequence plot
master["resnum"] = master["mutation"].str[1:-1].astype(int)

# hotspot tier colours
def tier_color(ddg_val):
    if ddg_val >= 2.0:   return "#d62728"   # hot spot
    if ddg_val >= 1.0:   return "#ff7f0e"   # warm spot
    if ddg_val >= 0.0:   return "#aec7e8"   # negligible
    return "#2ca02c"                         # anti-hotspot

colors = [tier_color(v) for v in ddg["ddG_binding"]]


# ── Figure 1: ΔΔG_binding bar chart ───────────────────────────────────────────
NOISE_SIGMA = 0.17   # stdev of WT interaction energies across paired runs
NOISE_BAND  = 2 * NOISE_SIGMA   # ±2σ = ±0.34 kcal/mol

fig, ax = plt.subplots(figsize=(10, 5))
# Noise band first so bars render on top
ax.axhspan(-NOISE_BAND, NOISE_BAND, color="grey", alpha=0.12, zorder=0, label=f"WT noise (±2σ = ±{NOISE_BAND:.2f})")
bars = ax.bar(ddg["mutation"], ddg["ddG_binding"], color=colors, edgecolor="black", linewidth=0.6, zorder=2)
ax.axhline( NOISE_BAND, color="grey", linestyle=":", linewidth=0.8, zorder=1)
ax.axhline(-NOISE_BAND, color="grey", linestyle=":", linewidth=0.8, zorder=1)
ax.axhline(2.0, color="#d62728", linestyle="--", linewidth=1.0)
ax.axhline(1.0, color="#ff7f0e", linestyle="--", linewidth=1.0)
ax.axhline(0.0, color="black",   linestyle="-",  linewidth=0.5)
ax.set_xlabel("Mutation", fontsize=12)
ax.set_ylabel("ΔΔG binding (kcal/mol)", fontsize=12)
ax.set_title("CC12.1–RBD Interface: Alanine Scanning ΔΔG (AnalyseComplex)", fontsize=13)
legend_patches = [
    mpatches.Patch(color="#d62728", label="Hot spot  (≥2.0)"),
    mpatches.Patch(color="#ff7f0e", label="Warm spot (1.0–2.0)"),
    mpatches.Patch(color="#aec7e8", label="Negligible (0–1.0)"),
    mpatches.Patch(color="#2ca02c", label="Anti-hotspot (<0)"),
    mpatches.Patch(color="grey",    alpha=0.3, label=f"WT noise ±2σ (±{NOISE_BAND:.2f} kcal/mol)"),
]
ax.legend(handles=legend_patches, fontsize=9, loc="upper left")
plt.xticks(rotation=45, ha="right")
plt.tight_layout()
plt.savefig(os.path.join(PLOTS, "01_ddG_binding_bar.png"), dpi=150)
plt.close()
print("Saved: 01_ddG_binding_bar.png")


# ── Figure 2: Stacked energy component bar chart ──────────────────────────────
# 5 mechanistic terms + "Other" catch-all so bars sum to true total.
# True total overlaid as black diamonds.
components = {
    "Sidechain H-bond":  "ddG_sidechain_hbond",
    "VdW":               "ddG_vdw",
    "Electrostatics":    "ddG_electrostatics",
    "Solvation polar":   "ddG_solvation_polar",
    "Solv. hydrophobic": "ddG_solvation_hydrophobic",
}
comp_colors = ["#4e79a7", "#f28e2b", "#e15759", "#76b7b2", "#59a14f"]

# "Other" = total minus the 5 named terms (captures VdW clashes, entropy, BB hbond, etc.)
named_sum = sum(master[col] for col in components.values())
master["ddG_other"] = master["ddG_interaction_energy"] - named_sum

fig, ax = plt.subplots(figsize=(12, 5))
x = np.arange(len(master))
bottoms_pos = np.zeros(len(master))
bottoms_neg = np.zeros(len(master))

all_components = list(components.items()) + [("Other terms", "ddG_other")]
all_colors     = comp_colors + ["#bdbdbd"]

for (label, col), color in zip(all_components, all_colors):
    vals = master[col].values
    pos_vals = np.where(vals > 0, vals, 0)
    neg_vals = np.where(vals < 0, vals, 0)
    ax.bar(x, pos_vals, bottom=bottoms_pos, label=label, color=color, edgecolor="white", linewidth=0.4)
    ax.bar(x, neg_vals, bottom=bottoms_neg, color=color, edgecolor="white", linewidth=0.4)
    bottoms_pos += pos_vals
    bottoms_neg += neg_vals

# Overlay true total as black diamonds
ax.scatter(x, master["ddG_interaction_energy"],
           color="black", marker="D", s=45, zorder=5, label="Total ΔΔG")

ax.axhline(0, color="black", linewidth=0.7)
ax.set_xticks(x)
ax.set_xticklabels(master["mutation"], rotation=45, ha="right")
ax.set_xlabel("Mutation", fontsize=12)
ax.set_ylabel("ΔΔG contribution (kcal/mol)", fontsize=12)
ax.set_title("Energy Component Breakdown per Alanine Substitution\n"
             "(bars sum to total; ◆ = true ΔΔG binding)", fontsize=12)
ax.legend(fontsize=9, bbox_to_anchor=(1.01, 1), loc="upper left")
plt.tight_layout()
plt.savefig(os.path.join(PLOTS, "02_energy_components_stacked.png"), dpi=150)
plt.close()
print("Saved: 02_energy_components_stacked.png")


# ── Figure 3: BuildModel vs AnalyseComplex scatter ────────────────────────────
merged = bm[["mutation", "ddG_total"]].merge(ddg[["mutation", "ddG_binding"]], on="mutation")

fig, ax = plt.subplots(figsize=(7, 6))
ax.scatter(merged["ddG_total"], merged["ddG_binding"],
           color=[tier_color(v) for v in merged["ddG_binding"]],
           s=80, edgecolors="black", linewidth=0.6, zorder=3)

for _, row in merged.iterrows():
    ax.annotate(row["mutation"],
                (row["ddG_total"], row["ddG_binding"]),
                textcoords="offset points", xytext=(6, 3), fontsize=8)

lim_min = min(merged["ddG_total"].min(), merged["ddG_binding"].min()) - 0.5
lim_max = max(merged["ddG_total"].max(), merged["ddG_binding"].max()) + 0.5
ax.plot([lim_min, lim_max], [lim_min, lim_max], "k--", linewidth=0.8, alpha=0.4, label="y = x")
ax.axhline(0, color="grey", linewidth=0.5)
ax.axvline(0, color="grey", linewidth=0.5)
ax.set_xlabel("ΔΔG fold stability — BuildModel (kcal/mol)", fontsize=12)
ax.set_ylabel("ΔΔG binding — AnalyseComplex (kcal/mol)", fontsize=12)
ax.set_title("Fold Stability vs Binding ΔΔG: These Are Different Quantities\n"
             "(points off the diagonal = fold effect ≠ binding effect)", fontsize=11)
ax.legend(fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(PLOTS, "03_buildmodel_vs_analysecomplex.png"), dpi=150)
plt.close()
print("Saved: 03_buildmodel_vs_analysecomplex.png")


# ── Figure 4: Heatmap (mutations × energy terms) ──────────────────────────────
heatmap_cols = {
    "Interact.\nenergy":   "ddG_interaction_energy",
    "SC\nH-bond":          "ddG_sidechain_hbond",
    "BB\nH-bond":          "ddG_backbone_hbond",
    "VdW":                 "ddG_vdw",
    "Electro-\nstatics":   "ddG_electrostatics",
    "Solv.\npolar":        "ddG_solvation_polar",
    "Solv.\nhydrophob":    "ddG_solvation_hydrophobic",
    "VdW\nclashes":        "ddG_vdw_clashes",
    "Entropy\nSC":         "ddG_entropy_sidechain",
}

heat_df = master.set_index("mutation")[[v for v in heatmap_cols.values()]]
heat_df.columns = list(heatmap_cols.keys())

vmax = max(abs(heat_df.values.min()), abs(heat_df.values.max()))
fig, ax = plt.subplots(figsize=(14, 6))
sns.heatmap(heat_df, annot=True, fmt=".2f", cmap="RdBu_r",
            center=0, vmin=-vmax, vmax=vmax,
            linewidths=0.5, linecolor="white",
            annot_kws={"size": 11},
            cbar_kws={"label": "ΔΔG (kcal/mol)", "shrink": 0.8},
            ax=ax)
ax.set_title("Heatmap: ΔΔG per Energy Term (AnalyseComplex)", fontsize=14)
ax.set_ylabel("")
ax.set_yticklabels(ax.get_yticklabels(), fontsize=11)
ax.set_xticklabels(ax.get_xticklabels(), rotation=0, fontsize=11)
plt.tight_layout()
plt.savefig(os.path.join(PLOTS, "04_energy_heatmap.png"), dpi=150)
plt.close()
print("Saved: 04_energy_heatmap.png")


# ── Figure 5: Sequence lollipop ───────────────────────────────────────────────
# Manual x-offsets (points) for labels that crowd in the 495–502 cluster
label_xoffsets = {
    "Y495A": -18, "G496A": -8, "Q498A": 4, "N501A": 14, "G502A": 24,
}

fig, ax = plt.subplots(figsize=(12, 6))
for _, row in master.iterrows():
    color = tier_color(row["ddG_interaction_energy"])
    ax.vlines(row["resnum"], 0, row["ddG_interaction_energy"],
              color=color, linewidth=1.5)
    ax.scatter(row["resnum"], row["ddG_interaction_energy"],
               color=color, s=70, edgecolors="black", linewidth=0.6, zorder=3)
    xoff = label_xoffsets.get(row["mutation"], 0)
    yoff = 10 if row["ddG_interaction_energy"] >= 0 else -16
    ax.annotate(row["mutation"],
                (row["resnum"], row["ddG_interaction_energy"]),
                textcoords="offset points",
                xytext=(xoff, yoff),
                ha="center", fontsize=9, rotation=0)

ax.axhspan(-NOISE_BAND, NOISE_BAND, color="grey", alpha=0.12, zorder=0)
ax.axhline( NOISE_BAND, color="grey", linestyle=":", linewidth=0.8, zorder=1)
ax.axhline(-NOISE_BAND, color="grey", linestyle=":", linewidth=0.8, zorder=1)
ax.axhline(2.0, color="#d62728", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axhline(1.0, color="#ff7f0e", linestyle="--", linewidth=0.8, alpha=0.7)
ax.axhline(0.0, color="black",   linestyle="-",  linewidth=0.5)

# Mark A475 — interface residue that cannot be scanned (already Ala in WT)
ax.scatter(475, 0, marker="^", color="none", edgecolors="#888888",
           s=80, linewidth=1.2, zorder=4)
ax.annotate("A475\n(unscanned)", (475, 0),
            textcoords="offset points", xytext=(0, 10),
            ha="center", fontsize=8, color="#666666", style="italic")

ax.set_xlabel("RBD Residue Position", fontsize=12)
ax.set_ylabel("ΔΔG binding (kcal/mol)", fontsize=12)
ax.set_title("CC12.1–RBD Hot Spots by Sequence Position", fontsize=13)
ax.set_xlim(408, 510)
legend_patches = [
    mpatches.Patch(color="#d62728", label="Hot spot  (≥2.0)"),
    mpatches.Patch(color="#ff7f0e", label="Warm spot (1.0–2.0)"),
    mpatches.Patch(color="#aec7e8", label="Negligible"),
    mpatches.Patch(color="#2ca02c", label="Anti-hotspot"),
    mpatches.Patch(color="grey",    alpha=0.3, label=f"WT noise ±2σ (±{NOISE_BAND:.2f} kcal/mol)"),
]
ax.legend(handles=legend_patches, fontsize=9)
plt.tight_layout()
plt.savefig(os.path.join(PLOTS, "05_sequence_lollipop.png"), dpi=150)
plt.close()
print("Saved: 05_sequence_lollipop.png")

print("\nAll plots saved to 04_analysis/plots/")
