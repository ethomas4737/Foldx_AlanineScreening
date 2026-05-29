# CC12.1 Antibody–RBD FoldX Alanine Scanning

Computational alanine scanning of the CC12.1 antibody / SARS-CoV-2 RBD interface using FoldX 5.1.

---

## Project Structure

```
project_root/
├── 01_structures/
│   ├── original/          CC12.1_RBD_clean.pdb       — starting structure (pre-FoldX)
│   ├── repaired/          CC12_Repair.pdb             — RepairPDB output (use for all runs)
│   ├── mutants/           CC12_Repair_1–12.pdb        — 12 Ala mutant structures (BuildModel)
│   └── wt_paired/         WT_CC12_Repair_1–12.pdb     — paired WT structures (BuildModel)
│
├── 02_foldx_results/
│   ├── buildmodel/
│   │   ├── individual_list.txt    — 12 mutations in FoldX notation (e.g. NA487A;)
│   │   ├── Dif_CC12_Repair.fxout  — ΔΔG fold stability per mutation
│   │   ├── Average_CC12_Repair.fxout
│   │   └── Raw_CC12_Repair.fxout
│   └── analysecomplex/
│       ├── wt_analysis/           — AnalyseComplex on repaired WT complex
│       ├── mutant_analysis/       — AnalyseComplex on 12 Ala mutants
│       └── wt_paired_analysis/    — AnalyseComplex on 12 paired WT structures
│
├── 03_data/
│   ├── binding_ddg.csv            — ΔΔG_binding per mutation (key results table)
│   ├── foldx_summary.csv          — BuildModel ΔΔG with energy component breakdown
│   └── analysecomplex_master.csv  — Full AnalyseComplex decomposition (78 columns)
│
├── 04_analysis/
│   ├── plot_alanine_scan.py       — Generates all 5 figures
│   ├── plots/                     — Output figures (PNG)
│   │   ├── 01_ddG_binding_bar.png
│   │   ├── 02_energy_components_stacked.png
│   │   ├── 03_buildmodel_vs_analysecomplex.png
│   │   ├── 04_energy_heatmap.png
│   │   └── 05_sequence_lollipop.png
│   ├── notebooks/
│   └── reports/
│       └── alanine_scan_report.md — Full interpretation of results
│
├── foldx_20261231                 — FoldX 5.1 binary (not included; see Dependencies)
└── molecules/                     — FoldX rotamer library
```

---

## Complex Description

| Chain | Identity | Residues |
|-------|----------|----------|
| A | SARS-CoV-2 RBD | 334–528 |
| H | CC12.1 antibody heavy chain | 1–216 |
| L | CC12.1 antibody light chain | 1–215 |
| Z | RBD (second copy) | 334–528 |
| X | Heavy chain (second copy) | 1–216 |
| Y | Light chain (second copy) | 1–215 |

Chains Z/X/Y are a symmetric copy in the crystal structure. All mutations were made on chain A; AnalyseComplex analysed the **A vs H+L** interface only.

---

## Workflow

### Step 1 — RepairPDB
```bash
./foldx_20261231 --command=RepairPDB --pdb=CC12.1_RBD_clean.pdb
```
Output: `CC12_Repair.pdb` → move to `01_structures/repaired/`

### Step 2 — BuildModel (alanine scanning)
```bash
./foldx_20261231 --command=BuildModel --pdb=CC12_Repair.pdb \
    --mutant-file=individual_list.txt
```
Output: `CC12_Repair_1–12.pdb` (mutants) and `WT_CC12_Repair_1–12.pdb` (paired WTs)

### Step 3 — AnalyseComplex on WT
```bash
./foldx_20261231 --command=AnalyseComplex --pdb=CC12_Repair.pdb \
    --analyseComplexChains=A,HL
```

### Step 4 — AnalyseComplex on all mutants and paired WTs
```bash
for i in $(seq 1 12); do
    ./foldx_20261231 --command=AnalyseComplex \
        --pdb=CC12_Repair_${i}.pdb --analyseComplexChains=A,HL
    ./foldx_20261231 --command=AnalyseComplex \
        --pdb=WT_CC12_Repair_${i}.pdb --analyseComplexChains=A,HL
done
```

### Step 5 — Parse results and plot
```bash
pip3 install matplotlib seaborn pandas
python3 04_analysis/plot_alanine_scan.py
```

---

## Key Results

ΔΔG_binding = G_bind(mutant) − G_bind(paired WT). Positive = residue contributes to binding.

| Mutation | ΔΔG binding | Classification |
|----------|-------------|----------------|
| N487A | +2.57 kcal/mol | **Hot spot** |
| N501A | +1.84 kcal/mol | Warm spot |
| K417A | +1.56 kcal/mol | Warm spot |
| L455A | +1.30 kcal/mol | Warm spot |
| G476A | +0.50 kcal/mol | Borderline |
| G416A | +0.39 kcal/mol | Within noise |
| S459A, G496A, Q498A | < ±0.34 | Within noise (n=1) |
| Y495A, T415A, G502A | negative | Anti-hotspot |

Full interpretation in `04_analysis/reports/alanine_scan_report.md`.

---

## Important Caveats

- **n = 1 replicate per mutation.** Results with |ΔΔG| < 0.34 kcal/mol (±2σ WT noise) should not be interpreted.
- **BuildModel ΔΔG ≠ binding ΔΔG.** Do not use `Dif_CC12_Repair.fxout` to rank binding hot spots — use `binding_ddg.csv` (from AnalyseComplex) instead.
- **Gly→Ala mutations** (G416A, G476A, G496A, G502A) report backbone packing effects, not sidechain contacts.
- **A475 cannot be scanned** — it is already Ala in the WT sequence but is an interface residue.
- Only the **RBD side** of the interface was scanned. Antibody hot spots are not characterised here.

---

## Dependencies

### FoldX

FoldX is free for academic use but must be downloaded directly from the developers — it cannot be redistributed here.

1. Register and download FoldX 5.1 from [foldxsuite.crg.eu](https://foldxsuite.crg.eu/)
2. Place the binary in the project root and rename it `foldx_20261231` (or adjust the commands in the workflow above to match your binary name)
3. Make it executable: `chmod +x foldx_20261231`

The `molecules/` rotamer library included in this repo is distributed alongside FoldX and must be in the same directory as the binary when running commands.

### Python

- Python ≥ 3.9
- matplotlib, seaborn, pandas (`pip3 install matplotlib seaborn pandas`)
