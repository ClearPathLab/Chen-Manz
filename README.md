# Chen-Manz Structure Checker

A batch checker for MOF crystal structures (CIF files) that flags misbonded or isolated atoms. It is intended for data cleaning before high-throughput computational screening.

## Acknowledgement

This code is developed from the structure-checking implementation in the **CoRE MOF 2024 DB**:

> G. Zhao, L. M. Brabson, S. Chheda, *et al.*, Y. G. Chung. CoRE MOF DB: A curated experimental metal-organic framework database with machine-learned properties for integrated material-process screening. *Matter* **8** (6), 2025. https://doi.org/10.1016/j.matt.2025.102140

The checking method was originally proposed by Chen and Manz:

> T. Chen, T. A. Manz. Identifying misbonded atoms in the 2019 CoRE metal–organic framework database. *RSC Advances* **10** (45), 26944–26951, 2020. https://doi.org/10.1039/d0ra02498h

If you use this code, please cite both papers (BibTeX at the end of this file).

## Checks

Let *r<sub>a</sub>* and *r<sub>b</sub>* be the covalent radii of atoms *a* and *b*, and let *d* be their distance under periodic boundary conditions (minimum image convention). Two atoms are considered bonded if *d* ≤ *r<sub>a</sub>* + *r<sub>b</sub>*.

| Column | Meaning | Criterion |
|---|---|---|
| `overlapping` | Overlapping atoms | Any atom pair with *d* < 0.5 (*r<sub>a</sub>* + *r<sub>b</sub>*) |
| `isolated` | Isolated atom | Any atom with no bonded neighbors |
| `under_carbon` | Under-bonded carbon (often caused by missing H) | Sum of bond orders of a carbon < 3.3 |
| `over_carbon` | Over-bonded carbon | Sum of bond orders of a carbon ≥ 5.5 |
| `misplaced_hydro` | Misplaced hydrogen | An H atom that has both an N/O neighbor and a metal neighbor within *r<sub>a</sub>* + *r<sub>b</sub>* + 0.3 Å |

The carbon bond-order check only applies to a carbon atom whose neighbors are all among H, B, C, N, O, Cl and Br. Carbon atoms with any other neighbor are skipped.

The bond order of each bond is computed as BO = 10<sup>A·d + C</sup>, where A and C are element-specific parameters. The C–H bond order is capped at 1.25.

Elements missing from the covalent radius table default to 0.8 Å. Metalloids are not included in the metal list.

## Requirements

```bash
pip install ase numpy tqdm
```

## Usage

```bash
python chen_manz.py <cif_dir> <output_csv> [-j WORKERS]
```

| Argument | Description |
|---|---|
| `cif_dir` | Directory containing the CIF files; every `*.cif` file in it is processed |
| `output_csv` | Path of the output CSV file |
| `-j, --workers` | Number of parallel processes (default: 32) |

Example:

```bash
python chen_manz.py ../../dataset/cifs ../../dataset/5checkers/Chen-Manz.csv -j 32
```

## Output

The output is a CSV file with one row per structure:

```
name,overlapping,under_carbon,over_carbon,isolated,misplaced_hydro
ABAVIJ_clean,False,False,False,False,False
...
```

- `name` is the CIF file name without its extension.
- Each check is reported as `True` or `False`. If a structure cannot be read or processed, all of its checks are reported as `unknown`.
- Runs are resumable. If `output_csv` already exists, structures already listed in it are skipped, and new results are appended to the end of the file.

## Notes

- In an earlier version, the `misplaced_hydro` check compared atom indices against element symbols, so it always returned `False`. This has been fixed. If you have output from the old version, delete the old CSV and rerun the script.
- The distance matrix has size N × N, so cells with a very large number of atoms can use a lot of memory.
