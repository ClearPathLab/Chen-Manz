import argparse
import csv
import os
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

import numpy as np
from ase.io import read
from tqdm import tqdm

ATR = {
    'H': 0.38, 'Li': 0.86, 'Be': 0.53, 'B': 1.01, 'C': 0.88, 'N': 0.86, 'O': 0.89, 'F': 0.82,
    'Na': 1.15, 'Mg': 1.28, 'Al': 1.53, 'Si': 1.38, 'P': 1.28, 'S': 1.20, 'Cl': 1.17, 'K': 1.44,
    'Ca': 1.17, 'Sc': 1.62, 'Ti': 1.65, 'V': 1.51, 'Cr': 1.53, 'Mn': 1.53, 'Fe': 1.43, 'Co': 1.31,
    'Ni': 1.33, 'Cu': 1.31, 'Zn': 1.41, 'Ga': 1.40, 'Ge': 1.35, 'As': 1.39, 'Se': 1.40, 'Br': 1.39,
    'Rb': 1.65, 'Sr': 1.30, 'Y': 1.84, 'Zr': 1.73, 'Nb': 1.66, 'Mo': 1.57, 'Ru': 1.58, 'Rh': 1.63,
    'Pd': 1.68, 'Ag': 1.56, 'Cd': 1.56, 'In': 1.53, 'Sn': 1.64, 'Sb': 1.64, 'Te': 1.65, 'I': 1.58,
    'Cs': 1.85, 'Ba': 1.52, 'La': 1.91, 'Ce': 1.98, 'Pr': 1.75, 'Nd': 1.92, 'Sm': 1.89, 'Eu': 1.83,
    'Gd': 1.79, 'Tb': 1.82, 'Dy': 1.79, 'Ho': 1.63, 'Er': 1.80, 'Tm': 1.84, 'Yb': 1.80, 'Lu': 1.86,
    'Hf': 1.73, 'W': 1.33, 'Re': 1.29, 'Ir': 1.50, 'Pt': 1.66, 'Au': 1.68, 'Hg': 1.88, 'Pb': 1.72,
    'Bi': 1.72, 'Th': 1.97, 'U': 1.76, 'Np': 1.73, 'Pu': 1.71,
}
DEFAULT_R = 0.8

COEF = {
    'H': (-0.6093, 0.5927), 'B': (-2.2011, 3.4380), 'C': (-1.2685, 1.8855),
    'N': (-1.2680, 1.8401), 'O': (-1.0525, 1.5189), 'Cl': (-0.7621, 1.3723),
    'Br': (-0.8003, 1.5272),
}

METALS = {
    'Li', 'Be', 'Na', 'Mg', 'Al', 'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu',
    'Zn', 'Ga', 'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Cs',
    'Ba', 'La', 'Ce', 'Pr', 'Nd', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu', 'Hf',
    'W', 'Re', 'Ir', 'Pt', 'Au', 'Hg', 'Pb', 'Bi', 'Th', 'U', 'Np', 'Pu',
}

HEADER = ["name", "overlapping", "under_carbon", "over_carbon", "isolated", "misplaced_hydro"]


def chen_manz(path):
    name = Path(path).stem
    try:
        atoms = read(path)
        sym = np.array(atoms.get_chemical_symbols())
        r = np.array([ATR.get(s, DEFAULT_R) for s in sym])

        d = atoms.get_all_distances(mic=True)
        np.fill_diagonal(d, np.inf)
        rsum = r[:, None] + r[None, :]
        bonded = d <= rsum

        overlapping = bool((d < 0.5 * rsum).any())
        isolated = bool((~bonded.any(axis=1)).any())

        under_c = over_c = False
        for a in np.flatnonzero(sym == 'C'):
            nb = np.flatnonzero(bonded[a])
            if any(sym[b] not in COEF for b in nb):
                continue
            bo = 0.0
            for b in nb:
                A, C = COEF[sym[b]]
                v = 10 ** (A * d[a, b] + C)
                bo += min(v, 1.25) if sym[b] == 'H' else v
            under_c |= bo < 3.3
            over_c |= bo >= 5.5

        misplaced_h = False
        for a in np.flatnonzero(sym == 'H'):
            nb = set(sym[d[a] <= rsum[a] + 0.3])
            if nb & {'N', 'O'} and nb & METALS:
                misplaced_h = True
                break

        return name, overlapping, under_c, over_c, isolated, misplaced_h
    except Exception:
        return (name,) + ("unknown",) * 5


def load_done(csv_path):
    if not os.path.exists(csv_path) or os.path.getsize(csv_path) == 0:
        return None
    with open(csv_path, newline="", encoding="utf-8") as f:
        return {row["name"] for row in csv.DictReader(f)}


def main():
    p = argparse.ArgumentParser(description="Chen-Manz structure checker")
    p.add_argument("cif_dir", help="Path to CIFs")
    p.add_argument("output_csv", help="Output CSV Path")
    p.add_argument("-j", "--workers", type=int, default=32)
    args = p.parse_args()

    done = load_done(args.output_csv)
    paths = sorted(Path(args.cif_dir).glob("*.cif"))
    if done:
        paths = [x for x in paths if x.stem not in done]

    with ProcessPoolExecutor(args.workers) as ex, \
         open(args.output_csv, "a", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        if done is None:
            writer.writerow(HEADER)
        futures = [ex.submit(chen_manz, str(x)) for x in paths]
        for fut in tqdm(as_completed(futures), total=len(futures), desc="Processing CIFs"):
            writer.writerow(fut.result())
            f.flush()


if __name__ == "__main__":
    main()