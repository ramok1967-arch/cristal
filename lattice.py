"""
lattice.py — geometria 14 sieci Bravais (bez GUI).

Konwencje:
  - wiersze macierzy komórki M to wektory a, b, c w układzie kartezjańskim,
    r_kart = f @ M  (f — współrzędne ułamkowe)
  - a ∥ x, b w płaszczyźnie xy (standard krystalograficzny)
  - układ jednoskośny: oś unikalna b (kąt β ≠ 90°)
  - sieć R: komórka romboedryczna (a = b = c, α = β = γ)
"""
from dataclasses import dataclass
import numpy as np

# Przesunięcia centrujące (współrzędne ułamkowe)
CENTERING = {
    "P": [(0, 0, 0)],
    "C": [(0, 0, 0), (0.5, 0.5, 0)],
    "I": [(0, 0, 0), (0.5, 0.5, 0.5)],
    "F": [(0, 0, 0), (0.5, 0.5, 0), (0.5, 0, 0.5), (0, 0.5, 0.5)],
    "R": [(0, 0, 0)],          # komórka romboedryczna jest prymitywna
}

CENTERING_NAME = {
    "P": "prymitywna",
    "C": "centrowana na podstawach (C)",
    "I": "przestrzennie centrowana (I)",
    "F": "ściennie centrowana (F)",
    "R": "romboedryczna (prymitywna)",
}


@dataclass(frozen=True)
class BravaisType:
    key: str            # np. "cI"
    system: str         # nazwa układu
    centering: str      # P, C, I, F, R
    free: tuple         # które parametry są niezależne
    default: tuple      # (a, b, c, alpha, beta, gamma)

    @property
    def label(self):
        return f"{self.system} {self.centering}  ({self.key})"


# 14 sieci Bravais
BRAVAIS = [
    BravaisType("aP", "trójskośny",     "P", ("a", "b", "c", "alpha", "beta", "gamma"), (3.0, 3.6, 4.2, 80.0, 70.0, 100.0)),
    BravaisType("mP", "jednoskośny",    "P", ("a", "b", "c", "beta"),                  (3.0, 3.6, 4.2, 90.0, 110.0, 90.0)),
    BravaisType("mC", "jednoskośny",    "C", ("a", "b", "c", "beta"),                  (3.0, 3.6, 4.2, 90.0, 110.0, 90.0)),
    BravaisType("oP", "rombowy",        "P", ("a", "b", "c"),                          (3.0, 3.6, 4.2, 90.0, 90.0, 90.0)),
    BravaisType("oC", "rombowy",        "C", ("a", "b", "c"),                          (3.0, 3.6, 4.2, 90.0, 90.0, 90.0)),
    BravaisType("oI", "rombowy",        "I", ("a", "b", "c"),                          (3.0, 3.6, 4.2, 90.0, 90.0, 90.0)),
    BravaisType("oF", "rombowy",        "F", ("a", "b", "c"),                          (3.0, 3.6, 4.2, 90.0, 90.0, 90.0)),
    BravaisType("tP", "tetragonalny",   "P", ("a", "c"),                               (3.0, 3.0, 4.5, 90.0, 90.0, 90.0)),
    BravaisType("tI", "tetragonalny",   "I", ("a", "c"),                               (3.0, 3.0, 4.5, 90.0, 90.0, 90.0)),
    BravaisType("hR", "trygonalny",     "R", ("a", "alpha"),                           (3.5, 3.5, 3.5, 70.0, 70.0, 70.0)),
    BravaisType("hP", "heksagonalny",   "P", ("a", "c"),                               (3.0, 3.0, 4.9, 90.0, 90.0, 120.0)),
    BravaisType("cP", "regularny",      "P", ("a",),                                   (3.0, 3.0, 3.0, 90.0, 90.0, 90.0)),
    BravaisType("cI", "regularny",      "I", ("a",),                                   (3.0, 3.0, 3.0, 90.0, 90.0, 90.0)),
    BravaisType("cF", "regularny",      "F", ("a",),                                   (3.0, 3.0, 3.0, 90.0, 90.0, 90.0)),
]
BRAVAIS_BY_KEY = {bt.key: bt for bt in BRAVAIS}

PARAM_NAMES = ("a", "b", "c", "alpha", "beta", "gamma")


def constrain(bt: BravaisType, values: dict) -> tuple:
    """Z podanych wartości buduje pełny zestaw (a,b,c,α,β,γ) zgodny z układem."""
    v = dict(zip(PARAM_NAMES, bt.default))
    for k in bt.free:
        if k in values:
            v[k] = float(values[k])
    s = bt.system
    if s in ("tetragonalny", "heksagonalny"):
        v["b"] = v["a"]
    if s == "regularny":
        v["b"] = v["c"] = v["a"]
    if s == "trygonalny":
        v["b"] = v["c"] = v["a"]
        v["beta"] = v["gamma"] = v["alpha"]
    return tuple(v[k] for k in PARAM_NAMES)


def cell_matrix(a, b, c, alpha, beta, gamma) -> np.ndarray:
    """Macierz 3x3, wiersze = wektory a, b, c (kartezjańskie)."""
    al, be, ga = np.radians([alpha, beta, gamma])
    ca, cb, cg, sg = np.cos(al), np.cos(be), np.cos(ga), np.sin(ga)
    cy = (ca - cb * cg) / sg
    cz2 = 1.0 - cb ** 2 - cy ** 2
    if cz2 <= 1e-9:
        raise ValueError("Niefizyczny zestaw kątów (objętość komórki ≤ 0).")
    M = np.array([
        [a, 0.0, 0.0],
        [b * cg, b * sg, 0.0],
        [c * cb, c * cy, c * np.sqrt(cz2)],
    ])
    M[np.abs(M) < 1e-12] = 0.0
    return M


def cell_volume(M: np.ndarray) -> float:
    return float(abs(np.linalg.det(M)))


def lattice_nodes(centering: str, n=(1, 1, 1), tol=1e-6):
    """
    Węzły sieci w bloku n1×n2×n3 komórek konwencjonalnych.
    Zwraca (frac, kind): współrzędne ułamkowe (N,3) oraz indeks przesunięcia
    centrującego (0 = narożnik, 1.. = węzły centrujące).
    """
    n = np.asarray(n, dtype=int)
    offs = np.array(CENTERING[centering], dtype=float)
    ii, jj, kk = np.meshgrid(*(np.arange(0, ni + 1) for ni in n), indexing="ij")
    grid = np.stack([ii.ravel(), jj.ravel(), kk.ravel()], axis=1).astype(float)

    frac_all, kind_all = [], []
    for idx, o in enumerate(offs):
        f = grid + o
        keep = np.all((f >= -tol) & (f <= n + tol), axis=1)
        frac_all.append(f[keep])
        kind_all.append(np.full(keep.sum(), idx))
    frac = np.vstack(frac_all)
    kind = np.concatenate(kind_all)

    # usunięcie duplikatów (dla bezpieczeństwa)
    _, uniq = np.unique(np.round(frac, 6), axis=0, return_index=True)
    uniq.sort()
    return frac[uniq], kind[uniq]


def grid_segments(n=(1, 1, 1)):
    """
    Odcinki krawędzi wszystkich komórek bloku (współrzędne ułamkowe).
    Zwraca tablicę (S, 2, 3): pary końców odcinków.
    Linie siatki biegną przez cały blok — każda linia to jeden odcinek.
    """
    n1, n2, n3 = n
    segs = []
    for j in range(n2 + 1):
        for k in range(n3 + 1):
            segs.append(((0, j, k), (n1, j, k)))
    for i in range(n1 + 1):
        for k in range(n3 + 1):
            segs.append(((i, 0, k), (i, n2, k)))
    for i in range(n1 + 1):
        for j in range(n2 + 1):
            segs.append(((i, j, 0), (i, j, n3)))
    return np.array(segs, dtype=float)


def unit_cell_segments():
    """12 krawędzi pojedynczej komórki (0..1)."""
    return grid_segments((1, 1, 1))


def nodes_per_cell(centering: str) -> int:
    return len(CENTERING[centering])
