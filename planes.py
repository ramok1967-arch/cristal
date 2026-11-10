"""
planes.py — geometria płaszczyzn sieciowych (hkl), bez GUI.

Płaszczyzna (hkl) o numerze m to zbiór punktów o współrzędnych ułamkowych f:
    h·f1 + k·f2 + l·f3 = m      (m całkowite)
m = 1 -> płaszczyzna odcinająca a/h, b/k, c/l; m = 0 -> przechodzi przez początek.
Blok n1×n2×n3 w układzie ułamkowym jest prostopadłościanem [0,n1]×[0,n2]×[0,n3],
więc przycinanie robimy we współrzędnych ułamkowych, a dopiero wynik
przeliczamy na kartezjańskie (r = f @ M).
"""
import itertools
import numpy as np


def reciprocal(M: np.ndarray) -> np.ndarray:
    """Wiersze: a*, b*, c* (bez czynnika 2π)."""
    return np.linalg.inv(M).T


def metric(M: np.ndarray) -> np.ndarray:
    return M @ M.T


def d_spacing(M, hkl) -> float:
    h = np.asarray(hkl, float)
    Gs = np.linalg.inv(metric(M))
    return 1.0 / np.sqrt(h @ Gs @ h)


def normal_cart(M, hkl) -> np.ndarray:
    """Jednostkowa normalna płaszczyzny (hkl) w układzie kartezjańskim."""
    g = np.asarray(hkl, float) @ reciprocal(M)
    return g / np.linalg.norm(g)


def m_range(hkl, n):
    """Zakres numerów m płaszczyzn przecinających blok (włącznie)."""
    h = np.asarray(hkl, float)
    corners = np.array(list(itertools.product(*[(0, ni) for ni in n])), float)
    v = corners @ h
    return int(np.ceil(v.min() - 1e-9)), int(np.floor(v.max() + 1e-9))


def _box_edges(n):
    n1, n2, n3 = n
    c = np.array(list(itertools.product((0, n1), (0, n2), (0, n3))), float)
    edges = []
    for i, j in itertools.combinations(range(8), 2):
        if np.sum(c[i] != c[j]) == 1:
            edges.append((c[i], c[j]))
    return edges


def plane_polygon_frac(hkl, m, n, tol=1e-9):
    """
    Wielokąt przecięcia płaszczyzny h·f = m z blokiem (współrzędne ułamkowe).
    Zwraca (K,3) — wierzchołki nieuporządkowane; K < 3 -> brak przecięcia.
    """
    h = np.asarray(hkl, float)
    pts = []
    for p0, p1 in _box_edges(n):
        d = p1 - p0
        den = h @ d
        num = m - h @ p0
        if abs(den) < tol:
            if abs(num) < tol:          # krawędź leży w płaszczyźnie
                pts += [p0, p1]
            continue
        t = num / den
        if -tol <= t <= 1 + tol:
            pts.append(p0 + np.clip(t, 0, 1) * d)
    if not pts:
        return np.zeros((0, 3))
    pts = np.array(pts)
    _, idx = np.unique(np.round(pts, 8), axis=0, return_index=True)
    return pts[np.sort(idx)]


def order_polygon(P: np.ndarray, normal: np.ndarray) -> np.ndarray:
    """Porządkuje wierzchołki wielokąta płaskiego (kartezjańskie) wokół środka."""
    c = P.mean(axis=0)
    u = P[0] - c
    if np.linalg.norm(u) < 1e-12:
        u = np.cross(normal, [1, 0, 0])
        if np.linalg.norm(u) < 1e-6:
            u = np.cross(normal, [0, 1, 0])
    u /= np.linalg.norm(u)
    v = np.cross(normal, u)
    ang = np.arctan2((P - c) @ v, (P - c) @ u)
    return P[np.argsort(ang)]


def plane_polygon(M, hkl, m, n):
    """Uporządkowany wielokąt (kartezjański) lub None."""
    Pf = plane_polygon_frac(hkl, m, n)
    if len(Pf) < 3:
        return None
    P = Pf @ M
    return order_polygon(P, normal_cart(M, hkl))


def on_plane_mask(frac, hkl, ms, tol=1e-6):
    """Maska węzłów leżących na którejkolwiek z płaszczyzn h·f = m, m ∈ ms."""
    v = frac @ np.asarray(hkl, float)
    ms = np.asarray(list(ms), float)
    if ms.size == 0:
        return np.zeros(len(frac), bool)
    return np.any(np.abs(v[:, None] - ms[None, :]) < tol, axis=1)


def intercepts_text(hkl):
    """Odcinki na osiach dla m = 1, np. '1/2 a, ∞ b, 1 c'."""
    out = []
    for ind, ax in zip(hkl, "abc"):
        if ind == 0:
            out.append(f"∞ {ax}")
        elif abs(ind) == 1:
            out.append(f"{'-' if ind < 0 else ''}{ax}")
        else:
            out.append(f"{'-' if ind < 0 else ''}{ax}/{abs(ind)}")
    return ", ".join(out)


def fmt_hkl(hkl):
    """(1 -1 0) -> '(1 1̄ 0)' z nadkreśleniem."""
    s = []
    for v in hkl:
        s.append(f"{abs(v)}̄" if v < 0 else str(v))
    return "(" + " ".join(s) + ")"


# ----------------------------------------------------------------------
# wybór numeru m

def valid_m_range(hkl, n):
    """
    Zakres m, dla których płaszczyzna tnie blok wzdłuż wielokąta (≥3 wierzchołki).
    Pomija płaszczyzny dotykające bloku tylko w narożniku lub wzdłuż krawędzi.
    Zwraca (lo, hi) albo None.
    """
    lo, hi = m_range(hkl, n)
    ok = [m for m in range(lo, hi + 1) if len(plane_polygon_frac(hkl, m, n)) >= 3]
    return (min(ok), max(ok)) if ok else None


def default_m(hkl, n):
    """m = 1 (klasyczna płaszczyzna a/h, b/k, c/l), jeśli tnie blok;
    w przeciwnym razie płaszczyzna najbliższa środka bloku."""
    r = valid_m_range(hkl, n)
    if r is None:
        return None
    lo, hi = r
    if lo <= 1 <= hi:
        return 1
    mc = int(np.floor(np.asarray(hkl, float) @ (np.asarray(n, float) / 2) + 0.5))
    return int(np.clip(mc, lo, hi))


# ----------------------------------------------------------------------
# dwie płaszczyzny: oś pasa, kąt, linia przecięcia

def reduce_int(v):
    v = np.asarray(v, int)
    g = np.gcd.reduce(np.abs(v[v != 0])) if np.any(v) else 1
    return tuple(int(x) for x in v // g)


def zone_axis(h1, h2):
    """[uvw] = (h1 k1 l1) × (h2 k2 l2), skrócone przez NWD; (0,0,0) gdy równoległe."""
    return reduce_int(np.cross(np.asarray(h1, int), np.asarray(h2, int)))


def angle_planes(M, h1, h2) -> float:
    """Kąt między płaszczyznami (= kąt między normalnymi), w stopniach, 0..180."""
    Gs = np.linalg.inv(metric(M))
    a, b = np.asarray(h1, float), np.asarray(h2, float)
    c = (a @ Gs @ b) / np.sqrt((a @ Gs @ a) * (b @ Gs @ b))
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def angle_dirs(M, u1, u2) -> float:
    """Kąt między kierunkami [uvw] (tensor metryczny sieci prostej)."""
    G = metric(M)
    a, b = np.asarray(u1, float), np.asarray(u2, float)
    c = (a @ G @ b) / np.sqrt((a @ G @ a) * (b @ G @ b))
    return float(np.degrees(np.arccos(np.clip(c, -1, 1))))


def dir_length(M, uvw) -> float:
    u = np.asarray(uvw, float)
    return float(np.sqrt(u @ metric(M) @ u))


def intersection_segment_frac(h1, m1, h2, m2, n, tol=1e-9):
    """
    Odcinek przecięcia płaszczyzn h1·f = m1 i h2·f = m2 przycięty do bloku
    (współrzędne ułamkowe). Zwraca (2,3) albo None.
    """
    a, b = np.asarray(h1, float), np.asarray(h2, float)
    t = np.cross(a, b)                       # kierunek = [uvw] we wsp. ułamkowych
    if np.linalg.norm(t) < tol:
        return None
    A = np.vstack([a, b])
    f0 = np.linalg.lstsq(A, np.array([m1, m2], float), rcond=None)[0]
    smin, smax = -np.inf, np.inf
    for i in range(3):                       # metoda slabów
        lo, hi = 0.0, float(n[i])
        if abs(t[i]) < tol:
            if f0[i] < lo - 1e-7 or f0[i] > hi + 1e-7:
                return None
            continue
        s1, s2 = (lo - f0[i]) / t[i], (hi - f0[i]) / t[i]
        smin, smax = max(smin, min(s1, s2)), min(smax, max(s1, s2))
    if smax - smin < 1e-9:
        return None
    return np.array([f0 + smin * t, f0 + smax * t])


def fmt_uvw(uvw):
    s = [f"{abs(v)}̄" if v < 0 else str(v) for v in uvw]
    return "[" + " ".join(s) + "]"


# ----------------------------------------------------------------------
# wygaszenia systematyczne wynikające z centrowania sieci

EXT_RULE = {
    "P": "brak warunków",
    "R": "brak warunków (komórka prymitywna)",
    "C": "h + k = 2n",
    "I": "h + k + l = 2n",
    "F": "h, k, l jednakowej parzystości",
    "R_hex": "−h + k + l = 3n",
}


def is_allowed(centering, hkl) -> bool:
    """Czy refleks hkl jest dozwolony przez centrowanie sieci."""
    h, k, l = (int(x) for x in hkl)
    if centering == "C":
        return (h + k) % 2 == 0
    if centering == "I":
        return (h + k + l) % 2 == 0
    if centering == "F":
        return h % 2 == k % 2 == l % 2
    if centering == "R_hex":
        return (-h + k + l) % 3 == 0
    return True


def lattice_plane_multiple(centering, hkl, nmax=6) -> int:
    """
    Najmniejsze n, dla którego n·(hkl) jest dozwolone.
    Rodzina n·(hkl) to faktyczne płaszczyzny sieciowe (zawierające wszystkie węzły),
    o odstępie d_hkl / n. Dla refleksu dozwolonego n = 1.
    """
    for n in range(1, nmax + 1):
        if is_allowed(centering, [n * x for x in hkl]):
            return n
    return 1


def fmt_m(m) -> str:
    """Numer płaszczyzny: całkowity albo ułamek (np. 1/2)."""
    from fractions import Fraction
    f = Fraction(m).limit_denominator(12)
    return str(f.numerator) if f.denominator == 1 else f"{f.numerator}/{f.denominator}"


# ----------------------------------------------------------------------
# notacja Millera-Bravais (układ heksagonalny)

def _bar(v):
    return f"{abs(v)}̄" if v < 0 else str(v)


def hkil(hkl):
    h, k, l = hkl
    return (h, k, -(h + k), l)


def fmt_hkil(hkl):
    return "(" + " ".join(_bar(v) for v in hkil(hkl)) + ")"


def uvw_to_UVTW(uvw):
    """[uvw] (3 wskaźniki) -> [UVTW]: U=(2u−v)/3, V=(2v−u)/3, T=−(u+v)/3, W=w;
    wynik sprowadzony do najmniejszych liczb całkowitych."""
    u, v, w = (int(x) for x in uvw)
    return reduce_int([2 * u - v, 2 * v - u, -(u + v), 3 * w])


def fmt_UVTW(uvw):
    return "[" + " ".join(_bar(v) for v in uvw_to_UVTW(uvw)) + "]"


# ----------------------------------------------------------------------
# sieć R: przeliczanie wskaźników między ustawieniami (obverse)

P_R2H = np.array([[1, -1, 0], [0, 1, -1], [1, 1, 1]])      # h_H = P · h_R


def hkl_R_to_H(hkl):
    return tuple(int(x) for x in P_R2H @ np.asarray(hkl, int))


def hkl_H_to_R(hkl):
    """h_R = P⁻¹ · h_H; gdy wynik ułamkowy (refleks wygaszony w R_hex),
    zwraca najmniejszą całkowitą wielokrotność (×3) — ta sama orientacja płaszczyzny."""
    v = np.linalg.solve(P_R2H.astype(float), np.asarray(hkl, float))
    if not np.allclose(v, np.round(v), atol=1e-9):
        v = v * 3
    return tuple(int(x) for x in np.round(v))
