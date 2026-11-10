# Cristal — sieci Bravais (etapy 1–4, 3 płaszczyzny)

## Instalacja (jednorazowo)
    conda create -n cristal python=3.11 -y
    conda activate cristal
    pip install -r requirements.txt

## Uruchomienie
    conda activate cristal
    cd C:\AI\cristal
    python cristal.py

## Pliki
- lattice.py  — geometria (macierz komórki, centrowania, węzły bloku), bez GUI
- planes.py   — płaszczyzny (hkl): d_hkl, normalna, przycięcie do bloku
- cristal.py  — okno PyQt5 + PyVista

## Obsługa
- lewy przycisk myszy: obrót, kółko: zoom, Shift+LPM / środkowy: przesunięcie
- parametry niezależne dla danego układu są aktywne, pozostałe wymusza symetria
- jednoskośny: oś unikalna b; trygonalny R: komórka romboedryczna
- "wzdłuż a/b/c" włącza rzut równoległy

## Plan
1. [x] sieci + obrót + blok n1×n2×n3
2. [x] płaszczyzna (hkl) przycięta do bloku
3. [x] rodzina płaszczyzn, d_hkl
4. [x] dwie płaszczyzny: oś pasa [uvw], kąt, widok wzdłuż osi pasa
5. [ ] wygaszenia, notacja hkil (Miller-Bravais), sieć R w ustawieniu heksagonalnym

## Płaszczyzny
- płaszczyzna (hkl) nr m: h·x + k·y + l·z = m (współrzędne ułamkowe);
  m = 1 odcina a/h, b/k, c/l; m = 0 przechodzi przez początek
- spinbox m obejmuje tylko płaszczyzny tnące blok (bez dotykających narożnika/krawędzi)
- po zmianie hkl m dobiera się samo: m = 1, jeśli tnie blok, inaczej płaszczyzna
  najbliższa środka bloku
- "cała rodzina" rysuje wszystkie m przecinające blok
- indeksy dotyczą komórki konwencjonalnej; w sieciach centrowanych część
  węzłów leży między płaszczyznami o całkowitym m (np. (100) w cF) — patrz etap 5

## Dwie / trzy płaszczyzny
- kąt: cos φ = h1ᵀG*h2 / √(h1ᵀG*h1 · h2ᵀG*h2), G* — tensor metryczny sieci odwrotnej
- oś pasa [uvw] = (h1k1l1) × (h2k2l2) (prawo pasowe Weissa), skrócona przez NWD
- do 3 płaszczyzn; dla każdej pary: kąt, oś pasa, linia przecięcia
  (kolory linii: 1∧2 czerwona, 1∧3 pomarańczowa, 2∧3 niebieska);
  w trybie rodzin — wszystkie przecięcia
- trzy płaszczyzny: det[h1;h2;h3] = 0 -> leżą w jednym pasie (tautozonalne);
  w przeciwnym razie czarna kulka = punkt wspólny (tryb „jedna”)
- lista + "wzdłuż osi pasa" ustawia widok wzdłuż [uvw] w rzucie równoległym (jak oś strefy w TEM)
