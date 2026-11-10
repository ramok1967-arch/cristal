"""
cristal.py — przeglądarka 14 sieci Bravais z płaszczyznami (hkl) (etapy 1–4, 3 płaszczyzny).

Uruchomienie:
    python cristal.py

Wymaga: numpy, PyQt5, pyvista, pyvistaqt
"""
import sys
import numpy as np
import pyvista as pv
import vtk
from pyvistaqt import QtInteractor
from PyQt5 import QtCore, QtWidgets

import lattice as L
import planes as PL

PLANE_COLORS = ["#f2b705", "#8e44ad", "#16a085"]
COL_ONPLANE = "#ffd400"
COL_ZONE = "#e0101a"       # linia przecięcia (oś pasa)
# kolory linii przecięcia par płaszczyzn: (1,2), (1,3), (2,3)
PAIR_COLORS = {(0, 1): "#e0101a", (0, 2): "#ff7f0e", (1, 2): "#1f77b4"}
COL_POINT = "#000000"      # punkt wspólny trzech płaszczyzn

COL_CORNER = "#3b6fd8"     # węzły narożne
COL_CENTER = "#e0553a"     # węzły centrujące
COL_CELL = "#111111"       # komórka 000
COL_GRID = "#9a9a9a"       # pozostałe komórki bloku
COL_AXES = {"a": "#d62728", "b": "#2ca02c", "c": "#1f77b4"}

LABEL_FONT = 20            # rozmiar etykiet a, b, c [px]

PARAM_LABEL = {"a": "a [Å]", "b": "b [Å]", "c": "c [Å]",
               "alpha": "α [°]", "beta": "β [°]", "gamma": "γ [°]"}


class PlaneControl(QtWidgets.QGroupBox):
    """Panel jednej płaszczyzny (hkl)."""
    changed = QtCore.pyqtSignal()

    def __init__(self, title, color, hkl=(1, 1, 1), checked=False):
        super().__init__(title)
        self.color = color
        self.setCheckable(True)
        self.setChecked(checked)
        self.setStyleSheet(f"QGroupBox::title{{color:{color}; font-weight:bold;}}")
        gl = QtWidgets.QGridLayout(self)

        gl.addWidget(QtWidgets.QLabel("h k l"), 0, 0)
        self.sp_hkl = []
        for i, v in enumerate(hkl):
            s = QtWidgets.QSpinBox(); s.setRange(-9, 9); s.setValue(v)
            s.valueChanged.connect(self._hkl_changed)
            self.sp_hkl.append(s); gl.addWidget(s, 0, 1 + i)

        self.rb_one = QtWidgets.QRadioButton("jedna, nr m =")
        self.rb_all = QtWidgets.QRadioButton("cała rodzina")
        self.rb_one.setChecked(True)
        self.sp_m = QtWidgets.QSpinBox(); self.sp_m.setRange(-50, 50); self.sp_m.setValue(1)
        tip = ("Numer płaszczyzny w rodzinie (hkl):  h·x + k·y + l·z = m\n"
               "m = 0 — przez początek (0,0,0)\n"
               "m = 1 — odcinki a/h, b/k, c/l\n"
               "m — odległość od początku m·d_hkl\n"
               "Zakres obejmuje tylko płaszczyzny przecinające blok.")
        self.rb_one.setToolTip(tip); self.sp_m.setToolTip(tip)
        self._auto_m = True       # dobierz m automatycznie przy następnym rysowaniu
        gl.addWidget(self.rb_one, 1, 0, 1, 2); gl.addWidget(self.sp_m, 1, 2)
        gl.addWidget(self.rb_all, 1, 3)

        hl = QtWidgets.QHBoxLayout()
        hl.addWidget(QtWidgets.QLabel("krycie"))
        self.sl_op = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.sl_op.setRange(10, 100); self.sl_op.setValue(45)
        hl.addWidget(self.sl_op)
        gl.addLayout(hl, 2, 0, 1, 4)

        self.toggled.connect(lambda _: self.changed.emit())
        self.rb_one.toggled.connect(lambda _: self.changed.emit())
        self.sp_m.valueChanged.connect(lambda _: self.changed.emit())
        self.sl_op.valueChanged.connect(lambda _: self.changed.emit())

    def _hkl_changed(self, *_):
        self._auto_m = True       # nowe hkl -> dobierz m tak, żeby płaszczyzna była widoczna
        self.changed.emit()

    @property
    def hkl(self):
        return tuple(s.value() for s in self.sp_hkl)

    @property
    def active(self):
        return self.isChecked() and any(self.hkl)

    def ms(self, n):
        """Numery m płaszczyzn do narysowania (tylko przecinające blok)."""
        r = PL.valid_m_range(self.hkl, n)
        if r is None:
            return []
        lo, hi = r
        if self.rb_all.isChecked():
            return list(range(lo, hi + 1))
        return [self.sp_m.value()] if lo <= self.sp_m.value() <= hi else []

    def update_m(self, n):
        """Ustawia zakres spinboxa m; przy nowym hkl (lub m poza zakresem) dobiera m."""
        if not any(self.hkl):
            return
        r = PL.valid_m_range(self.hkl, n)
        if r is None:
            return
        lo, hi = r
        cur = self.sp_m.value()
        self.sp_m.blockSignals(True)
        self.sp_m.setRange(lo, hi)
        if self._auto_m or not (lo <= cur <= hi):
            self.sp_m.setValue(PL.default_m(self.hkl, n))
        self.sp_m.blockSignals(False)
        self._auto_m = False


class MainWindow(QtWidgets.QMainWindow):
    def __init__(self):
        super().__init__()
        self.setWindowTitle("Cristal — sieci Bravais")
        self.resize(1300, 850)

        central = QtWidgets.QWidget()
        self.setCentralWidget(central)
        lay = QtWidgets.QHBoxLayout(central)

        panel = QtWidgets.QWidget()
        scroll = QtWidgets.QScrollArea()
        scroll.setWidget(panel)
        scroll.setWidgetResizable(True)
        scroll.setFixedWidth(330)
        scroll.setHorizontalScrollBarPolicy(QtCore.Qt.ScrollBarAlwaysOff)
        lay.addWidget(scroll)
        pl = QtWidgets.QVBoxLayout(panel)

        self.plotter = QtInteractor(central)
        self.plotter.set_background("white")
        lay.addWidget(self.plotter.interactor, 1)

        # --- wybór sieci
        g = QtWidgets.QGroupBox("Sieć Bravais")
        gl = QtWidgets.QVBoxLayout(g)
        self.cb_lattice = QtWidgets.QComboBox()
        for bt in L.BRAVAIS:
            self.cb_lattice.addItem(bt.label, bt.key)
        self.cb_lattice.setCurrentIndex(12)   # cI na start
        gl.addWidget(self.cb_lattice)
        pl.addWidget(g)

        # --- parametry komórki
        g = QtWidgets.QGroupBox("Parametry komórki")
        fl = QtWidgets.QFormLayout(g)
        self.sp = {}
        for name in L.PARAM_NAMES:
            s = QtWidgets.QDoubleSpinBox()
            s.setDecimals(3)
            s.setKeyboardTracking(False)
            if name in ("a", "b", "c"):
                s.setRange(0.5, 50.0); s.setSingleStep(0.1)
            else:
                s.setRange(20.0, 160.0); s.setSingleStep(1.0)
            self.sp[name] = s
            fl.addRow(PARAM_LABEL[name], s)
        self.bt_default = QtWidgets.QPushButton("Wartości domyślne")
        fl.addRow(self.bt_default)
        pl.addWidget(g)

        # --- blok komórek
        g = QtWidgets.QGroupBox("Blok komórek  n₁ × n₂ × n₃")
        hl = QtWidgets.QHBoxLayout(g)
        self.sp_n = []
        for _ in range(3):
            s = QtWidgets.QSpinBox(); s.setRange(1, 8); s.setValue(1)
            self.sp_n.append(s); hl.addWidget(s)
        pl.addWidget(g)

        # --- wygląd
        g = QtWidgets.QGroupBox("Wygląd")
        vl = QtWidgets.QVBoxLayout(g)
        self.chk_grid = QtWidgets.QCheckBox("Krawędzie wszystkich komórek"); self.chk_grid.setChecked(True)
        self.chk_axes = QtWidgets.QCheckBox("Wektory a, b, c"); self.chk_axes.setChecked(True)
        self.chk_color = QtWidgets.QCheckBox("Wyróżnij węzły centrujące"); self.chk_color.setChecked(True)
        self.chk_ortho = QtWidgets.QCheckBox("Rzut równoległy"); self.chk_ortho.setChecked(False)
        for w in (self.chk_grid, self.chk_axes, self.chk_color, self.chk_ortho):
            vl.addWidget(w)
        hl = QtWidgets.QHBoxLayout()
        hl.addWidget(QtWidgets.QLabel("Rozmiar węzłów"))
        self.sl_size = QtWidgets.QSlider(QtCore.Qt.Horizontal)
        self.sl_size.setRange(2, 30); self.sl_size.setValue(2)
        hl.addWidget(self.sl_size)
        vl.addLayout(hl)
        pl.addWidget(g)

        # --- płaszczyzny (hkl)
        self.planes = [
            PlaneControl("Płaszczyzna 1", PLANE_COLORS[0], (1, 1, 1), checked=True),
            PlaneControl("Płaszczyzna 2", PLANE_COLORS[1], (1, 0, 0), checked=False),
            PlaneControl("Płaszczyzna 3", PLANE_COLORS[2], (0, 0, 1), checked=False),
        ]
        for pc in self.planes:
            pl.addWidget(pc)
            pc.changed.connect(self.on_planes_changed)
        self.chk_onplane = QtWidgets.QCheckBox("Podświetl węzły na płaszczyznach")
        self.chk_onplane.setChecked(True)
        pl.addWidget(self.chk_onplane)
        self.chk_zone = QtWidgets.QCheckBox("Linie przecięcia płaszczyzn (osie pasa)")
        self.chk_zone.setChecked(True)
        pl.addWidget(self.chk_zone)

        # --- widok
        g = QtWidgets.QGroupBox("Widok")
        gl2 = QtWidgets.QGridLayout(g)
        self.bt_va = QtWidgets.QPushButton("wzdłuż a")
        self.bt_vb = QtWidgets.QPushButton("wzdłuż b")
        self.bt_vc = QtWidgets.QPushButton("wzdłuż c")
        self.bt_iso = QtWidgets.QPushButton("izometryczny")
        self.bt_png = QtWidgets.QPushButton("Zapisz PNG…")
        gl2.addWidget(self.bt_va, 0, 0); gl2.addWidget(self.bt_vb, 0, 1); gl2.addWidget(self.bt_vc, 0, 2)
        gl2.addWidget(self.bt_iso, 1, 0, 1, 2); gl2.addWidget(self.bt_png, 1, 2)
        self.cb_zone = QtWidgets.QComboBox()
        self.bt_vzone = QtWidgets.QPushButton("wzdłuż osi pasa")
        self.bt_vzone.setEnabled(False)
        gl2.addWidget(self.cb_zone, 2, 0, 1, 2)
        gl2.addWidget(self.bt_vzone, 2, 2)
        pl.addWidget(g)

        # --- informacje
        self.info = QtWidgets.QLabel()
        self.info.setWordWrap(True)
        self.info.setTextFormat(QtCore.Qt.RichText)
        self.info.setStyleSheet("QLabel{background:#f4f4f4; padding:6px; border:1px solid #ddd;}")
        pl.addWidget(self.info)
        pl.addStretch(1)

        # --- sygnały
        self.cb_lattice.currentIndexChanged.connect(self.on_lattice_changed)
        for s in self.sp.values():
            s.valueChanged.connect(self.on_param_changed)
        for s in self.sp_n:
            s.valueChanged.connect(lambda _: self.redraw(reset_camera=True))
        self.bt_vzone.clicked.connect(self.view_along_zone)
        for w in (self.chk_grid, self.chk_axes, self.chk_color, self.chk_onplane, self.chk_zone):
            w.toggled.connect(lambda _: self.redraw())
        self.chk_ortho.toggled.connect(self.on_ortho)
        self.sl_size.valueChanged.connect(lambda _: self.redraw())
        self.bt_default.clicked.connect(self.load_defaults)
        self.bt_va.clicked.connect(lambda: self.view_along(0))
        self.bt_vb.clicked.connect(lambda: self.view_along(1))
        self.bt_vc.clicked.connect(lambda: self.view_along(2))
        self.bt_iso.clicked.connect(self.view_iso)
        self.bt_png.clicked.connect(self.save_png)

        self._busy = False
        self.M = np.eye(3)
        self._zones = []          # [(etykieta, uvw)]
        self._point3 = None
        self._labels = []         # aktory etykiet a, b, c
        self.on_lattice_changed()

    # ------------------------------------------------------------------
    @property
    def bt(self) -> L.BravaisType:
        return L.BRAVAIS_BY_KEY[self.cb_lattice.currentData()]

    def n_block(self):
        return tuple(s.value() for s in self.sp_n)

    def set_params(self, params):
        self._busy = True
        for name, v in zip(L.PARAM_NAMES, params):
            self.sp[name].setValue(v)
        self._busy = False

    def load_defaults(self):
        self.set_params(self.bt.default)
        self.redraw(reset_camera=True)

    def on_lattice_changed(self, *_):
        bt = self.bt
        for name, s in self.sp.items():
            s.setEnabled(name in bt.free)
        self.set_params(bt.default)
        self.redraw(reset_camera=True)
        self.view_iso()

    def on_param_changed(self, *_):
        if self._busy:
            return
        vals = {k: s.value() for k, s in self.sp.items()}
        params = L.constrain(self.bt, vals)
        self.set_params(params)          # pokaż wartości wymuszone symetrią
        self.redraw()

    def on_planes_changed(self):
        self.redraw()

    def on_ortho(self, on):
        if on:
            self.plotter.enable_parallel_projection()
        else:
            self.plotter.disable_parallel_projection()
        self.plotter.render()

    # ------------------------------------------------------------------
    def redraw(self, reset_camera=False):
        bt = self.bt
        params = L.constrain(bt, {k: s.value() for k, s in self.sp.items()})
        try:
            M = L.cell_matrix(*params)
        except ValueError as e:
            self.info.setText(f"<b style='color:#c00'>{e}</b>")
            return
        self.M = M
        n = self.n_block()
        p = self.plotter
        p.clear_actors()
        self._labels = []

        # węzły
        frac, kind = L.lattice_nodes(bt.centering, n)
        xyz = frac @ M
        r = 0.01 * self.sl_size.value() * min(params[:3])

        # płaszczyzny aktywne: (kontrolka, hkl, lista m)
        active = []
        for pc in self.planes:
            pc.update_m(n)
            if pc.active:
                active.append((pc, pc.hkl, pc.ms(n)))

        onplane = np.zeros(len(frac), bool)
        if self.chk_onplane.isChecked():
            for _, hkl, ms in active:
                onplane |= PL.on_plane_mask(frac, hkl, ms)

        sphere = pv.Sphere(radius=r, theta_resolution=20, phi_resolution=20)
        sphere_big = pv.Sphere(radius=1.25 * r, theta_resolution=20, phi_resolution=20)
        if self.chk_color.isChecked():
            groups = [((kind == 0) & ~onplane, COL_CORNER, sphere),
                      ((kind > 0) & ~onplane, COL_CENTER, sphere)]
        else:
            groups = [(~onplane, COL_CORNER, sphere)]
        groups.append((onplane, COL_ONPLANE, sphere_big))
        for i, (mask, col, geom) in enumerate(groups):
            if mask.any():
                pts = pv.PolyData(xyz[mask])
                p.add_mesh(pts.glyph(geom=geom, scale=False, orient=False),
                           color=col, smooth_shading=True, specular=0.4, name=f"nodes{i}")

        # wielokąty płaszczyzn
        for j, (pc, hkl, ms) in enumerate(active):
            polys = [PL.plane_polygon(M, hkl, m, n) for m in ms]
            polys = [q for q in polys if q is not None]
            if not polys:
                continue
            pts, faces, off = [], [], 0
            for q in polys:
                pts.append(q)
                faces += [len(q)] + list(range(off, off + len(q)))
                off += len(q)
            mesh = pv.PolyData(np.vstack(pts), faces=np.array(faces))
            p.add_mesh(mesh, color=pc.color, opacity=pc.sl_op.value() / 100.0,
                       show_edges=False, lighting=False, name=f"plane{j}")
            p.add_mesh(mesh.extract_feature_edges(boundary_edges=True, feature_edges=False,
                                                  manifold_edges=False, non_manifold_edges=False),
                       color=pc.color, line_width=2.5, name=f"plane{j}_edge")

        # osie pasa: linie przecięcia każdej pary płaszczyzn
        idx_active = [k for k, pc in enumerate(self.planes) if pc.active]
        by_k = {k: (pc, hkl, ms) for k, (pc, hkl, ms) in zip(idx_active, active)}
        zones = []
        for (i, j) in PAIR_COLORS:
            if i not in by_k or j not in by_k:
                continue
            (_, h1, ms1), (_, h2, ms2) = by_k[i], by_k[j]
            uvw = PL.zone_axis(h1, h2)
            if not any(uvw):
                continue
            zones.append((f"{i + 1}∧{j + 1}  {PL.fmt_uvw(uvw)}", uvw))
            if self.chk_zone.isChecked():
                segs = []
                for m1 in ms1:
                    for m2 in ms2:
                        sg = PL.intersection_segment_frac(h1, m1, h2, m2, n)
                        if sg is not None:
                            segs.append(sg)
                    if len(segs) >= 600:
                        break
                if segs:
                    P = np.vstack(segs) @ M
                    p.add_mesh(pv.line_segments_from_points(P), color=PAIR_COLORS[(i, j)],
                               line_width=6, render_lines_as_tubes=True, name=f"zone{i}{j}")

        # punkt wspólny trzech płaszczyzn (tylko tryb „jedna”)
        self._point3 = None
        if len(active) == 3 and all(len(ms) == 1 for _, _, ms in active):
            H = np.array([hkl for _, hkl, _ in active], float)
            if abs(np.linalg.det(H)) > 1e-9:
                f = np.linalg.solve(H, np.array([ms[0] for _, _, ms in active], float))
                self._point3 = f
                if np.all(f >= -1e-7) and np.all(f <= np.array(n) + 1e-7):
                    p.add_mesh(pv.Sphere(radius=max(1.6 * r, 0.06 * min(params[:3])),
                                         center=f @ M),
                               color=COL_POINT, name="point3")

        # lista osi pasa do widoku
        cur = self.cb_zone.currentIndex()
        self._zones = zones
        self.cb_zone.blockSignals(True)
        self.cb_zone.clear()
        for lab, _ in zones:
            self.cb_zone.addItem(lab)
        if 0 <= cur < len(zones):
            self.cb_zone.setCurrentIndex(cur)
        self.cb_zone.blockSignals(False)
        self.bt_vzone.setEnabled(bool(zones))
        self.cb_zone.setEnabled(bool(zones))

        # krawędzie bloku
        if self.chk_grid.isChecked() and n != (1, 1, 1):
            segs = L.grid_segments(n).reshape(-1, 3) @ M
            p.add_mesh(pv.line_segments_from_points(segs), color=COL_GRID, line_width=1.5, name="grid")
        cell = L.unit_cell_segments().reshape(-1, 3) @ M
        p.add_mesh(pv.line_segments_from_points(cell), color=COL_CELL, line_width=3.5, name="cell")

        # wektory a, b, c
        if self.chk_axes.isChecked():
            for vec, lab in zip(M, "abc"):
                ln = np.linalg.norm(vec)
                arr = pv.Arrow(start=(0, 0, 0), direction=vec, scale=ln,
                               tip_length=0.18, tip_radius=0.05, shaft_radius=0.018)
                p.add_mesh(arr, color=COL_AXES[lab], name=f"ax_{lab}")
            # etykiety jako billboardy 3D: w odróżnieniu od napisów 2D (add_point_labels)
            # nie znikają przy zapisie PNG w powiększonej rozdzielczości
            for pos, lab in zip(M * 1.12, "abc"):
                ta = vtk.vtkBillboardTextActor3D()
                ta.SetInput(lab)
                ta.SetPosition(*pos)
                tp = ta.GetTextProperty()
                tp.SetFontSize(LABEL_FONT); tp.SetBold(True); tp.SetColor(0, 0, 0)
                tp.SetJustificationToCentered(); tp.SetVerticalJustificationToCentered()
                p.add_actor(ta, name=f"ax_label_{lab}")
                self._labels.append(ta)

        if reset_camera:
            p.reset_camera()
        p.render()
        self.update_info(bt, params, M, len(frac), active, frac)

    def update_info(self, bt, params, M, nnodes, active=(), frac=None):
        a, b, c, al, be, ga = params
        V = L.cell_volume(M)
        z = L.nodes_per_cell(bt.centering)
        txt = (
            f"<b>{bt.system}</b>, {L.CENTERING_NAME[bt.centering]}<br>"
            f"symbol: <b>{bt.key}</b><br>"
            f"a={a:.3f}  b={b:.3f}  c={c:.3f} Å<br>"
            f"α={al:.2f}°  β={be:.2f}°  γ={ga:.2f}°<br>"
            f"V = {V:.3f} Å³<br>"
            f"węzłów / komórkę: {z}  (V<sub>prym</sub> = {V / z:.3f} Å³)<br>"
            f"węzłów w bloku: {nnodes}"
        )
        for pc, hkl, ms in active:
            d = PL.d_spacing(M, hkl)
            nrm = PL.normal_cart(M, hkl)
            non = int(PL.on_plane_mask(frac, hkl, ms).sum()) if frac is not None else 0
            txt += (
                f"<hr><b style='color:{pc.color}'>{PL.fmt_hkl(hkl)}</b><br>"
                f"d<sub>hkl</sub> = {d:.4f} Å<br>"
                f"odcinki (m=1): {PL.intercepts_text(hkl)}<br>"
                f"normalna (kart.): [{nrm[0]:.3f}, {nrm[1]:.3f}, {nrm[2]:.3f}]<br>"
                f"płaszczyzn w bloku: {len(ms)}"
                f"{'  (m = ' + ', '.join(map(str, ms)) + ')' if 0 < len(ms) <= 8 else ''}<br>"
                + (f"m = {ms[0]}: położenie m·d = {ms[0] * d:+.4f} Å<br>" if len(ms) == 1 else "")
                + f"węzłów na nich: {non}"
            )
        idx_active = [k for k, pc in enumerate(self.planes) if pc.active]
        by_k = {k: hkl for k, (_, hkl, _) in zip(idx_active, active)}
        for (i, j), col in PAIR_COLORS.items():
            if i not in by_k or j not in by_k:
                continue
            h1, h2 = by_k[i], by_k[j]
            phi = PL.angle_planes(M, h1, h2)
            txt += (f"<hr><b>{i + 1}∧{j + 1}: {PL.fmt_hkl(h1)} ∧ {PL.fmt_hkl(h2)}</b><br>"
                    f"kąt φ = <b>{phi:.3f}°</b>  (180°−φ = {180 - phi:.3f}°)<br>")
            uvw = PL.zone_axis(h1, h2)
            if any(uvw):
                t = np.asarray(uvw, float) @ M
                t /= np.linalg.norm(t)
                txt += (f"oś pasa: <b style='color:{col}'>{PL.fmt_uvw(uvw)}</b><br>"
                        f"|{PL.fmt_uvw(uvw)}| = {PL.dir_length(M, uvw):.4f} Å<br>"
                        f"kierunek (kart.): [{t[0]:.3f}, {t[1]:.3f}, {t[2]:.3f}]")
            else:
                txt += "płaszczyzny równoległe — oś pasa nieokreślona"
        if len(active) == 3:
            H = np.array([hkl for _, hkl, _ in active], int)
            txt += "<hr><b>trzy płaszczyzny</b><br>"
            if round(np.linalg.det(H)) == 0:
                # det = 0 -> normalne współpłaszczyznowe -> wspólny pas (tautozonalne)
                axes = [PL.zone_axis(H[a], H[b]) for a, b in ((0, 1), (0, 2), (1, 2))]
                ax = next((u for u in axes if any(u)), None)
                txt += ("leżą w jednym pasie" + (f" — oś {PL.fmt_uvw(ax)}" if ax else "")
                        + "<br>(det = 0, normalne współpłaszczyznowe)")
            else:
                txt += f"nie leżą w jednym pasie (det = {round(np.linalg.det(H))})<br>"
                if self._point3 is not None:
                    f = self._point3 + 0.0
                    f[np.abs(f) < 5e-4] = 0.0
                    inside = np.all(f >= -1e-7) and np.all(f <= np.array(self.n_block()) + 1e-7)
                    txt += (f"punkt wspólny (ułamk.): ({f[0]:.3f}, {f[1]:.3f}, {f[2]:.3f})"
                            + ("" if inside else " — poza blokiem"))
        self.info.setText(txt)

    # ------------------------------------------------------------------
    def _center(self):
        return (np.array(self.n_block()) / 2.0) @ self.M

    def view_along(self, i):
        """Patrz wzdłuż wektora a, b lub c."""
        self._view_dir(self.M[i], self.M[(i + 2) % 3])   # „góra”: c dla a/b, b dla c

    def view_along_zone(self):
        """Patrz wzdłuż wybranej osi pasa [uvw] (jak w TEM wzdłuż osi strefy)."""
        k = self.cb_zone.currentIndex()
        if not (0 <= k < len(self._zones)):
            return
        v = np.asarray(self._zones[k][1], float) @ self.M
        dn = v / np.linalg.norm(v)
        # „góra” = wektor a/b/c najbardziej prostopadły do osi pasa
        k = int(np.argmin([abs(dn @ (r / np.linalg.norm(r))) for r in self.M]))
        self._view_dir(v, self.M[k])

    def _view_dir(self, v, up):
        d = v / np.linalg.norm(v)
        up = up - (up @ d) * d
        if np.linalg.norm(up) < 1e-6:
            up = np.array([0, 0, 1.0])
        ctr = self._center()
        dist = 4.0 * np.linalg.norm(np.array(self.n_block()) @ self.M)
        # rzut równoległy — inaczej warstwy wzdłuż osi rozjeżdżają się w perspektywie
        if not self.chk_ortho.isChecked():
            self.chk_ortho.setChecked(True)
        self.plotter.camera_position = [tuple(ctr + d * dist), tuple(ctr), tuple(up / np.linalg.norm(up))]
        self.plotter.reset_camera()
        self.plotter.render()

    def view_iso(self):
        ctr = self._center()
        d = np.array([1.0, 0.7, 0.55]); d /= np.linalg.norm(d)
        dist = 4.0 * np.linalg.norm(np.array(self.n_block()) @ self.M)
        self.plotter.camera_position = [tuple(ctr + d * dist), tuple(ctr), (0, 0, 1)]
        self.plotter.reset_camera()
        self.plotter.render()

    def save_png(self):
        fn, _ = QtWidgets.QFileDialog.getSaveFileName(self, "Zapisz obraz", f"{self.bt.key}.png", "PNG (*.png)")
        if not fn:
            return
        scale = 2
        # billboardy mają rozmiar w pikselach -> przy powiększonym zapisie zwiększ czcionkę
        for ta in self._labels:
            ta.GetTextProperty().SetFontSize(LABEL_FONT * scale)
        try:
            self.plotter.screenshot(fn, transparent_background=False, scale=scale)
        finally:
            for ta in self._labels:
                ta.GetTextProperty().SetFontSize(LABEL_FONT)
            self.plotter.render()

    def closeEvent(self, e):
        self.plotter.close()
        super().closeEvent(e)


def main():
    app = QtWidgets.QApplication(sys.argv)
    w = MainWindow()
    w.show()
    sys.exit(app.exec_())


if __name__ == "__main__":
    main()
