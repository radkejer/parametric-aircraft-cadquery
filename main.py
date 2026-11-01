#!/usr/bin/env python3
"""Parametric aircraft geometry generator using CadQuery."""

from pathlib import Path
import math
import numpy as np
import yaml
import cadquery as cq

# ---------- Helpers ----------
def deg2rad(d): return math.radians(d)

def naca4(code: str, n=101):
    code = str(code).zfill(4)
    m, p, t = int(code[0])/100, int(code[1])/10, int(code[2:])/100
    beta = np.linspace(0, np.pi, n)
    x = (1 - np.cos(beta)) / 2
    yt = 5*t*(0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)
    yc = np.zeros_like(x)
    dyc = np.zeros_like(x)
    if m > 0 and p > 0:
        mask = x < p
        yc[mask]  = m/p**2 * (2*p*x[mask] - x[mask]**2)
        yc[~mask] = m/(1-p)**2 * ((1-2*p) + 2*p*x[~mask] - x[~mask]**2)
        dyc[mask]  = 2*m/p**2 * (p - x[mask])
        dyc[~mask] = 2*m/(1-p)**2 * (p - x[~mask])
    th = np.arctan(dyc)
    xu, yu = x - yt*np.sin(th), yc + yt*np.cos(th)
    xl, yl = x + yt*np.sin(th), yc - yt*np.cos(th)
    return np.column_stack([np.concatenate([xu[::-1], xl[1:]]),
                            np.concatenate([yu[::-1], yl[1:]])])

def get_airfoil(spec, n=101):
    if Path(spec).suffix.lower() == ".dat":
        pts = np.loadtxt(spec, skiprows=1)
        return pts
    return naca4(spec, n)

def make_wing(p):
    S, b, taper = p["S"], p["b"], p["taper"]
    c_root = 2*S / (b*(1+taper))
    c_tip  = taper * c_root
    y_pos = np.linspace(0, b/2, p.get("n_sections", 7))
    sections = []
    for y in y_pos:
        eta = y / (b/2)
        c = c_root*(1-eta) + c_tip*eta
        af = get_airfoil(p.get("airfoil_root", "2412"))
        if "airfoil_tip" in p:
            af = (1-eta)*af + eta*get_airfoil(p["airfoil_tip"])
        af = af * c
        # simple placement (sweep + dihedral + incidence)
        pts = [cq.Vector(x + y*math.tan(deg2rad(p.get("sweep_le",0))),
                         y,
                         z + y*math.tan(deg2rad(p.get("dihedral",0))))
               for x,z in af]
        # rotate for incidence/twist (simplified)
        wire = cq.Workplane("XZ").polyline([(v.x, v.z) for v in pts]).close()
        sections.append(wire.translate((0, y, 0)).val())
    loft = cq.Workplane().add(sections[0])
    for s in sections[1:]:
        loft = loft.add(s)
    solid = loft.loft(ruled=False)
    return solid.mirror("XZ", union=True)

def make_fuselage(p):
    L, W, H = p["length"], p["max_width"], p["max_height"]
    xs = np.linspace(0, L, 9)
    sections = []
    for x in xs:
        r = 0.3 + 0.7 * math.sin(math.pi * x / L)**0.8   # simple shape
        ell = cq.Workplane("YZ").ellipse(r*W/2, r*H/2).translate((x + p.get("x0",0), 0, 0))
        sections.append(ell.val())
    loft = cq.Workplane().add(sections[0])
    for s in sections[1:]:
        loft = loft.add(s)
    return loft.loft(ruled=False)

def build_aircraft(cfg):
    assy = cq.Assembly()
    assy.add(make_wing(cfg["wing"]), name="wing")
    assy.add(make_fuselage(cfg["fuselage"]), name="fuselage")
    # HT & VT reuse wing logic with simplified params
    ht = make_wing({**cfg["ht"], "airfoil_root": cfg["ht"].get("airfoil","0012"),
                    "n_sections": 5, "twist": 0})
    assy.add(ht, name="ht")
    vt = make_wing({**cfg["vt"], "airfoil_root": cfg["vt"].get("airfoil","0012"),
                    "n_sections": 5, "dihedral": 0, "twist": 0})
    vt = vt.rotate((0,0,0), (1,0,0), 90)
    assy.add(vt, name="vt")
    return assy

def main():
    cfg = yaml.safe_load(open("example_aircraft.yaml"))
    print("Building aircraft...")
    assy = build_aircraft(cfg)
    Path("output").mkdir(exist_ok=True)
    assy.save("output/aircraft.step", exportType="STEP")
    cq.exporters.export(assy.toCompound(), "output/aircraft.stl")
    print("Exported → output/aircraft.step and .stl")

if __name__ == "__main__":
    main()
