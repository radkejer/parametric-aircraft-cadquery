# Parametric Aircraft Geometry with CadQuery

Generate 3D aircraft models (wing, fuselage, horizontal & vertical tail) from a simple YAML parameter file using [CadQuery](https://cadquery.readthedocs.io/).

## Features
- Parametric wing (area, span, taper, sweep, dihedral, twist, incidence)
- NACA 4-digit airfoils + support for custom `.dat` files
- Simple lofted fuselage
- Horizontal and vertical tails
- Exports STEP and STL

## Quick Start
```bash
pip install -r requirements.txt
python main.py
