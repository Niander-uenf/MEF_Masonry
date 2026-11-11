"""mef_masonry — thermal and thermo-mechanical finite elements for block
masonry (BCEPQ, structural and infill clay units, concrete), with round holes,
rectangular cells or horizontal holes, grout, dry or mortar joints,
in pure Python
(NumPy, SciPy, Matplotlib).

Modules:
  parameters  properties, geometry and test protocol (with their sources)
  mesh        2D section with holes (O-grid mesh), rectangular cells or vertical section
              of horizontal holes, and 3D prism by extrusion (W6, H8, W15), with
              grout, mortar and horizontal holes course by course
  elements    fast element matrices (T3 triangle and W6 wedge)
  elements3d  generic isoparametric elements W6, H8 and W15
  nonlinear   Mazars damage in the unit, mortar and grout; dry joints with gaps and friction;
              transient creep (Anderberg and Thelandersson) during heating under load
  wall        wall under load and fire: corotational beams with fibre sections and transient creep
  wall_data   curves of the BCEPQ prism (30 to 1000 °C) used as the wall law
  sensitivity global sensitivity analysis of the wall under ISO 834 (Morris method)
  library     units, materials, mortars and grouts with their sources (fill in the project)
  project     project of the program (schema, validation, execution)
  visual      export to the WebGL post-processor
  thermal     transient conduction with convection and radiation
  mechanical  3D thermo-elasticity and 2D generalised plane strain
  curves      furnace and ISO 834 curves
  studies     joint calibration, 3D UFV test, section studies
  output      VTK (ParaView), CSV and figures
  legacy      compatibility with version 1.x (Portuguese names of modules, functions and
              project fields)
"""
__version__ = "2.0.0"

from . import legacy as _legacy

_legacy.install()   # import mef_masonry.projeto etc. (version 1.x names) still works
