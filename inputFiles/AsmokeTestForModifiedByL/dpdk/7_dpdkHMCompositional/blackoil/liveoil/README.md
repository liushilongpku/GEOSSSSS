# Three-phase black-oil dual-continuum HM smoke case

This case is a shortened black-oil version of
`co2_brine_wag_with_mechanical`. It keeps the original four-layer VTK mesh,
partial matrix/fracture dual-continuum configuration, gravity, injection and
production locations, and sequential poromechanics coupling.

The fluid model is now `BlackOilFluid` with three phases and three components:

- phases: oil, gas, water;
- PVT tables: `tables/pvto_bo.txt`, `tables/pvtg_norv_bo.txt`, and
  `tables/pvtw_bo.txt`;
- independent matrix and fracture relative-permeability models;
- three-phase capillary-pressure models on both continua;
- gas and water injection into the fracture reservoir region;
- oil-rich initial composition in both continua.

The purpose is a fast execution check, not a history-matched reservoir study.
The simulation ends after `86400 s` and uses a one-day forced time step. The
case is expected to complete in seconds to minutes on one rank, depending on
the build and hardware.

Run from this directory:

```bash
../../../../../../build/bin/geosx -v -i input.xml -o /tmp/geos_blackoil_deadoil_validate
../../../../../../build/bin/geosx -i input.xml -o /tmp/geos_blackoil_deadoil_run
```

The second command is the actual smoke calculation. A successful run should
report converged `dualFlow`, `solidMech`, and `poroSolver` iterations, zero
time-step cuts, and a final `Finished at` line. VTK output is written to the
specified output directory.
