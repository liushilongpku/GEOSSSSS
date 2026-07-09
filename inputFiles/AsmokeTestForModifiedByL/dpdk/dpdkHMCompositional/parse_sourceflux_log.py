#!/usr/bin/env python3
"""
Parse a GEOS run log into a per-timestep CSV of SourceFlux injection rates.

GEOS's built-in <SourceFluxStatistics writeCSV="1"> overwrites its CSV every
event (SourceFluxStatistics.cpp opens the file in truncate mode with a per-call
local buffer), so the on-disk CSV only ever keeps the *last* stats snapshot.
The full history, however, is printed to the console at logLevel="2".

This script reconstructs the time series from that console log.

Usage:
    geosx -i dpdkHMCompositional_2d_dual_poromechanics_cyclicSourceFlux.xml > run.log 2>&1
    python3 parse_sourceflux_log.py run.log [out.csv]

Edit FLUXES below if you rename the SourceFlux entries.
"""
import re
import sys

# (flux name in the deck) -> (output column, which component of the [co2, water] pair)
FLUXES = {
    'injCo2M':   ('co2_rate_matrix_mol_s',   0),
    'injWaterM': ('water_rate_matrix_mol_s', 1),
    'injCo2F':   ('co2_rate_frac_mol_s',     0),
    'injWaterF': ('water_rate_frac_mol_s',   1),
}

RATE_RE = re.compile(r'\[\s*([-\d.eE]+)\s*,\s*([-\d.eE]+)\s*\]')
TIME_RE = re.compile(r'Time:\s*([\d.eE+]+)\s*s')


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    log = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else 'sourceflux_injection_timeseries.csv'

    cols = [c for c, _ in FLUXES.values()]
    rows = {}
    cur_t = None
    with open(log) as f:
        for line in f:
            m = TIME_RE.search(line)
            if m:
                cur_t = float(m.group(1))
                continue
            if cur_t is None or 'all_regions' not in line:
                continue
            for flux, (col, comp) in FLUXES.items():
                if flux in line:
                    br = RATE_RE.findall(line)
                    if br:  # last bracket on the line is the rate column
                        rows.setdefault(cur_t, {c: 0.0 for c in cols})[col] = float(br[-1][comp])

    with open(out, 'w') as o:
        o.write('time_s,' + ','.join(cols) + ',phase\n')
        for t in sorted(rows):
            r = rows[t]
            gas = any(r[c] for c in cols if 'co2' in c)
            wat = any(r[c] for c in cols if 'water' in c)
            phase = 'gas' if gas else ('water' if wat else 'none')
            o.write(f'{t:g},' + ','.join(f'{r[c]:g}' for c in cols) + f',{phase}\n')
    print(f'wrote {out} ({len(rows)} rows)')


if __name__ == '__main__':
    main()
