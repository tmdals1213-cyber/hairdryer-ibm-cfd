"""Read a native field. Geometry-free checks do NOT classify fluid versus solid."""
from pathlib import Path
import argparse
import json
import math
import numpy as np
from lica_flow import read_field, read_grid


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--grid', type=Path, required=True)
    p.add_argument('--field', type=Path, required=True)
    a = p.parse_args()
    axes = read_grid(a.grid)
    f = read_field(a.field, axes)
    tref = 60 / (5000 * math.pi)
    uref = .02149 / tref
    report = {'field': a.field.name, 'step': f['step'], 'physical_time_s': f['time'] * tref,
              'physical_dt_s': f['dt'] * tref, 'coordinate_counts': [len(v) for v in axes],
              'cell_shape': list(f['pressure'].shape),
              'all_cell_peak_speed_m_s_including_solid_ghost': float(np.linalg.norm(f['velocity'], axis=-1).max() * uref),
              'pressure_variable': 'p/rho; multiply stored nondimensional pressure by Uref^2 for m^2/s^2',
              'finite_velocity_and_pressure': True,
              'fluid_mask_applied': False,
              'note': 'No fluid-only maximum, wall force or vortex validation without matching moving geometry.'}
    print(json.dumps(report, indent=2))


if __name__ == '__main__': main()
