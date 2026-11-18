"""Validate separately supplied private assets and prepare a NEW run; never solve."""
from pathlib import Path
import argparse
import hashlib
import json
import math
import shutil
import struct

ROOT = Path(__file__).resolve().parents[1]


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024 * 1024), b''): h.update(block)
    return h.hexdigest()


def check_function(path, expected):
    with path.open('rb') as f:
        raw = f.read(68)
    if len(raw) != 68: raise ValueError('Truncated geometry header')
    magic, nx, ny, nz, *metrics = struct.unpack('<8s3i6d', raw)
    if magic != b'LICAFN02' or min(nx, ny, nz) < 2:
        raise ValueError('Invalid LICAFN02 geometry')
    if not all(math.isfinite(v) for v in metrics) or min(metrics[3:]) <= 0:
        raise ValueError('Invalid geometry coordinates/spacing')
    if path.stat().st_size != 68 + 4 * nx * ny * nz:
        raise ValueError('Geometry size mismatch')
    digest = sha(path)
    if digest != expected:
        raise ValueError('Geometry hash mismatch: this is a case-specific solver, not an arbitrary STL runner')
    return {'sha256': digest, 'shape': [nx, ny, nz]}


def restart_header(path):
    records = []
    with path.open('rb') as f:
        for _ in range(2):
            marker = f.read(4)
            if len(marker) != 4: raise ValueError('Missing restart record')
            length, = struct.unpack('<i', marker)
            if length != 32: raise ValueError('Unexpected restart precision/record format')
            data = f.read(length)
            if len(data) != length or f.read(4) != marker:
                raise ValueError('Corrupt restart header')
            records.append(data)
    nx, ny, nz, re = struct.unpack('<qqqd', records[0])
    step, segment_step, time, dt = struct.unpack('<qqdd', records[1])
    grid_dims = tuple(map(int, (ROOT / 'setup/grid.dat').read_text().split()[:3]))
    if (nx, ny, nz) != grid_dims: raise ValueError('Restart grid mismatch')
    if not all(math.isfinite(v) for v in (re, time, dt)) or time < 0 or dt <= 0:
        raise ValueError('Invalid restart time/scales')
    expected_re = (.02149 * 5000 * math.pi / 60) * .02149 / 1.5e-5
    if abs(re - expected_re) > 1e-6 * expected_re: raise ValueError('Restart Reynolds scale mismatch')
    return {'step': step, 'segment_step': segment_step, 'time_s': time * 60 / (5000 * math.pi),
            'dt_s': dt * 60 / (5000 * math.pi), 'sha256': sha(path)}


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--fixed', type=Path, required=True)
    p.add_argument('--rotor', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--profile', choices=['startup', 'continuation'], default='startup')
    p.add_argument('--restart', type=Path)
    p.add_argument('--steps', type=int, default=20)
    p.add_argument('--output-every', type=int, default=10)
    p.add_argument('--end-time-s', type=float)
    a = p.parse_args()
    if not 1 <= a.steps <= 20000 or a.output_every < 1:
        p.error('Use 1..20000 steps per segment and a positive output interval')
    if a.profile == 'continuation' and a.restart is None:
        p.error('continuation requires an accepted --restart field')
    if a.end_time_s is not None and (not math.isfinite(a.end_time_s) or a.end_time_s <= 0):
        p.error('--end-time-s must be finite and positive')
    output = a.output.resolve()
    if output.exists(): p.error('Output must be a new directory')
    config = json.loads((ROOT / 'setup/case.json').read_text())
    paths = [a.fixed.resolve(strict=True), a.rotor.resolve(strict=True)]
    if any(len(str(q)) > 1000 or '\n' in str(q) for q in paths): p.error('Geometry path too long/invalid')
    geometry = {name: check_function(path, config['required_private_geometry_sha256'][name])
                for name, path in zip(['fixed.fn', 'rotor.fn'], paths)}
    restart = restart_header(a.restart) if a.restart else None
    if a.profile == 'continuation' and restart['time_s'] < .1473299:
        p.error('Archived larger-step profile starts at 0.14733 s or later')
    start_s = restart['time_s'] if restart else 0.
    if a.end_time_s is not None and a.end_time_s <= start_s:
        p.error('End time must be after restart time')
    dt, cfl = (1e-5, .1) if a.profile == 'startup' else (.0002, .5)
    tref = config['scales']['time_s']
    lines = (ROOT / 'setup/lica.in').read_text().splitlines()
    q = lines[1].split(); q[2] = str(int(restart is not None)); lines[1] = ' '.join(q)
    lines[3] = f'{a.steps} {a.output_every} 100 1'
    lines[5] = f'2 {dt / tref:.17g} {cfl:.17g}'
    lines[14] = 'initial.fld' if restart else 'none'
    output.mkdir(parents=True)
    for name in ['grid.dat', 'rotor.in', 'rotation_scale.in', 'tip_seal.in']:
        shutil.copy2(ROOT / 'setup' / name, output / name)
    (output / 'lica.in').write_text('\n'.join(lines) + '\n', encoding='ascii')
    (output / 'geometry.in').write_text('\n'.join(str(q) for q in paths) + '\n', encoding='utf-8')
    if restart: shutil.copy2(a.restart, output / 'initial.fld')
    if a.end_time_s is not None:
        (output / 'run_limit.in').write_text(f'{a.end_time_s / tref:.17g}\n', encoding='ascii')
    receipt = {'status': 'PREPARED_NOT_LAUNCHED', 'profile': a.profile, 'steps': a.steps,
               'dt_cap_s': dt, 'CFL_target': cfl, 'target_rpm': 500, 'end_time_s': a.end_time_s,
               'geometry': geometry, 'restart': restart,
               'fresh_start_not_bitwise_historical_restart': restart is None}
    (output / 'setup_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    print('PREPARED_NOT_LAUNCHED: ' + str(output))


if __name__ == '__main__': main()
