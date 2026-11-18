"""Run one prepared segment in the foreground; no automatic retry or RPM changes."""
from pathlib import Path
import argparse
import datetime
import hashlib
import json
import os
import subprocess
import time


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run', type=Path, required=True)
    p.add_argument('--exe', type=Path, required=True)
    p.add_argument('--threads', type=int, default=8)
    p.add_argument('--runtime-dir', type=Path)
    a = p.parse_args()
    run, exe = a.run.resolve(strict=True), a.exe.resolve(strict=True)
    if not 1 <= a.threads <= 8: p.error('Use 1..8 threads; larger resources need a separate reviewed setup')
    for name in ['setup_receipt.json', 'lica.in', 'grid.dat', 'geometry.in', 'rotor.in',
                 'rotation_scale.in', 'tip_seal.in']:
        if not (run / name).is_file(): p.error('Missing prepared input: ' + name)
    if (run / 'run_receipt.json').exists() or (run / 'solver.log').exists() or list(run.glob('fld[0-9]*')):
        p.error('Run already contains output; prepare a new restart directory')
    if (run / 'STOP_AFTER_STEP').exists(): p.error('Existing STOP_AFTER_STEP request; do not launch')
    lock = run / 'run.lock'
    try:
        with lock.open('x') as f: f.write(str(os.getpid()))
    except FileExistsError: p.error('Run is locked; inspect the owner, do not start a duplicate')
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS=str(a.threads), OMP_DYNAMIC='FALSE', OMP_STACKSIZE='64M',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1')
    if a.runtime_dir: env['PATH'] = str(a.runtime_dir.resolve(strict=True)) + os.pathsep + env.get('PATH', '')
    receipt = {'status': 'STARTING', 'started_utc': datetime.datetime.now(datetime.timezone.utc).isoformat(),
               'OMP_threads': a.threads, 'binary_sha256': hashlib.sha256(exe.read_bytes()).hexdigest(),
               'automatic_retry': False}
    def save():
        (run / 'run_receipt.json').write_text(json.dumps(receipt, indent=2) + '\n')
    start = time.monotonic()
    child = None
    try:
        save()
        with (run / 'solver.log').open('w') as log:
            flags = subprocess.CREATE_NEW_PROCESS_GROUP | subprocess.CREATE_NO_WINDOW if os.name == 'nt' else 0
            child = subprocess.Popen([str(exe)], cwd=run, env=env, stdout=log, stderr=subprocess.STDOUT,
                                     creationflags=flags, start_new_session=os.name != 'nt')
            receipt.update(status='RUNNING', solver_pid=child.pid); save()
            print(f'Running one segment with {a.threads} OpenMP threads. For a saved stop, create {run / "STOP_AFTER_STEP"}', flush=True)
            while True:
                try:
                    code = child.wait(timeout=2)
                    break
                except subprocess.TimeoutExpired: pass
                except KeyboardInterrupt:
                    (run / 'STOP_AFTER_STEP').touch()
                    print('Requested a saved stop; waiting for the solver to finish its step.', flush=True)
        receipt.update(status='EXITED' if code == 0 else 'FAILED_NO_RESTART', exit_code=code,
                       wall_seconds=time.monotonic() - start)
        save()
    except Exception as exc:
        receipt.update(status='LAUNCH_OR_RUN_ERROR', error=str(exc)); save()
        raise
    finally:
        # Keep the lock if an unexpectedly detached child is still running.
        if child is None or child.poll() is not None: lock.unlink(missing_ok=True)
    raise SystemExit(code)


if __name__ == '__main__': main()
