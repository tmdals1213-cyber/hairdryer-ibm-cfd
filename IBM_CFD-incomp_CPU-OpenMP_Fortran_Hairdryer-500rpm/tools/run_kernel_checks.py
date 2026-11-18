"""Test the unmodified production velocity module against small synthetic walls."""
from pathlib import Path
import argparse
import hashlib
import json
import os
import shutil
import subprocess
from build_reference import ROOT, FLAGS, compiler_env


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--output', type=Path, default=ROOT / 'build_kernel')
    args = parser.parse_args()
    compiler, env = compiler_env(args.compiler)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    source = ROOT / 'src/mod_ibm_surface_velocity.f90'
    for file in [ROOT / 'tests/synthetic_dependencies.f90', source,
                 ROOT / 'tests/test_actual_velocity_module.f90']:
        shutil.copy2(file, output)
    executable = output / ('test_kernel.exe' if os.name == 'nt' else 'test_kernel')
    cmd = [compiler, *FLAGS, '-O0', '-g', '-fcheck=all', '-fbacktrace', '-frecursive',
           'synthetic_dependencies.f90', 'mod_ibm_surface_velocity.f90',
           'test_actual_velocity_module.f90', '-o', str(executable)]
    with (output / 'build.log').open('w', encoding='utf-8') as log:
        subprocess.run(cmd, cwd=output, env=env, stdout=log,
                       stderr=subprocess.STDOUT, check=True, timeout=120)
    checks = []
    for arg, expected in [('', 0), ('failstop', 1)]:
        result = subprocess.run([str(executable)] + ([arg] if arg else []),
                                cwd=output, env=env, capture_output=True, text=True, timeout=30)
        text = result.stdout + result.stderr
        (output / ('checks.log' if not arg else 'expected_failstop.log')).write_text(text, encoding='utf-8')
        if expected == 0:
            passed = result.returncode == 0 and 'ALL_KERNEL_CHECKS_PASS' in text
        else:
            passed = result.returncode != 0 and 'lacks bounded visible fluid support' in text
        checks.append({'mode': arg or 'positive_and_rejection_checks',
                       'exit_code': result.returncode, 'passed': passed})
        if not passed:
            raise SystemExit('CHECK FAILED: ' + text)
    report = {'status': 'PASS', 'production_module_sha256': hashlib.sha256(source.read_bytes()).hexdigest(),
              'real_production_module_compiled_unchanged': True,
              'synthetic_geometry_only': True, 'CFD_solver_launched': False,
              'checks': checks, 'compile_command': cmd}
    (output / 'result.json').write_text(json.dumps(report, indent=2), encoding='utf-8')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
