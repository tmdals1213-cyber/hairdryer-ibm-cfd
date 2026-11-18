"""Compile the supplied reference source; NEVER launch a CFD calculation."""
from pathlib import Path
import argparse
import json
import os
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[1]
NAMES = ['mod_common', 'mod_flowarray', 'mod_stl_sdf', 'mod_rotor_motion',
         'mod_ibm_surface_velocity', 'mod_ibm_coupled_targets',
         'mod_ibm_pressure_neumann', 'mod_movingIBM', 'mod_poiss_mgbicg',
         'lica_sgs', 'lica_movingIBM', 'lica_sphere', 'poiss_mgbicgstab',
         'ambient_boundaries', 'lica']
FLAGS = ['-std=legacy', '-fopenmp', '-fdefault-integer-8', '-fdefault-real-8',
         '-fdefault-double-8', '-ffree-line-length-none', '-O2']


def compiler_env(name):
    found = shutil.which(name)
    if not found:
        raise SystemExit('gfortran not found; specify --compiler with its full path')
    compiler = str(Path(found).resolve())
    env = os.environ.copy()
    env['PATH'] = str(Path(compiler).parent) + os.pathsep + env.get('PATH', '')
    env.update(OMP_NUM_THREADS='1', OMP_DYNAMIC='FALSE', OPENBLAS_NUM_THREADS='1')
    return compiler, env


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compiler', default='gfortran')
    parser.add_argument('--output', type=Path, default=ROOT / 'build_reference')
    args = parser.parse_args()
    compiler, env = compiler_env(args.compiler)
    output = args.output.resolve()
    output.mkdir(parents=True, exist_ok=False)
    commands = []
    with (output / 'build.log').open('w', encoding='utf-8') as log:
        for name in NAMES:
            shutil.copy2(ROOT / 'src' / (name + '.f90'), output)
            flags = FLAGS + (['-fmax-stack-var-size=1048576'] if name in
                            ['lica', 'lica_sgs', 'poiss_mgbicgstab'] else ['-frecursive'])
            cmd = [compiler, *flags, '-c', name + '.f90', '-o', name + '.o']
            commands.append(cmd)
            subprocess.run(cmd, cwd=output, env=env, stdout=log,
                           stderr=subprocess.STDOUT, check=True, timeout=180)
        executable = 'lica_reference.exe' if os.name == 'nt' else 'lica_reference'
        cmd = [compiler, '-fopenmp']
        if os.name == 'nt':
            cmd.append('-Wl,--stack,268435456')
        cmd += [name + '.o' for name in NAMES] + ['-o', executable]
        commands.append(cmd)
        subprocess.run(cmd, cwd=output, env=env, stdout=log,
                       stderr=subprocess.STDOUT, check=True, timeout=180)
    result = {'status': 'COMPILE_PASS', 'solver_launched': False,
              'compiler': subprocess.check_output([compiler, '--version'], env=env,
                                                 text=True).splitlines()[0],
              'commands': commands, 'executable': str(output / executable)}
    (output / 'result.json').write_text(json.dumps(result, indent=2), encoding='utf-8')
    print('COMPILE_PASS; no CFD solver was launched: ' + str(output / executable))


if __name__ == '__main__':
    main()
