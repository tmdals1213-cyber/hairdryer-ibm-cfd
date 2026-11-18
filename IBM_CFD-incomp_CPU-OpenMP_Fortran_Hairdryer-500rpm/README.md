# IBM · incompressible CFD · CPU/OpenMP · Fortran · 500 rpm

## 1. Build and test (no private data needed)

Requirements: Python 3.10+ and gfortran with OpenMP. Windows GCC 14.2 was used; the helpers also contain a POSIX build path, which has not been tested here. Run commands from this case folder.

```console
python tools/verify_release.py
python tools/run_kernel_checks.py --compiler gfortran
python tools/build_reference.py --compiler gfortran
```

If needed, replace `gfortran` with the full executable path. Build folders must be new. The kernel test compiles the **actual unchanged velocity module** against synthetic walls; it tests interpolation, wall ownership, blocked/solid donors, and the expected missing-support stop. It does not run the hairdryer solver.

## 2. Prepare a private case

Obtain the authorized `fixed.fn` and original-gap `rotor.fn` files separately. Their required SHA-256 values are in `setup/case.json`. No CAD-to-function conversion or proprietary geometry download is hidden in this package.

```console
python tools/prepare_case.py --fixed /private/fixed.fn --rotor /private/rotor.fn --output runs/startup --steps 20
```

This validates geometry and creates inputs only. The default is a short **fresh quiescent-start setup**, 10 µs cap, CFL 0.1; it is not the exact archived restart trajectory. To launch deliberately:

```console
python tools/run_case.py --run runs/startup --exe build_reference/lica_reference.exe --threads 8
```

On Linux use `build_reference/lica_reference`. Keep the compiler's runtime directory on PATH, or pass `--runtime-dir`. The runner is foreground, uses a per-run lock, writes a receipt, never loops/retries, and will not overwrite a previous run. Monitor `run_progress.dat`, `pressure_stage_audit.dat`, `fluid_wall_stencils.dat` and `solver.log`. Create an empty `STOP_AFTER_STEP` file in the run folder to request a saved stop.

For a separately supplied accepted 500 rpm checkpoint at or beyond 0.14733 s:

```console
python tools/prepare_case.py --fixed /private/fixed.fn --rotor /private/rotor.fn --restart /private/fld014733 --profile continuation --output runs/continue --steps 20000 --end-time-s 5.002
```

The continuation profile retains the archived 200 µs cap, CFL 0.5, strict pressure criteria and save-every-10-steps choice. Header compatibility cannot establish that an arbitrary checkpoint is physically acceptable: inspect its provenance and diagnostics first. Exact historical replay also needs the original accepted checkpoints, which are excluded. No automatic RPM escalation or production supervisor is included.

## 3. Inspect results and share

`results/completed_500rpm.json` records the accepted terminal field. The videos are in `media/`; the summary was cut from the approved landscape result sequence. See `docs/VISUALIZATION.md` for the rendering conventions.

For local fields only, install `pip install -r post/requirements.txt`, then:

```console
python post/inspect_field.py --grid runs/startup/grid.dat --field runs/startup/fld000020
```

This reports units/time/finite values, not a fluid-only speed or force certificate. `post/lica_flow.py` also provides nonuniform-grid Q, λ₂ and vorticity utilities; a matching geometry mask is necessary before interpreting wall-adjacent derivatives.

Copy `docs/AGENT_PROMPT_KO.txt` into the receiving agent. For a different solver, port the algorithm and tests, not the case coordinates.
