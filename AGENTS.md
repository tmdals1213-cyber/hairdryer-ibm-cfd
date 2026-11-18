# Agent handoff

Work in `IBM_CFD-incomp_CPU-OpenMP_Fortran_Hairdryer-500rpm`. Read its README, `setup/case.json`, and `docs/AGENT_PROMPT_KO.txt` first.
This archive is the completed 500 rpm version. Do not silently substitute newer high-RPM code or label coarse results as validated LES.
Preserve Fortran source bytes and original copyright notices. Build into a new ignored directory.
No private CAD/STL/FN, raw fields, checkpoints, internal media, credentials or machine-specific paths may be committed. Use ignored `private/` and `runs/` directories locally.
Verify `MANIFEST.json` with `python tools/verify_release.py`; compile and run the synthetic kernel tests before porting. Never auto-start a long CFD run, change RPM, or overwrite a checkpoint based solely on this file. Use the user's current authorized scope.
The geometry function and repair region are case-specific. Port algorithms with explicit tests, not just coordinate constants. Missing valid fluid donors must be reported, not hidden by speed clipping.
