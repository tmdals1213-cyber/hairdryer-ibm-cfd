# Source notices and provenance

`src/lica.f90` retains the original 2017 LICA notice, including:

> Copyright(c) 2017 All rights reserved by Haecheon Choi
>
> Noncommercial use with copyrighted mark

The full contributor list and contact information remain in that source file. Original notices in other source files, including the John Burkardt attribution in `lica_sphere.f90`, are also retained. The 15 Fortran sources are byte-for-byte copies of the completed 500 rpm build inputs.

This repository does not replace those notices with MIT/GPL, grant additional upstream rights, or assert an unrestricted open-source license. Obtain any permissions needed for uses outside the supplied notice from the relevant rights holders. Public availability is not a new license grant. No third-party compiler or runtime binary is redistributed.

Project-specific work added the composite moving/fixed geometry interface, quiescent pressure openings, rotor-local wall/fluid interpolation and pressure-donor handling, diagnostics, packaging helpers and presentation. See `docs/NUMERICS.md` and `results/source_provenance.json` in the case folder for the exact scope and source hashes. These case-specific modifications are not represented as an upstream LICA release.

Public media were prepared from the completed numerical demonstration. Non-blade product geometry is generalized; private CAD and internal-report assets are excluded. These are research demonstration results, not a product endorsement or performance certification.
