# What this snapshot contains

| Item | Completed 500 rpm case |
|---|---|
| Method | Staggered nonuniform tensor-product Cartesian grid; moving IBM |
| Physics | Incompressible, isothermal; ν = 1.5×10⁻⁵ m²/s; initially quiescent |
| Motion | −500 rpm about +X; 24 ms cosine ramp; long duct inlet → short end outlet |
| Geometry | Original rotor/hub and rotating shaft; fixed casing; six experimental pins removed; handle excluded from CFD |
| Grid | 105×97×97 face-coordinate entries, 104×96×96 = 958,464 finite-volume cells |
| Resolution | Minimum spacing ≈0.363–0.369 mm; original radial gap ≈0.0697 mm, not resolved |
| Time | Explicit RK3 convection; Crank–Nicolson/ADI viscous treatment; pressure projection each stage |
| Pressure | Nonperiodic MG-BiCGStab; both relative 10⁻⁸ and scaled absolute 10⁻⁷ criteria |
| SGS input | Vreman (`ILES=1`, `ISGSMOD=3`) retained; not a claim of developed or resolved turbulence |
| Compute | 8 OpenMP threads; no MPI/GPU backend in this package |

The reference scales remain tied to 5,000 rpm, while `rotation_scale.in=0.1` gives actual 500 rpm. Do not change only the RPM label or confuse nondimensional time with seconds. `setup/case.json` records the conversion scales.

Ambient openings use pressure-projected normal velocity, permit reverse flow, set incoming tangential velocity to zero, and use zero normal derivatives for outgoing tangential velocity. The pressure face treatment includes an inflow kinetic correction (`p/ρ = −0.5 min(u_n,0)^2` in nondimensional form). The central exterior is solid; lateral openings are enabled in the end reservoirs. This is not six unrestricted velocity-zero faces.

## Wall/fluid interpolation repair

Within the scoped rotor region (x ∈ [−66.75,−59] mm, r ≤ 11.1 mm), unsafe interpolation stencils are replaced using the **actual wall intersection B** and a visible native **fluid donor F**:

`u_G = (1 + α) u_B − α u_F`, with `α = |G−B| / |B−F| ≤ 1`.

The actual wall velocity is Ω×r for a rotating owner and zero for a stationary owner. G is a ghost/forcing location, not the wall. The search uses 26 directions at one or two index offsets, rejects non-fluid donors and rays crossing a second solid, and stops if no valid support exists. It does not clamp the solved fluid speed to the blade-tip speed. Safe existing stencils and regions outside the repair scope retain the original treatment.

The coupled-target module separates predictor values and interpolated targets. The scoped pressure treatment initializes newly uncovered fluid pressure from visible current/previous fluid donors and reconstructs old-pressure gradients with a zero-normal-gradient constraint and a rank-adaptive tangential fit. Rank loss is counted; missing donors fail. The pressure correction operator is unchanged: this is **not** a globally rederived fluid-only Neumann Poisson solver.

Integration points: `lica_movingIBM.f90` → `rebuild_fluid_wall_stencils`; `mod_ibm_surface_velocity.f90` → wall intersections/donors; `mod_ibm_coupled_targets.f90` → predictor/targets; `mod_ibm_pressure_neumann.f90` → pressure handling. The synthetic kernel checks exercise the real velocity module, not a rewritten formula.

## Evidence and limits

The archived run reached 5.00166393 s at 500 rpm and exited normally at the requested checkpoint. Terminal fluid-cell peak ≈0.56585 m/s; tip speed ≈0.55896 m/s; outlet area-mean axial speed ≈0.007629 m/s and flux ≈4.650×10⁻⁶ m³/s. These are a single recorded section/field, not an uncertainty bound or product prediction. Blade-tip speed is not a universal upper bound on fluid speed.

The presentation intentionally uses a coarse mesh for tractable demonstration time. At least four cells across the tip gap were requested as a minimum planning target; mesh-dependence testing is separate. This archive does not meet that gap-resolution target. Wall stress, torque, tip leakage, statistical convergence and experimental similarity remain unqualified. Post-projection wall interpolation error is nonzero and recorded in the terminal summary.

The November 16–18 high-RPM convection/ghost/dissipation experiments are **not** included; they did not produce these videos. Preserve this version when comparing or migrating the completed case.
