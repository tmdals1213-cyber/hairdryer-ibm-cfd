# Public video handoff

The 37-second portrait and landscape movies are byte-identical to the approved November 16 red-blade editions. The short summary extracts the result sequence from the landscape movie (frame 991 onward) at 720p. No new CFD field, timing change or particle trajectory was invented for this release.

Maintain the black background, translucent generalized body, gray rotor and **one fixed material blade marked red**. The red identity rotates with the rotor. It is not a contour or a screen-space highlight, and periodic blade patterns may still appear to reverse perceptually. Geometry is not mirrored to correct that illusion.

The upper views show overall flow and a rotor close-up. Lower views show instantaneous 3D streamlines and separate speed/kinematic-pressure sections. Particle colors retain their fixed logarithmic speed scale. Stored fields start at 10.6 ms and end at 5.001664 s; no artificial zero-time field is inserted. The result sequence uses approximately 4× slow motion and 60 fps.

The opening Q-derived visual uses a fixed threshold 647.347 s⁻², cosmetic 0.30 mm smoothing and speed colors 0–0.6 m/s. It is a presentation surface, not a quantitative vortex boundary or proof of turbulence. Future quantitative postprocessing must retain native geometry/masks, fixed physical thresholds, nonuniform-grid gradients and explicit wall-stencil exclusion.

The packaged native field reader can compute Q, λ₂ and vorticity from saved velocity fields; no new solve is required for instantaneous diagnostics. Time-resolved particle accuracy is limited by the saved cadence. Derivative stencils touching solid cells must be masked using geometry at the matching rotor phase. Temporal RMS is not automatically a converged TKE statistic.

Recreating the exact 3D film requires the private authorized geometry, saved fields and original tracked-particle/render assets; these are deliberately absent. This package supplies the approved films and portable native-field utilities, not an asset-free recreation of the full editing pipeline. For new cases use the same visual conventions, preserve upstream particle identity through the rotor, and distinguish streamlines from time-integrated pathlines.
