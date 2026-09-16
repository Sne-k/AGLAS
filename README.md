# AGLAS

Gust load response analysis of a flexible aircraft wing, in MATLAB.

A slender wing is idealised as a clamped-root cantilever beam and discretised
with finite elements. The model computes its free-vibration characteristics,
its transient response to atmospheric gusts, the resulting bending stress and
strain fields, the factor of safety, the wall thickness needed to meet a
safety target, and the aeroelastic stability boundary.

Everything here is plain MATLAB. There are no toolbox dependencies: the
unsteady aerodynamics needs only `besselj` and `bessely`, and the turbulence
model is synthesised directly rather than through a filter design function.
The code also runs unmodified under GNU Octave.

## Quick start

```matlab
run_aglas                    % verification suite, then all 12 stages
run_aglas('main')            % main pipeline only        (c1p ... c6p)
run_aglas('supplementary')   % supplementary only        (c1m ... c6m)
run_aglas('tests')           % verification suite only
run_aglas('both', false)     % skip the verification suite
```

A full run takes about two minutes and writes `.mat` results to `data/`,
figures to `results/` and an animated GIF to `animation/`.

To change the wing, edit `src/common/aglas_config.m`. Nothing else hard-codes
a physical parameter.

## Results

Baseline wing, 10 m semi-span, 1.5 m chord, aluminium 2024-T3 box spar, under
a 10 m/s one-minus-cosine gust of 2 s duration at 50 m/s.

| Quantity | Value |
|---|---|
| First bending frequency | 2.100 Hz |
| Second bending frequency | 13.149 Hz |
| First torsion frequency | 21.859 Hz |
| Peak tip deflection | 0.666 m, 6.7 % of span |
| Peak tip twist | 0.222 deg |
| Dynamic amplification factor | 1.041 |
| Peak root bending moment | 142.1 kN m |
| Peak root stress | 169.4 MPa, 49 % of yield |
| Minimum factor of safety | 2.037 |
| Wall thickness for a factor of safety of 2.5 | 11.42 mm, from 9.00 mm |
| Flutter speed | 314.8 m/s at 13.15 Hz |
| Divergence speed | 399.8 m/s |

Figures are in `results/`: mode shapes, gust response, V-g and V-f flutter
diagrams, stress and strain fields, safety zones, mesh convergence and the
thickness sizing study.

## The model

**Structure.** Euler-Bernoulli beam elements with cubic Hermite interpolation
and consistent mass matrices, 20 elements by default. Each node carries
vertical deflection, bending rotation and torsional twist. Bending and torsion
are elastically uncoupled, because twist is measured about the elastic axis,
but are inertially coupled through the offset between the elastic axis and the
section centre of gravity. That coupling is what makes bending-torsion flutter
possible; a bending-only beam cannot flutter at any speed.

**Section.** Thin-walled closed rectangular box spar. Bending inertia from the
standard hollow-rectangle formula, torsion constant from the Bredt-Batho
single-cell relation. A solid rectangular section is also available.

**Modal reduction.** The generalised eigenproblem is solved and the shapes are
normalised so that `Phi' * M * Phi = I`. Mass normalisation is what makes a
modal damping matrix `diag(2*zeta*omega)` mean 2 % of critical damping, and
makes `diag(omega^2)` the correct modal stiffness.

**Gusts.** The discrete gust is the CS-25.341 one-minus-cosine shape over a
full cosine cycle. Continuous turbulence is synthesised from the von Karman
vertical-velocity spectrum by random-phase Fourier synthesis, under a fixed
seed, so the realisation is reproducible and its spectrum is correct by
construction rather than approximated by a filter.

**Aerodynamics.** Strip theory. Gust loads use the Kussner indicial function
for gust penetration, evaluated by an exact recurrence rather than numerical
quadrature. Flutter uses Theodorsen's `C(k)` computed from Hankel functions,
solved by the p-k method with the aerodynamic matrix split into stiffness and
damping parts.

**Stress recovery.** Curvature is the exact second derivative of the element
interpolation, using both the translation and the rotation degrees of freedom,
so the clamped root, where a cantilever carries its largest moment, is
evaluated directly.

## Repository layout

```
run_aglas.m                     master runner for every stage
src/common/                     shared, verified physics library
src/main_pipeline/    c1p..c6p  wing, structure, gust response,
                                flutter, animation, stress
src/supplementary_pipeline/
                      c1m..c6m  gust models, mesh convergence, response
                                comparison, stress field, strain, sizing
tests/                          verification suite, 76 checks
data/                           generated .mat results (not tracked)
results/                        generated figures
animation/                      generated GIF (not tracked)
```

The supplementary pipeline reads the structural model and gust response from
the main pipeline, so run the main pipeline first or use `run_aglas('both')`.

## Verification

```matlab
addpath('tests'); run_all_tests
```

76 checks, all against closed-form references rather than previously recorded
output. A suite that only compared against stored numbers would have locked in
the defects this code was written to remove.

| Check | Reference | Agreement |
|---|---|---|
| Cantilever bending frequencies | Euler-Bernoulli characteristic equation | 5e-6 % on mode 1 |
| Static tip deflection | `q*L^4/(8*EI)` | 4e-12 relative |
| Tip slope | `q*L^3/(6*EI)` | 4e-12 relative |
| Tip twist | `m*L^2/(2*GJ)` | 6e-15 relative |
| Root bending moment | `q*L^2/2` | 0.042 %, second-order convergent |
| Torsion constant | Bredt-Batho | exact |
| Theodorsen `C(k)` | published tables, k = 0.1 to 1.0 | 8e-5 absolute |
| von Karman spectrum | its own -5/3 inertial slope | -1.6666 against -1.6667 |
| Quasi-steady gust load | `0.5*rho*U*c*cl_alpha*w_g` | exact |
| Divergence speed | `GJ*(pi/2L)^2/(c*e*cl_alpha)` | 0.026 % |
| Assembled mass | `m*L` | machine precision |

The suite also checks the properties whose absence caused the original
defects: that the total applied load does not change when the mesh is refined,
that peak stress lands at the root, that the response scales linearly with
gust amplitude, that the measured decay rate recovers the prescribed damping
ratio, that the gust returns to zero with zero slope, and that a bending-only
model reports no flutter rather than a silent `NaN`.

`run_all_tests` returns the failure count, so it can gate a build:

```sh
matlab -batch "addpath('tests'); exit(run_all_tests() > 0)"
```

## Known limitations

These are properties of the model, not defects. The code warns at run time
when a result falls outside the range where its assumptions hold.

- **Incompressible aerodynamics.** Theodorsen theory assumes incompressible
  flow and is reliable below about Mach 0.3. The computed flutter speed is
  Mach 0.93 at sea level, so it is indicative rather than a certifiable
  number. A compressible unsteady method, or at minimum a Prandtl-Glauert
  correction, would be needed to quote it.
- **Linear aerodynamics.** A 10 m/s gust at 50 m/s is an 11.3 degree change of
  incidence, beyond the roughly 10 degrees over which lift stays proportional
  to incidence. The computed loads are therefore an upper bound; a real
  aerofoil would begin to stall.
- **Strip theory.** Each spanwise station is treated as a two-dimensional
  section. There is no induced-flow correction and no tip relief, so loads
  near the tip are overestimated.
- **Linear beam theory.** Deflections are assumed small. The code warns above
  10 % of span; the baseline case reaches 6.7 %.
- **Prismatic, unswept, untapered.** No sweep, taper, twist or spanwise
  property variation.
- **Open loop.** The name AGLAS originally stood for Adaptive Gust Load
  Alleviation System. No control system is implemented; this is the
  uncontrolled structural response only.

## References

- Theodorsen, T. *General Theory of Aerodynamic Instability and the Mechanism
  of Flutter.* NACA Report 496, 1935.
- Bisplinghoff, R.L., Ashley, H. and Halfman, R.L. *Aeroelasticity.*
  Addison-Wesley, 1955.
- Hodges, D.H. and Pierce, G.A. *Introduction to Structural Dynamics and
  Aeroelasticity.* Cambridge University Press, 2nd ed., 2011.
- Wright, J.R. and Cooper, J.E. *Introduction to Aircraft Aeroelasticity and
  Loads.* Wiley, 2nd ed., 2014.
- EASA CS-25.341, discrete gust and continuous turbulence design criteria.
- MIL-HDBK-1797, von Karman and Dryden atmospheric turbulence models.
