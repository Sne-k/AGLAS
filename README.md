# AGLAS

Adaptive Gust Load Alleviation System. Aeroelastic analysis and gust load
alleviation control for a flexible aircraft wing, in MATLAB.

A slender wing is idealised as a clamped-root cantilever beam and discretised
with finite elements. The project computes its free-vibration characteristics,
its transient response to atmospheric gusts, the resulting stress and strain
fields, the factor of safety, the wall thickness needed to meet a safety
target, and the aeroelastic stability boundary. It then designs, tunes and
simulates a closed-loop controller that cuts the peak gust load in half.

Everything is plain MATLAB with no toolbox dependencies. The unsteady
aerodynamics needs only `besselj` and `bessely`; the turbulence model is
synthesised directly rather than through a filter design function; and the
Riccati, Lyapunov and controllability routines are written out rather than
taken from the Control System Toolbox. The analysis code also runs unmodified
under GNU Octave. Only the Simulink model needs Simulink.

## Quick start

Open-loop analysis, twelve stages plus the verification suite:

```matlab
run_aglas                    % verification suite, then all 12 stages
run_aglas('main')            % main pipeline only        (c1p ... c6p)
run_aglas('supplementary')   % supplementary only        (c1m ... c6m)
run_aglas('tests')           % verification suite only
```

`run_aglas('main')` and `run_aglas('supplementary')` still run the full
verification suite first and abort if any check fails; pass `false` as a
second argument to skip it, e.g. `run_aglas('main', false)`. The
supplementary pipeline reads the structural model and gust response written
by the main pipeline (`c3m`, `c4m`), so `run_aglas('supplementary')` on its
own only works after `run_aglas('main')` has been run at least once — for a
single self-contained run, use `run_aglas` or `run_aglas('both')`.

Closed-loop control, in Simulink:

```matlab
cd simulink
aglas_sim_setup      % design the controller
build_gla_model      % generate and open aglas_gla.slx
compare_simulink     % run it and check against the MATLAB reference
```

A full open-loop run takes about two minutes and writes `.mat` results to
`data/`, figures to `results/` and an animated GIF to `animation/`.

To change the wing, edit `src/common/aglas_config.m`. Nothing else hard-codes
a physical parameter.

## Results

Baseline wing: 10 m semi-span, 1.5 m chord, aluminium 2024-T3 box spar, under a
10 m/s one-minus-cosine gust of 2 s duration at 50 m/s.

### Structure and loads

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

### Gust load alleviation

A 25 % chord trailing-edge surface over 6.0 to 9.5 m of the semi-span, with a
20 degree deflection limit and a 150 deg/s rate limit.

The open-loop baseline in this table is 145.3 kN m rather than the 142.1 kN m
above. The analysis pipeline applies the Kussner gust-penetration lag but
models the wing in still air; the default control plant drops the Kussner lag
and instead includes quasi-steady aerodynamic feedback from the wing's own
motion, which is essential for a controller to be testable off its design
condition. The 2 % difference between them is the net of those two choices.

| Controller | Peak root moment | Reduction | Peak deflection |
|---|---|---|---|
| Open loop | 145.3 kN m | — | — |
| LQG | 68.4 kN m | 52.9 % | 17.9 deg |
| MRAC augmented | 71.8 kN m | 50.6 % | 17.2 deg |
| Policy search | 58.6 kN m | **59.6 %** | 20.0 deg |

With the unsteady control plant (`cfg.aero.plant_model = 'unsteady'`, see
below), which carries the Kussner lag as well as the wing's own circulation lag,
the open-loop peak is 141.2 kN m, within 0.6 % of the analysis pipeline, and
the ranking is unchanged:

| Controller | Unsteady, incompressible | Unsteady, compressible |
|---|---|---|
| Open loop | 141.2 kN m | 142.6 kN m |
| LQG | 69.4 kN m, 50.8 % | 69.4 kN m, 51.3 % |
| MRAC augmented | 71.7 kN m, 49.2 % | 71.8 kN m, 49.6 % |
| Policy search | 59.2 kN m, **58.1 %** | 59.8 kN m, **58.0 %** |

Under sustained von Karman turbulence the LQG cuts the root moment RMS from
30.0 to 14.1 kN m, a 53 % reduction, consistent with the discrete gust result.

The policy-search controller is the best of the three, but it buys that with
14 % more control effort and sits on the deflection limit where LQG leaves two
degrees unused. See the section on learned control below.

Figures are in `results/`.

## The model

**Structure.** Euler-Bernoulli beam elements with cubic Hermite interpolation
and consistent mass matrices, 20 elements by default. Each node carries
vertical deflection, bending rotation and torsional twist. Bending and torsion
are elastically uncoupled, because twist is measured about the elastic axis,
but are inertially coupled through the offset between that axis and the section
centre of gravity. That coupling is what makes bending-torsion flutter
possible; a bending-only beam cannot flutter at any speed.

**Aerodynamics.** Strip theory throughout. Gust loads use the Kussner indicial
function for gust penetration, evaluated by an exact recurrence. Flutter uses
Theodorsen's `C(k)` from Hankel functions, solved by the p-k method with
eigenvector-based mode tracking. The control plant uses the same strip theory
in the time domain, in one of two forms. The default is the quasi-steady limit,
C(k) = 1. The unsteady form realises C(k) with a two-state indicial lag (R.T.
Jones' fit to the Wagner function, or the compressible fit of Leishman and
Nguyen, AIAA J. 28(5), 1990), adds the apparent-mass terms and the Kussner gust
lag, and optionally applies Prandtl-Glauert compressibility. Its stability
boundary reproduces the p-k flutter result to 1 % (312 m/s and 13.21 Hz,
against 314.8 m/s and 13.15 Hz); with compressibility the boundary falls to
264 m/s.

One consequence is worth stating because it changes how the wing should be
thought about: **aerodynamic damping dominates structural damping.** The first
bending damping ratio is 0.30 in flight, not the 0.02 structural value. The
wing is far better damped in the air than the structural model alone suggests,
which limits how much any controller can add.

**Control.** An 8-state aeroservoelastic plant, three mass-normalised modes
plus a second-order actuator, or 18 states with the unsteady aerodynamic lags. LQG with an observer that estimates the gust as
a state alongside the structure, which is what turns the law from a purely
reactive loop into one with disturbance feedforward.

## Repository layout

```
run_aglas.m                     master runner for the open-loop pipeline
src/common/                     shared, verified library: structures,
                                aerodynamics, flutter, and control
src/main_pipeline/    c1p..c6p  wing, structure, gust response,
                                flutter, animation, stress
src/supplementary_pipeline/
                      c1m..c6m  gust models, mesh convergence, response
                                comparison, stress field, strain, sizing
simulink/                       closed-loop model, built by script
tests/                          verification suite, 102 checks
docs/control-design-log.md      full control design record, including dead ends
data/                           generated .mat results
results/                        generated figures
```

## Verification

```matlab
addpath('tests'); run_all_tests
```

The same suite runs on GitHub Actions for every push to `main` and every pull
request (`.github/workflows/tests.yml`).

102 checks. The 76 physics checks are all against closed-form references rather
than previously recorded output; 7 check the unsteady control plant against
Theodorsen's `C(k)` and the p-k flutter boundary; the 14 control checks cover the Riccati,
Lyapunov and controllability routines against closed forms, and pin the
closed-loop results in the table above so a change to the controller cannot
silently move them; 5 cover the gust feedforward path. A suite that only compared against stored numbers would have locked in
the defects this code was written to remove.

| Check | Reference | Agreement |
|---|---|---|
| Cantilever bending frequencies | Euler-Bernoulli characteristic equation | 5e-6 % on mode 1 |
| Static tip deflection | `q*L^4/(8*EI)` | 4e-12 relative |
| Tip twist | `m*L^2/(2*GJ)` | 6e-15 relative |
| Root bending moment | `q*L^2/2` | 0.042 %, second-order convergent |
| Theodorsen `C(k)` | published tables, k = 0.1 to 1.0 | 8e-5 absolute |
| von Karman spectrum | its own -5/3 inertial slope | -1.6666 against -1.6667 |
| Divergence speed | `GJ*(pi/2L)^2/(c*e*cl_alpha)` | 0.026 % |
| Riccati solver | double integrator, `K = [1, sqrt(3)]` | residual 3e-15 |
| Lyapunov solver | scalar `A = -1, Q = 2`, `P = 1` | exact |
| LQG / MRAC / policy search | peak root moment reductions in the table above | reproduced to 0.1 % |

The Simulink model is checked separately against the MATLAB reference by
`compare_simulink`, and agrees to **0.003 %** on peak root moment and peak
command, and 0.007 % RMS over the waveform.

## Known limitations

These are properties of the model, not defects. The code warns at run time
when a result falls outside the range where its assumptions hold.

- **Incompressible aerodynamics.** Theodorsen theory is reliable below about
  Mach 0.3. The computed flutter speed is Mach 0.93 at sea level, so it is
  indicative rather than certifiable.
- **Linear aerodynamics.** A 10 m/s gust at 50 m/s is an 11.3 degree change of
  incidence, beyond the roughly 10 degrees over which lift stays proportional
  to incidence. Computed loads are an upper bound.
- **Strip theory.** No induced-flow correction and no tip relief, so loads near
  the tip are overestimated.
- **Linear beam theory.** Deflections assumed small. The code warns above 10 %
  of span; the baseline reaches 6.7 %.
- **Prismatic, unswept, untapered.** No sweep, taper or spanwise variation.
- **The default control plant is quasi-steady.** It omits the `C(k)` phase lag
  and loses stability at about 216 m/s, well below the 314.8 m/s flutter speed.
  Use it for design well below that boundary. The unsteady plant reproduces the
  p-k boundary to 1 %, but the design-log experiments were measured with the
  quasi-steady one.
- **Compressibility is Prandtl-Glauert only.** The unsteady plant's
  compressible option scales lift slopes and lag rates with
  `beta = sqrt(1 - M^2)`, a linear subsonic correction that refuses to run
  above Mach 0.9. Its 264 m/s boundary is Mach 0.78, where transonic effects
  already matter; the apparent-mass terms stay incompressible.

## On learned control

LQG is exactly optimal, so there should be nothing to learn. It is optimal for
a different problem though: it minimises an integral quadratic cost on a linear
plant, while the structure is sized by the single largest moment and the real
loop has hard deflection and rate limits that LQR cannot represent.

Searching directly over four scales on the LQG gain, scored on the real
objective with saturation included, closes that gap and gives a 14 % further
reduction. It is not overfitted: tuned on one gust, it still wins by 9 to 21 %
on other gust durations, other amplitudes, continuous turbulence and an
off-design flight speed.

This is direct policy search, the derivative-free end of reinforcement
learning, and it is labelled as such rather than as deep RL. Deep RL would be
the wrong tool: the plant is known, low order and linear apart from two
saturations, so learning from scratch would spend orders of magnitude more
evaluations rediscovering what a Riccati equation gives in closed form. The
useful move is to start from that solution and search only over what it cannot
model.

## On the "adaptive" in AGLAS

The project is named for an adaptive system, and the answer turns out to depend
on where the adaptation sits.

**Adaptive feedback (MRAC) does not earn its complexity.** MRAC is fully
implemented, as a Lyapunov-based augmentation of the LQG with normalisation,
sigma-modification and a parameter bound. It was evaluated against speed
variation, mass variation, control effectiveness loss, discrete gusts and
sustained turbulence. In every case it matched or underperformed the fixed-gain
LQG it augments. That is not a tuning failure: fixed-gain LQG is stable from
35 to 248 m/s and performs *better* at 80 m/s than at its 50 m/s design point,
because higher dynamic pressure hands it more control authority.

**Adaptive feedforward from a gust probe does, off the design point.** Following
Wildschek et al. (AIAA 2006-6054) and Zeng et al. (JGCD 2010), a probe 10 m
ahead of the wing feeds an FIR filter that adapts online by filtered-x LMS
(`feedforward_design`, `ctrl.ff` in `closed_loop_sim`). In 240 s of
von Karman turbulence with a 10 m/s gust arriving at 232 s, adaptation running:

| Condition | Controller | Turbulence RMS | Gust peak | Saturation |
|---|---|---|---|---|
| 50 m/s | LQG, r = 0.3 | 20.55 kN m | 62.9 kN m | 0.0 % |
| 50 m/s | LQG, r = 0.01 (fixed, aggressive) | -60 % | 52.0 kN m | 6.8 % |
| 50 m/s | LQG + adaptive feedforward | -57 % | 52.3 kN m | 6.3 % |
| 80 m/s | LQG, r = 0.3 | 19.12 kN m | 59.8 kN m | 0.1 % |
| 80 m/s | LQG, r = 0.01 (fixed, aggressive) | -70 % | 16.7 kN m | 29.9 % |
| 80 m/s | LQG + fixed feedforward, scheduled on q | -49 % | 39.6 kN m | 1.4 % |
| 80 m/s | LQG + adaptive feedforward | **-85 %** | **7.4 kN m** | 0.9 % |

At the design point adaptation only matches a more aggressively tuned fixed
LQG. At 80 m/s it beats every fixed alternative tried, including a feedforward
filter redesigned for that speed, while barely touching the deflection stop.
Both ingredients are needed: with the probe moved to the wing it falls to
-71 %, and with a constant reference in place of the probe to -59 %. Under
60 % control effectiveness loss nothing helps, since every option is
saturation-bound. The result holds with sensor noise (up to 5 times the
estimator's design levels) and with 0.5 m/s RMS of probe noise.

Two cautions. First, the mechanism is not a better-learned filter: with
leakage the weights end near zero, and the normalised LMS step behaves as a
high-gain, error-driven correction on root strain. Its stability is shown only
by simulation, not proved. Second, the simulation has a perfect plant apart
from the off-design changes, no sensor delay, and an ideal probe position; any
of those could erode the margin.

`docs/control-design-log.md` has the full evidence, every measurement and every
dead end.

## References

- Theodorsen, T. *General Theory of Aerodynamic Instability and the Mechanism
  of Flutter.* NACA Report 496, 1935.
- Bisplinghoff, R.L., Ashley, H. and Halfman, R.L. *Aeroelasticity.*
  Addison-Wesley, 1955.
- Hodges, D.H. and Pierce, G.A. *Introduction to Structural Dynamics and
  Aeroelasticity.* Cambridge University Press, 2nd ed., 2011.
- Wright, J.R. and Cooper, J.E. *Introduction to Aircraft Aeroelasticity and
  Loads.* Wiley, 2nd ed., 2014.
- Lavretsky, E. and Wise, K.A. *Robust and Adaptive Control.* Springer, 2013.
- EASA CS-25.341, discrete gust and continuous turbulence design criteria.
- MIL-HDBK-1797, von Karman and Dryden atmospheric turbulence models.
