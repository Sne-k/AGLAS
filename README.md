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

The open-loop baseline in this table is 145.2 kN m rather than the 142.1 kN m
above. The two are different models and both are correct. The analysis pipeline
applies the Kussner gust-penetration lag but models the wing in still air; the
control plant drops the Kussner lag and instead includes quasi-steady
aerodynamic feedback from the wing's own motion, which is essential for a
controller to be testable off its design condition. The 2 % difference between
them is the net of those two choices.

| Controller | Peak root moment | Reduction | Peak deflection |
|---|---|---|---|
| Open loop | 145.2 kN m | — | — |
| LQG | 68.4 kN m | 52.9 % | 17.9 deg |
| MRAC augmented | 71.4 kN m | 50.8 % | 17.3 deg |
| Policy search | 58.6 kN m | **59.6 %** | 20.0 deg |

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
eigenvector-based mode tracking. The control plant uses the quasi-steady limit
of the same model.

One consequence is worth stating because it changes how the wing should be
thought about: **aerodynamic damping dominates structural damping.** The first
bending damping ratio is 0.30 in flight, not the 0.02 structural value. The
wing is far better damped in the air than the structural model alone suggests,
which limits how much any controller can add.

**Control.** An 8-state aeroservoelastic plant, three mass-normalised modes
plus a second-order actuator. LQG with an observer that estimates the gust as
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
tests/                          verification suite, 76 checks
docs/control-design-log.md      full control design record, including dead ends
data/                           generated .mat results
results/                        generated figures
```

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
| Tip twist | `m*L^2/(2*GJ)` | 6e-15 relative |
| Root bending moment | `q*L^2/2` | 0.042 %, second-order convergent |
| Theodorsen `C(k)` | published tables, k = 0.1 to 1.0 | 8e-5 absolute |
| von Karman spectrum | its own -5/3 inertial slope | -1.6666 against -1.6667 |
| Divergence speed | `GJ*(pi/2L)^2/(c*e*cl_alpha)` | 0.026 % |
| Riccati solver | double integrator, `K = [1, sqrt(3)]` | exact |

The Simulink model is checked separately against the MATLAB reference by
`compare_simulink`, and agrees to **0.003 %** on peak root moment.

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
- **The control plant is quasi-steady.** It omits the `C(k)` phase lag and so
  loses stability at about 377 m/s where the unsteady solution puts flutter at
  314.8 m/s. Use it for design well below that boundary and `flutter_pk` for
  anything touching the stability margin.

## On learned control

LQG is exactly optimal, so there should be nothing to learn. It is optimal for
a different problem though: it minimises an integral quadratic cost on a linear
plant, while the structure is sized by the single largest moment and the real
loop has hard deflection and rate limits that LQR cannot represent.

Searching directly over four scales on the LQG gain, scored on the real
objective with saturation included, closes that gap and gives a 14 % further
reduction. It is not overfitted: tuned on one gust, it still wins by 11 to 21 %
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

The project is named for an adaptive system, so the finding that matters most
is this one: **on this aircraft, adaptation does not earn its complexity.**

MRAC is fully implemented, as a Lyapunov-based augmentation of the LQG with
normalisation, sigma-modification and a parameter bound. It was evaluated
against speed variation, mass variation, control effectiveness loss, discrete
gusts and sustained turbulence. In every case it matched or underperformed the
fixed-gain LQG it augments.

The reason is not a tuning failure. Fixed-gain LQG is stable from 35 to
200 m/s and performs *better* at 80 m/s than at its 50 m/s design point,
because higher dynamic pressure hands it more control authority. There is no
performance deficit for adaptation to recover. Reporting otherwise would
require choosing a flight condition outside the usable envelope, or weakening
the baseline being compared against.

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
