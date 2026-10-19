# AGLAS control design log

A record of how the gust load alleviation controller was arrived at, including
the things that did not work. Kept because the dead ends carry as much
information as the final design, and because a report written from the final
code alone would silently repeat several of the mistakes below.

Every number here came from an actual run. Nothing is estimated.

---

## 1. What had to be built before any control was possible

The open-loop model had a disturbance input and no control input. There was no
control surface anywhere in the project, so gust load alleviation was not
merely unimplemented, it was inexpressible.

**Control effector.** A trailing-edge flap, 25 % chord, covering 6.0 to 9.5 m
of the semi-span, snapped to element boundaries so the work-equivalent nodal
load vector stays exact. Thin-aerofoil theory with a plain flap gives

| Quantity | Value |
|---|---|
| dCl/ddelta | 3.826 /rad, 60.9 % of 2*pi |
| dCm_c4/ddelta | -0.6495 /rad |
| Lift per radian | 8789 N/m |
| Moment about the elastic axis per radian | -919.5 N m/m |

The 61 % figure is the textbook value for a quarter-chord flap. Note the net
moment about the elastic axis is **negative** even though the flap lift acts
ahead of that axis: the quarter-chord couple dominates. That is the mechanism
of aeroelastic effectiveness loss, and it is what eventually produces control
reversal (section 7).

**Aeroelastic feedback.** The plant originally had no dependence on airspeed
at all. Its A matrix was structure plus actuator, with aerodynamics appearing
only in the input matrices. That makes the poles independent of flight
condition, so no controller could be tested off its design point in any
meaningful way. Quasi-steady strip theory was added, taken as the k to 0 limit
of the same Theodorsen model the flutter analysis uses.

Consequence, and it is a large one: aerodynamic damping dominates structural
damping. First bending damping goes from the 0.02 structural value to **0.30**
at 50 m/s. A hand check of 0.5*rho*U*c*cl_alpha acting through the mode shape
gives 0.279, so the implementation is right. The wing is far better damped in
flight than the structural model alone suggests, and that materially limits how
much a controller can add.

**Plant summary.** 8 states: 3 mass-normalised modes (6) plus a second-order
15 Hz actuator (2). Two measurements, root strain and tip acceleration. Two
performance outputs, root bending moment and tip deflection. Fully
controllable and observable, though the Hautus margins show the flap's
authority is concentrated almost entirely in mode 1 (B_modal = [2154, 14, -128]);
mode 2 has a node inside the flap span, so the load largely self-cancels.

---

## 2. Toolbox independence

`lqr`, `care`, `lyap`, `ctrb` and `obsv` are all Control System Toolbox. The
rest of the project has no toolbox dependencies and that property was kept.

- `care_solve` solves the algebraic Riccati equation by the Hamiltonian
  eigenvector method. Verified against two closed forms: the scalar case
  A=0, B=1, Q=1, R=1 gives P=1, K=1 exactly; the double integrator gives
  K = [1, 1.7320508] against the exact [1, sqrt(3)].
- `lyap_solve` vectorises to a single n^2 linear solve. Residual 1e-16.
- `ctrb_obsv_rank` builds normalised Krylov matrices and also reports Hautus
  minimum singular values, which are better conditioned than a rank test when
  the plant mixes 2 Hz structural modes with a 15 Hz actuator.

---

## 3. Two bugs that produced convincing wrong answers

Both are worth recording because both produced output that looked plausible.

### 3.1 An observer with no gust model

The first LQG used a Kalman filter over the structural states only. The gust
appears directly in the tip accelerometer through the feedthrough term Dw, and
the filter had no way to represent it, so it explained that signal as
structural motion and produced a badly wrong estimate.

The symptom was unmistakable once looked for: **output feedback appeared to
beat full-state feedback**, 51.8 % load reduction against 7.1 %. That is
impossible. Full-state LQR is the optimum for this cost and nothing measuring
less can do better. The "better" result was a mis-estimating filter driving the
surface into its stops, where the resulting bang-bang action happened to
suppress the peak.

Fix: model the gust as a first-order Markov state and estimate it alongside the
structure. That also improves the regulator, because the Riccati solution then
produces a gain on the gust state, which is disturbance feedforward. For load
alleviation that term does most of the work; a purely reactive loop can never
get ahead of a disturbance arriving on the structure's own timescale.

### 3.2 An observer too fast to integrate

After the fix, LQG still beat the ideal, now by 97.3 % against 29.1 %.
Diagnosis by direct comparison of estimate against truth showed estimator
states reaching **1e300**. It was not a control result at all; it was numerical
overflow.

Cause: the process-to-measurement noise ratio placed the fastest observer pole
at **3340 Hz**. Classical RK4 is stable only while |lambda|*dt stays below
about 2.78, and at dt = 2e-4 s that product was 26. The integration diverged
while still reporting a plausible-looking peak.

| gust_psd | observer bandwidth | verdict at dt = 2e-4 s |
|---|---|---|
| 1e-4 | 10.5 Hz | too slow to be useful |
| 1e-2 | 68.5 Hz | usable |
| **0.05** | **150 Hz** | **chosen, 7x the highest mode** |
| 1 | 668 Hz | marginal |
| 25 | 3340 Hz | diverges |

Two guards were added so this cannot pass silently again: an explicit RK4
stability check against the fastest closed-loop pole, and a divergence test on
the states each step. Both raise an error naming the cause.

With those in place, LQG lands at 29.0 % against the ideal's 29.1 %, estimator
relative error 0.8 %. That ordering is what a correct implementation must show.

---

## 4. Preliminary tuning

Single knob, `r_command`, the price on commanded deflection. Open-loop peak
root moment 145.3 kN m.

| r_command | Peak moment | Reduction | Peak deflection | Saturated |
|---|---|---|---|---|
| 0.001 | 72.2 kN m | 50.3 % | 20.1 deg | 85.4 % |
| 0.01 | 64.2 kN m | 55.8 % | 20.1 deg | 14.5 % |
| 0.1 | 63.9 kN m | 56.0 % | 20.0 deg | 10.5 % |
| **0.3** | **68.4 kN m** | **52.9 %** | **17.9 deg** | **0 %** |
| 1 | 108.1 kN m | 25.6 % | 8.7 deg | 0 % |
| 10 | 140.4 kN m | 3.4 % | 1.1 deg | 0 % |

**r = 0.3 is the design point**: best reduction available without touching a
limit. The two cheaper-control rows beat it by about three points, but only by
riding the deflection stop for 10 to 15 % of the event, which leaves no margin
for a larger gust.

The r = 0.001 row is the instructive one. Asking for more authority than the
surface has stops helping: the loop sits on a limit 85 % of the time and does
worse than r = 0.3. LQR is optimal for the unconstrained linear problem and has
no representation of the actuator limits, so it happily designs a gain that
saturates, and a saturated loop is no longer the loop that was designed.

*Table re-measured after the actuator-limit alignment between Simulink and the
`.m` reference. An earlier version of this table, measured before that fix,
reported r = 0.001 as six times worse than no control (-536 %) and r = 0.3 as
55.6 %; neither reproduces with the current limiter.*

---

## 5. MRAC: implemented, and it does not earn its place here

Architecture: adaptive augmentation of the working LQG, not a replacement.
`u = -K*xhat + theta'*phi`, reference model is the nominal closed loop,
Lyapunov update with normalisation, sigma-modification and a norm bound.

### 5.1 Scaling

The first configuration used gamma = 5e3 and pinned theta at its bound within
the first gust, making the loop **far worse** than the LQG it augmented
(8.2 % reduction against 52.9 %). The drive term scales as
`Gamma*|phi|*|B'*P_lyap*e|`, and with Bryson-scaled weights that product is of
order |phi| ~ 5 and |P_lyap*B| ~ 0.76. Order 5 is the correct scale, three
orders of magnitude lower.

At the design point, where there is nothing to adapt to, any adaptation only
adds control activity, and the faster it adapts the worse the result, until
the parameter bound starts to limit the damage:

| gamma | Peak moment | vs LQG |
|---|---|---|
| 0.05 | 68.8 kN m | +0.6 % |
| 0.5 | 71.8 kN m | +5.0 % |
| 5 | 135.3 kN m | +97.8 % |
| 50 | 118.1 kN m | +72.6 % |

### 5.2 Off-nominal, which is the case adaptation exists for

Controller designed at 50 m/s and nominal mass, then flown elsewhere:

| Condition | Open loop | LQG | MRAC | MRAC benefit |
|---|---|---|---|---|
| Nominal 50 m/s | 145.3 | 68.4 | 71.8 | -5.0 % |
| 35 m/s | 101.2 | 65.3 | 66.2 | -1.5 % |
| 80 m/s | 236.4 | 61.0 | 178.6 | -193 % |
| Fuel burnt, -11 % mass | 144.8 | 68.3 | 72.0 | -5.4 % |
| 80 m/s and fuel burnt | 235.8 | 61.0 | 223.1 | -266 % |

Control effectiveness loss, controller still believing it has full authority:

| Lambda | LQG | MRAC | MRAC benefit |
|---|---|---|---|
| 1.00 | 68.4 | 71.8 | -5.0 % |
| 0.60 | 95.0 | 93.9 | +1.1 % |
| 0.40 | 112.7 | 112.3 | +0.3 % |
| 0.15 | 133.0 | 133.1 | -0.0 % |

**MRAC does not help on this aircraft.** Two independent reasons:

1. **There is no performance deficit to recover.** Fixed-gain LQG is stable
   from 35 to 248 m/s and only loses stability at 250 m/s. At 80 m/s it
   performs *better* than at its design point, 74 % reduction against 53 %,
   because higher dynamic pressure hands it more control authority. The
   uncertainty makes the plant easier to control, not harder.

2. **A discrete gust gives adaptation no time to converge.** The load peak
   arrives about one second in, comparable to the structure's own settling
   time. Under effectiveness loss the parameter needed is of order
   `K*(1-Lambda)/Lambda ~ 0.15`; theta only reached 0.047 before the event was
   over. Adaptation needs persistent excitation over many time constants and a
   one-shot gust does not provide it.

This is a result, not a failure to tune. Reporting MRAC as beneficial here
would require either choosing a flight condition outside the usable envelope or
quietly weakening the LQG baseline it is compared against.

### 5.3 Sustained turbulence: the open question, now answered

The obvious objection to the discrete-gust result is that adaptation never had
time to work. Under sustained von Karman turbulence the excitation is
persistent and the parameters have tens of seconds to converge, which is the
condition MRAC theory actually assumes. Thirty seconds at 2.10 m/s RMS, root
bending moment RMS in kN m:

| Lambda | Open loop | LQG | Best MRAC | MRAC benefit |
|---|---|---|---|---|
| 1.00 | 30.02 | 14.11 | 14.53 | -3.0 % |
| 0.40 | 30.02 | 20.67 | 21.13 | -2.2 % |

Best MRAC is over gamma = 0.05, 0.5 and 5. **MRAC is worse here too**,
including under a 60 % control effectiveness loss, though by a few percent
rather than the 14 to 23 % an earlier measurement recorded with a different
set of adaptation rates. So the
timescale argument in section 5.2 was not the whole story: giving adaptation
fifteen times longer does not rescue it.

The likelier explanation is that the adaptive law is driven by the error
between the plant and a reference model that only sees the gust through the
*estimated* gust state. Under broadband turbulence that error is dominated by
disturbance the reference model cannot reproduce, not by model mismatch. The
adaptation then chases disturbance-induced error, sigma-modification spends its
effort pulling the parameters back, and the net result is extra control
activity with no uncertainty being cancelled.

For the record the LQG itself performs consistently across both disturbance
types: 53 % RMS reduction under turbulence against 53 % peak reduction under
the discrete gust.

---

## 6. Where fixed-gain LQG actually fails

Designed at 50 m/s, applied across the envelope:

| Airspeed | Closed loop | max Re(eig), closed loop | max Re(eig), open loop |
|---|---|---|---|
| 35 m/s | stable | -3.03 | -2.87 |
| 80 m/s | stable | -1.92 | -4.81 |
| 200 m/s | stable | -3.63 | -1.63 |
| 240 m/s | stable | -3.13 | **+3.25** |
| 250 m/s | **unstable** | > 0 | > 0 |
| 300 m/s | **unstable** | +24.86 | +13.31 |

The failure boundary sits at 250 m/s, well above any realistic operating speed
for this wing. Between 216 and 248 m/s the LQG is stabilising a plant that is
itself unstable, but that plant is quasi-steady and its own boundary at 216 m/s
is an artefact of the missing circulation lag (section 8), so none of this
region should be trusted.

---

## 7. Control reversal, an incidental finding

Static control effectiveness, steady root moment per radian of flap:

| Airspeed | dM/ddelta |
|---|---|
| 50 m/s | +239,705 N m/rad |
| 240 m/s | +3,518,211 N m/rad (peak) |
| 280 m/s | +2,802,846 N m/rad |
| **320 m/s** | **-1,656,622 N m/rad, reversed** |

Effectiveness peaks near 240 m/s then collapses through zero at about 300 m/s.
Beyond that the surface works backwards: the nose-down twist it induces
overwhelms the lift it generates. This is genuine aeroelastic control reversal,
it falls out of the model rather than being put in, and it is another reason
nothing above 200 m/s should be trusted to a fixed-sign control law. Standard
MRAC cannot rescue it either, since the Lyapunov argument assumes the sign of
the control effectiveness is known.

---

## 8. Quasi-steady versus unsteady, a caution

The default control plant is quasi-steady and therefore omits the C(k) phase
lag that drives classical flutter. It loses stability at about **216 m/s** in a
17.7 Hz mode, whereas the unsteady p-k solution puts flutter at **314.8 m/s and
13.1 Hz**. Use the control plant for design well below that boundary and
`flutter_pk` for anything touching the stability margin itself.

An earlier version of this section reported the quasi-steady boundary as
377 m/s and called the over-prediction expected. That number came from a
defect: the plant's aerodynamic damping was taken as Im Q(k)/k at k = 1e-5,
but Im C(k) ~ k ln k near zero, so the division left a ln(1e-5) term that made
the torsion mode about 27 times better damped than C = 1 strip theory gives.
The quasi-steady plant now uses the exact C = 1 limit. At the 50 m/s design
point the fix moved the results in sections 4 and 5 by less than 0.5 %; the
tables above are re-measured with it.

`cfg.aero.plant_model = 'unsteady'` realises C(k) in the time domain with a
two-state indicial lag, adds apparent mass and the Kussner gust lag, and its
stability boundary reproduces p-k to 1 %: 312 m/s and 13.21 Hz. With
Prandtl-Glauert compressibility the boundary falls to 264 m/s, Mach 0.78.

---

## 9. Simulink model

Generated by script rather than committed as a binary, so it can be reviewed in
a diff and regenerated when the design matrices change. Fixed-step ode4 at
2e-4 s, because the rate limiter is a stateful nonlinearity that a
reject-and-retry solver can corrupt.

Written without a Simulink licence, so all block dimensions were first checked
only against the design matrices. A block diagram can be wrong in ways that
raise no error: a swapped Mux order, a Selector off by one, a sign on a Sum.
`compare_simulink.m` runs the model against the verified `.m` path on
identical inputs and reports the difference.

**Since executed.** With MRAC on at the nominal condition, the model matches
the `.m` reference to 0.003 % on peak root moment (71.81 kN m in both) and on
peak command (17.21 deg), 0.007 % RMS over the waveform and 0.008 % worst
pointwise. Rerun `compare_simulink` after any change to the design matrices.

---

## 10. Status of the original three-part request

| Task | Status |
|---|---|
| Selection of control algorithm | Done. LQG with a gust-estimating observer, evidence in sections 3 to 5. |
| Algorithm implementation | Done, toolbox-free, in `src/common`. |
| Preliminary tuning | Done, section 4. Design point r = 0.3, 52.9 % load reduction. |
| MRAC | Implemented and evaluated. Does not help on this aircraft; section 5. |
| Simulink model | Built and verified against the `.m` reference to 0.003 %; section 9. |
| Reinforcement learning | Done as direct policy search; section 11. |
| Adaptive feedforward | Implemented and evaluated; helps off-design with gust preview; section 12. |
| Lidar preview | Swept 0 to 200 m; about 0.1 s of warning suffices on this wing; section 13. |

---

## 11. Learned control: direct policy search

### 11.1 Why there is anything to learn

LQG is exactly optimal, so at first sight there is nothing to search for. It is
optimal for a specific problem though, and that problem is not the one being
solved:

1. It minimises an **integral quadratic cost**. The structure is sized by the
   **single largest** root bending moment in the event. Different objectives,
   different optima.
2. It is derived for a **linear** plant. The real loop has a deflection limit
   and a slew limit, and LQR has no representation of either. Section 4 shows
   how badly that bites.

Deep reinforcement learning would be the wrong tool. The plant is known, low
order, and linear apart from two saturations; learning a value function from
scratch would spend orders of magnitude more evaluations rediscovering a
controller the Riccati equation already gives in closed form. The useful move
is to start from that controller and search only over what it cannot model.

### 11.2 Method

Four multiplicative scales on blocks of the LQG gain, starting at `[1 1 1 1]`,
which is the LQG solution itself, so the search can never do worse than its
starting point. Compass search: try plus and minus one step per coordinate,
move to the best improvement, halve the step when none improves. Deterministic,
derivative-free, no toolbox. Each evaluation is a full simulation of the
saturated closed loop, scored on peak root bending moment.

### 11.3 Results

| r_command | LQG | Tuned | Improvement | LQG saturation | Scales |
|---|---|---|---|---|---|
| 0.3 | 68.42 | **58.65** | +14.3 % | 0 % | [1.8, 0, 0.64, 1] |
| 0.03 | 64.40 | 59.26 | +8.0 % | 13.5 % | [0, 0.79, 0, 2.2] |
| 0.01 | 64.17 | 59.26 | +7.7 % | 14.5 % | [0, 0.85, 0.14, 1.9] |

Best is 58.65 kN m, a **59.6 %** reduction from the 145.3 kN m open-loop peak,
against LQG's 52.9 %. Policy search beats LQG at every control weight tried, so
the gain is not simply a matter of choosing a better `r`.

The tuned solution at the design weight is structurally interesting: it scales
modal position feedback **up** by 1.8 and turns modal rate feedback **off**
entirely. That is not a fit to one waveform, it is a different control law.

### 11.4 Does it generalise, or is it overfitted to one gust?

Tuned on the 2.0 s design gust, then evaluated on disturbances it never saw:

| Validation case | LQG | Tuned | Change |
|---|---|---|---|
| Gust duration 1.0 s | 70.29 | 60.84 | +13.4 % |
| Gust duration 3.0 s | 68.07 | 58.77 | +13.7 % |
| Gust duration 4.0 s | 67.95 | 58.92 | +13.3 % |
| Gust amplitude 6 m/s | 41.05 | 34.50 | +16.0 % |
| Turbulence, peak | 35.49 | 32.16 | +9.4 % |
| Turbulence, RMS | 14.11 | 12.20 | +13.5 % |
| 80 m/s, off design | 60.99 | 48.33 | +20.8 % |

The turbulence rows use the 30 s, 2.10 m/s RMS record of section 5.3.

It generalises across disturbance shape, amplitude, disturbance *class* and
flight condition. It even wins on RMS under turbulence, which is closer to the
integral criterion LQG optimises than to the peak criterion the search used.

### 11.5 The honest caveat

The improvement is not free. At the design point:

| | Peak moment | Peak deflection | RMS command | Position saturation |
|---|---|---|---|---|
| LQG | 68.42 kN m | 17.93 deg | 6.32 deg | 0 % |
| Tuned | 58.65 kN m | 20.01 deg | 7.22 deg | 3.9 % |

The tuned policy uses **14 % more control effort** and sits right on the
20 degree deflection limit, where LQG left two degrees unused. Part of the gain
is simply spending authority the LQG's `R` weight was holding back.

That is still a real finding rather than an artefact. `R` is an indirect proxy:
there is no way to tell an LQR "use all the authority available and no more",
and section 4 shows that trying to force it there by lowering `R` gains about
three points at best, spends 10 to 15 % of the event on the stop, and then
gets worse again as `R` falls further, because the linear design has no idea
the limit exists. Policy search beats every one of those rows.
Policy search optimises the real objective under the real constraint and finds
the boundary directly. The cost is that the tuned design operates with no
deflection margin, which is a trade a designer should make deliberately rather
than inherit.

---

## 12. Adaptive feedforward from a gust probe

Section 5 tested adaptation in the feedback path and found nothing to gain. The
gust load alleviation literature puts it somewhere else: Wildschek et al.
(AIAA 2006-6054), Zeng et al. (JGCD 33(3), 2010) and Wang, Li and Da Ronch
(AIAA 2015-2243) adapt a *feedforward* filter driven by a gust measurement taken
ahead of the wing. That is the case tested here.

### 12.1 Set-up

A probe 10 m ahead of the wing gives 0.2 s of warning at 50 m/s and 0.125 s at
80 m/s. An FIR filter of 100 taps at 200 Hz maps the probe history to a surface
command added to the LQG command ahead of the actuator limits. The filter is
either fixed, at the Wiener solution trained on a separate 120 s turbulence
record at 50 m/s, or adapted online by normalised filtered-x LMS
(mu = 0.003, optionally leakage 1e-4 per update) on measured root strain. The
secondary-path model is always the nominal 50 m/s closed loop, never updated.

Each case is 240 s of von Karman turbulence, a different realisation from the
training record, with a 10 m/s, 2 s one-minus-cosine gust arriving at 232 s.
RMS is taken over 120 to 232 s, after adaptation has had two minutes, and the
gust peak over the 6 s from 232 s with adaptation still running.

### 12.2 Results

| Condition | Controller | Turbulence RMS | vs LQG | Gust peak | Saturation |
|---|---|---|---|---|---|
| 50 m/s | LQG | 20.55 kN m | | 62.9 kN m | 0.0 % |
| 50 m/s | LQG + fixed FF | 11.17 kN m | -46 % | 49.5 kN m | 20.7 % |
| 50 m/s | LQG + adaptive FF | 9.79 kN m | -52 % | 52.0 kN m | 7.1 % |
| 50 m/s | LQG + adaptive FF, leak | 8.81 kN m | -57 % | 52.3 kN m | 6.3 % |
| 80 m/s | LQG | 19.12 kN m | | 59.8 kN m | 0.1 % |
| 80 m/s | LQG + fixed FF | 14.55 kN m | -24 % | 46.0 kN m | 18.2 % |
| 80 m/s | LQG + fixed FF, q-scheduled | 9.83 kN m | -49 % | 39.6 kN m | 1.4 % |
| 80 m/s | LQG + FF redesigned at 80 m/s | 8.88 kN m | -54 % | 27.4 kN m | 16.5 % |
| 80 m/s | LQG + adaptive FF | 6.86 kN m | -64 % | 23.8 kN m | 1.4 % |
| 80 m/s | LQG + adaptive FF, leak | 2.81 kN m | **-85 %** | **7.4 kN m** | 0.9 % |
| Effectiveness 0.4 | LQG | 30.42 kN m | | 99.1 kN m | 1.9 % |
| Effectiveness 0.4 | best feedforward option | 24.45 kN m | -20 % | 103.2 kN m | 47.2 % |

### 12.3 Is it the adaptation?

Three alternative explanations were tested and ruled out as the whole story.

1. **The baseline is merely conservative.** r = 0.3 was chosen to keep the
   discrete gust off the stops, not for turbulence. A fixed LQG at r = 0.01
   gives -60 % at 50 m/s, which matches adaptation there, so at the design
   point adaptation adds nothing. At 80 m/s the same LQG gives -70 % with 30 %
   saturation, against -85 % with 0.9 % for adaptation.
2. **The fixed filter just has the wrong gain off-design.** Mostly true:
   scaling the nominal filter by 0.5 at 80 m/s reproduces the redesigned
   filter exactly (8.18 kN m on the shorter record), and scheduling it on
   dynamic pressure gives -49 %. Adaptation still beats both.
3. **It is feedback in disguise, not feedforward.** Partly true. With leakage
   the weights end near zero (|w| = 0.007 against 0.57 for the Wiener filter),
   so the benefit does not come from a converged filter: each normalised LMS
   step moves the command by roughly mu*e/|S|, a high-gain, error-driven
   correction on root strain. But the probe matters. Moving it to the wing
   gives -71 % at 80 m/s, and replacing it with a constant reference gives
   -59 %.

### 12.4 Sensor noise

With `cfg.sensor.simulate_noise`, white noise at the estimator's design levels
(2e-6 strain, 0.02 m/s^2 acceleration) on the measurements, and at five times
those levels, every 80 m/s result above holds to within 0.1 kN m. Root strain
at these loads is about 3e-4, so even 5 times the design noise is a few percent
of the signal. With probe noise of 0.5 m/s RMS the adaptive result stays at
-85 % (gust peak 6.0 kN m) while the q-scheduled fixed filter degrades from
-49 % to -40 %.

### 12.5 Robustness to delay, unmodelled modes and gust coherence

Same set-up as 12.1, on a 180 s record (adaptation over 0 to 90 s, RMS over 90
to 172 s, gust at 172 s). `cfg.control.delay` delays the command reaching the
plant while the estimator assumes it acted at once. Unmodelled modes are
simulated by running a 6- or 8-mode wing against the 3-mode design. Coherence
a mixes the probe signal as a*w + sqrt(1-a^2)*w', with w' an independent
turbulence realisation, supplied through `ff.ref`.

Turbulence RMS against the r = 0.3 LQG, and the gust peak:

| Case | 50 m/s LQG r=0.01 | 50 m/s adaptive FF | 80 m/s LQG r=0.01 | 80 m/s q-scheduled FF | 80 m/s adaptive FF |
|---|---|---|---|---|---|
| Nominal | -61 % | -58 % | -71 %, 16.5 kN m | -49 %, 36.6 kN m | **-86 %, 6.0 kN m** |
| Delay 5 ms | **+32 %** | -58 % | -66 % | -48 % | -84 % |
| Delay 10 ms | +72 % | -58 % | -17 %, 53.1 kN m | -47 % | -83 %, 8.1 kN m |
| Delay 20 ms | +127 % | -57 % | +228 % | -44 % | -83 %, 8.2 kN m |
| Delay 40 ms | +205 % | -57 % | +378 % | -34 % | -65 %, 19.2 kN m |
| 6-mode wing | -57 % | -54 % | -63 % | -47 % | -77 %, 10.9 kN m |
| 8-mode wing | -56 % | -53 % | -62 % | -47 % | -76 %, 12.9 kN m |
| Coherence 0.9 | | -51 % | | -39 %, 29.9 kN m | -69 %, 32.5 kN m |
| Coherence 0.7 | | -28 % | | -24 %, 28.5 kN m | -41 %, 40.4 kN m |
| Coherence 0.5 | | -26 % | | -11 %, 30.4 kN m | -17 %, 34.6 kN m |

The r = 0.3 LQG and the q-scheduled filter move by a few percent at most under
delay; the aggressive LQG that matched adaptation at the design point is the
one that cannot tolerate it, and becomes worse than the baseline at 5 ms.

**Delay** therefore strengthens the case for adaptive feedforward: with any
realistic latency it beats every fixed alternative at both speeds.
**Unmodelled modes** cost it about ten points at 80 m/s without changing the
ranking. **Coherence** is the real limit. It still has the lowest turbulence
RMS of the probe-based options at every coherence tried, but below 1 its
discrete-gust peak is no better than the fixed filter's, and a zero-delay
aggressive LQG beats it.

The coherence model is pessimistic for a 10 m probe: it removes correlation
uniformly across all frequencies, whereas over that distance the gust field
decorrelates mainly at high frequency and mainly across the span. A
frequency-dependent spanwise coherence model would be the next refinement.

### 12.6 What this does and does not show

It shows that, in this simulation, adaptation earns its place off the design
point when it sits in a feedforward path with gust preview, and that the
advantage survives realistic sensor and probe noise. It overturns the
unqualified version of section 5's conclusion; section 5 stands for MRAC.

It does not show that the adaptive loop is stable in general: stability rests
on simulation over the cases tried, including those in 12.5. Of the model
errors tested, imperfect probe-to-wing coherence is the one that erodes the
margin.

---

## 13. How much gust preview is worth: a lidar look-ahead sweep

Khalil and Fezans (AIAA 2019-0822) and Fezans et al. (CEAS Aeronaut. J., 2019)
design load alleviation around a Doppler lidar measuring the gust 30 to 200 m
ahead. Section 12 used a 10 m probe. This section sweeps the look-ahead
distance to see how much further preview is worth on this wing.

### 13.1 Set-up

As in 12.5: 180 s of turbulence with a gust at 172 s, LQG at r = 0.3, fixed
(Wiener) and adaptive feedforward, the FIR filter sized to span the preview
plus 0.5 s (up to 900 taps). `probe_signal` produces the measurement: the gust
advanced by d/U, mixed frequency by frequency with an independent field
through the coherence gamma(f) = exp(-C*f*d/U). C = 0 is Taylor's frozen
turbulence; C = 1 and 4 let the gust evolve between sensor and wing. At 80 m/s
the fixed filter is trained at 80 m/s; an earlier run trained it at 50 m/s,
which put its timing wrong by d*(1/50 - 1/80), 1.5 s at 200 m, and is not
reported.

### 13.2 Results

Turbulence RMS against LQG (20.93 kN m at 50 m/s, 19.46 kN m at 80 m/s):

| Look-ahead | 0 m | 10 m | 25 m | 50 m | 100 m | 200 m |
|---|---|---|---|---|---|---|
| Preview at 50 m/s | 0 | 0.20 s | 0.50 s | 1.0 s | 2.0 s | 4.0 s |
| 50 m/s fixed, C = 0 | -22 % | -47 % | -47 % | -47 % | -46 % | -46 % |
| 50 m/s adaptive, C = 0 | -55 % | -57 % | -57 % | -56 % | -55 % | -53 % |
| 50 m/s fixed, C = 1 | -22 % | -42 % | -39 % | -36 % | -31 % | -25 % |
| 50 m/s adaptive, C = 1 | -55 % | -54 % | -51 % | -49 % | -43 % | -43 % |
| 50 m/s fixed, C = 4 | -22 % | -37 % | -30 % | -25 % | -19 % | -14 % |
| 50 m/s adaptive, C = 4 | -55 % | -49 % | -40 % | -37 % | -39 % | -36 % |
| 80 m/s fixed, C = 0 | -12 % | -55 % | -56 % | -55 % | -54 % | -54 % |
| 80 m/s adaptive, C = 0 | -73 % | **-85 %** | -80 % | -74 % | -77 % | -70 % |
| 80 m/s fixed, C = 1 | -12 % | -52 % | -46 % | -40 % | -35 % | -29 % |
| 80 m/s adaptive, C = 1 | -73 % | -77 % | -69 % | -61 % | -60 % | -44 % |
| 80 m/s adaptive, C = 4 | -73 % | -68 % | -56 % | -12 % | -6 % | -51 % |

The adaptive gust peak at 80 m/s and C = 0 rises from 6.9 kN m at 10 m to
21.9 kN m at 200 m. At C = 4 the adaptive filter becomes erratic at 50 and
100 m, spending over 20 % of the time on the stop.

**Even with a perfect, frozen-turbulence measurement, preview beyond about
10 m buys nothing**, and as soon as the gust evolves in transit every further
metre costs performance. The longer filters the adaptive scheme needs also
learn worse. On this wing the best sensor is the closest one that still gives
about 0.1 s of warning.

### 13.3 Why preview runs out so early

Three checks at 50 m/s, fixed feedforward, C = 0:

| Look-ahead | 0 m | 5 m | 10 m | 25 m | 50 m | 100 m |
|---|---|---|---|---|---|---|
| Baseline (from 13.2) | -22 % | | -47 % | -47 % | -47 % | -46 % |
| Actuator rate cut to 30 deg/s | -17 % | | -36 % | -36 % | -36 % | -36 % |
| Unsteady plant (Wagner, Kussner lags) | -47 % | -54 % | -54 % | -54 % | -54 % | -54 % |
| No actuator limits at all | -99 % | -100 % | -99 % | -99 % | -99 % | -100 % |

The first guess, that a fast actuator needs little warning, was wrong: a five
times slower actuator still levels off at 10 m. The last row gives the actual
reason. With unlimited deflection, a perfect gust measurement taken *at* the
wing cancels 99 % of the load, because flap and gust act on the same strips
through the same dynamics. The binding constraint is deflection authority (the
surface sits on its 20 degree stop about 20 % of the time), and warning time
cannot create authority. Preview only needs to cover the actuator and
aerodynamic lags, about 0.1 s, which the unsteady plant confirms with the
benefit complete by 5 m.

This is partly a property of the model. The wing is clamped, so there is no
rigid-body pitch or plunge, and the published lidar benefits come largely
from whole-aircraft response, slower spoiler-type effectors and rate-limited
manoeuvre loads. A free-flying model would be needed to see them.
