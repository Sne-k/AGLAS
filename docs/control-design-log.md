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
root moment 145.2 kN m.

| r_command | Peak moment | Reduction | Peak deflection | Saturated |
|---|---|---|---|---|
| 0.001 | 931.1 kN m | **-536 %** | 20.2 deg | 90.6 % |
| 0.01 | 73.1 kN m | 50.1 % | 20.1 deg | 16.6 % |
| 0.1 | 69.3 kN m | 52.6 % | 20.0 deg | 8.6 % |
| **0.3** | **64.9 kN m** | **55.6 %** | **18.1 deg** | **0 %** |
| 1 | 104.0 kN m | 29.0 % | 9.3 deg | 0 % |
| 10 | 138.5 kN m | 5.4 % | 1.6 deg | 0 % |

**r = 0.3 is the design point**: best reduction available without touching a
limit.

The r = 0.001 row is the most instructive in the table. Asking for more
authority than the surface has does not merely stop helping, it makes the load
**six times worse than doing nothing**. LQR is optimal for the unconstrained
linear problem and has no representation of the actuator limits, so it happily
designs a gain that saturates, and a saturated loop is no longer the loop that
was designed. Any claim that an optimal controller is safe by construction
should be read against this row.

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

At the design point, where there is nothing to adapt to, the degradation is a
clean function of adaptation rate, which is exactly what theory predicts:

| gamma | max abs theta | Peak moment | vs LQG |
|---|---|---|---|
| 0.05 | 0.0012 | 68.7 kN m | +0.5 % |
| 0.5 | 0.0101 | 71.4 kN m | +4.5 % |
| 5 | 0.0273 | 87.8 kN m | +28.4 % |
| 50 | 0.0478 | 101.0 kN m | +47.8 % |

### 5.2 Off-nominal, which is the case adaptation exists for

Controller designed at 50 m/s and nominal mass, then flown elsewhere:

| Condition | Open loop | LQG | MRAC | MRAC benefit |
|---|---|---|---|---|
| Nominal 50 m/s | 145.2 | 68.4 | 71.4 | -4.5 % |
| 35 m/s | 101.2 | 65.2 | 66.0 | -1.1 % |
| 80 m/s | 236.1 | 60.9 | 70.0 | -15.0 % |
| Fuel burnt, -11 % mass | 144.7 | 68.3 | 71.7 | -5.0 % |
| 80 m/s and fuel burnt | 235.4 | 60.8 | 213.8 | -251 % |

Control effectiveness loss, controller still believing it has full authority:

| Lambda | LQG | MRAC | MRAC benefit |
|---|---|---|---|
| 1.00 | 68.4 | 71.5 | -4.5 % |
| 0.60 | 95.0 | 94.0 | +1.0 % |
| 0.40 | 112.7 | 112.3 | +0.3 % |
| 0.15 | 132.9 | 132.2 | +0.5 % |

**MRAC does not help on this aircraft.** Two independent reasons:

1. **There is no performance deficit to recover.** Fixed-gain LQG is stable
   from 35 to 200 m/s and only loses stability at 240 m/s. At 80 m/s it
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

**Open question, not yet answered.** Whether MRAC pays off under *sustained*
von Karman turbulence, where the excitation is persistent and the parameters
have tens of seconds to converge, is the obvious follow-up. A 40 s run at
Lambda = 0.4 was started and killed before it finished. It remains untested and
nothing should be claimed about it.

---

## 6. Where fixed-gain LQG actually fails

Designed at 50 m/s, applied across the envelope:

| Airspeed | Closed loop | max Re(eig) |
|---|---|---|
| 35 to 200 m/s | stable | about -3.6 |
| 240 m/s | **unstable** | +15.8 |
| 300 m/s | **unstable** | +47.1 |

The failure boundary sits between 200 and 240 m/s, well above any realistic
operating speed for this wing and approaching the 314.8 m/s flutter speed where
the quasi-steady plant is no longer trustworthy anyway.

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

The control plant is quasi-steady and therefore omits the C(k) phase lag that
drives classical flutter. It loses stability at about **377 m/s** in a 6.5 Hz
oscillatory mode, whereas the unsteady p-k solution puts flutter at
**314.8 m/s and 13.1 Hz**. Quasi-steady theory over-predicting flutter speed is
the expected result, not a discrepancy to reconcile. Use the control plant for
design well below that boundary and `flutter_pk` for anything touching the
stability margin itself.

---

## 9. Simulink model

Generated by script rather than committed as a binary, so it can be reviewed in
a diff and regenerated when the design matrices change. Fixed-step ode4 at
2e-4 s, because the rate limiter is a stateful nonlinearity that a
reject-and-retry solver can corrupt.

**Not executed.** No Simulink licence was available where it was written. All
block dimensions were checked against the design matrices, but a block diagram
can be wrong in ways that raise no error: a swapped Mux order, a Selector off
by one, a sign on a Sum. `compare_simulink.m` runs the model against the
verified `.m` path on identical inputs and reports the difference. Run it before
trusting any model output.

---

## 10. Status of the original three-part request

| Task | Status |
|---|---|
| Selection of control algorithm | Done. LQG with a gust-estimating observer, evidence in sections 3 to 5. |
| Algorithm implementation | Done, toolbox-free, in `src/common`. |
| Preliminary tuning | Done, section 4. Design point r = 0.3, 55.6 % load reduction. |
| MRAC | Implemented and evaluated. Does not help on this aircraft; section 5. |
| Simulink model | Built, untested here, awaiting a run. |
| Reinforcement learning | Not started. |

On reinforcement learning: the honest framing is that LQG is already the exact
optimum for the linear quadratic problem, so a learned policy has nothing to
win on the nominal plant. The one place a learned policy could genuinely beat
it is the regime section 4 exposes, where saturation is active and the linear
design is no longer optimal. That argues for direct policy search over a
structured controller against the saturated closed loop, rather than deep
reinforcement learning, and it should be presented as such rather than as
reinforcement learning for its own sake.
