# Interactive Mass–Spring–Damper Simulator

This project is an interactive simulation of a **mass–spring–damper system** implemented in **Julia** using **DifferentialEquations.jl** for numerical integration and **GLMakie.jl** for real-time visualization.

The simulator allows users to adjust system parameters and observe how the dynamic response changes for different damping conditions.

## System Description

The mass–spring–damper system is a second-order mechanical system consisting of:

* mass **m** (kg)
* spring stiffness **k** (N/m)
* viscous damping coefficient **b** (N·s/m)

The mass oscillates vertically under the restoring force of the spring and the dissipative force of the damper.

## Equation of Motion

Applying Newton’s Second Law:

$$
m\ddot{x}(t) + b\dot{x}(t) + kx(t) = 0
$$

where $x(t)$ is the displacement from the equilibrium position.

## Transfer Function

Taking the Laplace transform with zero initial conditions:

$$
G(s) = \frac{X(s)}{F(s)} = \frac{1}{ms^2 + bs + k}
$$

Dividing the denominator by $m$ gives:

$$
G(s) = \frac{1/m}{s^2 + \frac{b}{m}s + \frac{k}{m}}
$$

This expression can be written in the standard form:

$$
G(s) = \frac{\omega_n^2}{s^2 + 2\zeta\omega_n s + \omega_n^2}
$$

where

$$
\omega_n = \sqrt{\frac{k}{m}}
$$

is the natural frequency, and

$$
\zeta = \frac{b}{2\sqrt{mk}}
$$

is the damping ratio.

## Damping Behavior

| Case              | Condition   | Behavior                                          |
| ----------------- | ----------- | ------------------------------------------------- |
| Underdamped       | $\zeta < 1$ | Oscillatory response with decaying amplitude      |
| Critically damped | $\zeta = 1$ | Fastest return to equilibrium without oscillation |
| Overdamped        | $\zeta > 1$ | Slow return to equilibrium without oscillation    |

## Features

* Real-time mass–spring–damper animation
* Adjustable system parameters:

  * mass **m**
  * damping **b**
  * stiffness **k**
  * initial displacement **x₀**
* Interactive GUI controls
* Live displacement plot
* Visualization of equilibrium and initial position

## Simulation Tools

* **Julia**
* **DifferentialEquations.jl**
* **GLMakie.jl**

Solver used:

$$
\text{Tsit5 (adaptive Runge–Kutta 4/5)}
$$

## Installation

Install the required packages in Julia:

```julia
using Pkg
Pkg.add("DifferentialEquations")
Pkg.add("GLMakie")
```

## Running the Simulation

Run the Julia script:

```julia
julia simulator.jl
```

The interactive GUI will open where system parameters can be adjusted using sliders.

## Example Parameter Sets

**Underdamped**

```text
m = 1
k = 1
b = 0.5
```

**Critically damped**

```text
m = 1
k = 1
b = 2
```

**Overdamped**

```text
m = 1
k = 1
b = 4
```

## Educational Purpose

This simulator is designed to help visualize the dynamics of second-order systems and understand the effect of damping on system behavior.

It can be used as a teaching tool for courses in:

* Dynamics
* Control Systems
* Mechanical Vibrations
* System Modeling
