# Differential Equation Solver

A comprehensive C++ library for solving **Ordinary Differential Equations (ODEs)** and **Partial Differential Equations (PDEs)** with Python bindings. This project implements multiple numerical integration methods with error analysis and convergence testing capabilities.

## Features

- **Multiple ODE Stepping Methods:**
  - Euler Method (1st order)
  - Midpoint Method (2nd order)
  - Runge–Kutta 3rd Order (RK3)
  - Runge–Kutta 4th Order (RK4)
  - Runge–Kutta 5th Order (RK5)

- **PDE Support:**
  - 1D, 2D, and 3D Heat Equation solver using Method of Lines
  - Flexible boundary conditions (Dirichlet and Neumann)

- **Error Analysis:**
  - Multiple error metrics: L2 norm, RMS, relative error, absolute error, max norm
  - Convergence testing with automatic order estimation
  - Error tracking throughout integration

- **Performance & Usability:**
  - High-performance C++ backend
  - Complete Python bindings via pybind11
  - Fixed-step and adaptive integration
  - Detailed result logging and visualization utilities

## System Requirements

### macOS
- **macOS 10.15** or later
- **Clang/LLVM compiler** (included with Xcode)
- **CMake 3.20** or later
- **Python 3.8** or later
- **pybind11** (automatically found if installed via Homebrew or pip)

### Required Tools
```bash
brew install cmake python@3.11 pybind11
```

## Build Instructions

### 1. Clone and Navigate to Repository
```bash
git clone https://github.com/KRTRibb/Differential-Equation-Solver.git
cd Differential-Equation-Solver
```

### 2. Set Up Python Virtual Environment (Recommended)
```bash
# Create virtual environment
python3 -m venv .venv

# Activate it
source .venv/bin/activate

# Install required Python packages
pip install --upgrade pip
pip install pybind11 matplotlib numpy scipy
```

### 3. Create Build Directory
```bash
mkdir -p build
cd build
```

### 4. Configure with CMake
```bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DPython3_FIND_VIRTUALENV=ONLY \
      -DPython3_ROOT_DIR="$(python3 -c 'import sys; print(sys.prefix)')" \
      ..
```

For M1/M2 Macs, you may also specify the architecture:
```bash
cmake -DCMAKE_BUILD_TYPE=Release \
      -DCMAKE_OSX_ARCHITECTURES=arm64 \
      -DPython3_FIND_VIRTUALENV=ONLY \
      -DPython3_ROOT_DIR="$(python3 -c 'import sys; print(sys.prefix)')" \
      ..
```

### 5. Build the Project
```bash
cmake --build . --config Release -j$(sysctl -n hw.ncpu)
```

This command:
- Builds all targets (executables and Python module)
- Uses all available CPU cores for parallel compilation
- Produces optimized release builds

### 6. Verify Build Success
The following targets will be created in the `build/` directory:
- `diffeq_demo` – C++ demonstration program
- `1d_heat_example` – 1D heat equation example
- `diffeqpy.so` – Python extension module (installed to site-packages)

## Usage

### C++ Examples

#### Running the Demo
```bash
cd build
./diffeq_demo
```

#### Running Heat Equation Example
```bash
./1d_heat_example
```

### Python Usage

After building, activate your virtual environment and import the module:

```python
import diffeqpy
import numpy as np

# Define an ODE system: dy/dt = -y
def f(t, y):
    return [-y[0]]

# Initial condition
y0 = [1.0]

# Create problem
prob = diffeqpy.IVPProblem(f, y0, 0.0)

# Choose a stepper
stepper = diffeqpy.RK4Stepper()

# Create solver and integrate
solver = diffeqpy.Solver(stepper)
result = solver.integrateFixedSteps(prob, t_end=1.0, h=0.01)

# Access results
print(f"Integration steps: {result.n_steps}")
print(f"Computation time: {result.total_time} seconds")
```

### Python Examples

Several example scripts are provided in the `python/` directory:

1. **plotter.py** – Compares different stepper methods with exact solutions
2. **predator_prey.py** – Solves the Lotka–Volterra predator-prey system
3. **plot_animate.py** – Animates ODE solution trajectories
4. **compare_steppers.py** – Analyzes convergence behavior
5. **2nd_order.py** – Solves a damped harmonic oscillator

Run any example:
```bash
# From the root directory with venv activated
python3 python/predator_prey.py
python3 python/plot_animate.py
```

## Project Structure

```
.
├── CMakeLists.txt              # Build configuration
├── include/                    # C++ header files
│   ├── solver.hpp              # Main solver interface
│   ├── problem.hpp             # IVP/BVP problem definitions
│   ├── types.hpp               # Type definitions (Vec, RHS)
│   ├── Steppers/               # Stepper implementations
│   │   ├── stepper.hpp         # Base stepper class
│   │   ├── euler_stepper.hpp
│   │   ├── midpoint_stepper.hpp
│   │   ├── rk3_stepper.hpp
│   │   ├── rk4_stepper.hpp
│   │   └── rk5_stepper.hpp
│   ├── pde/                    # PDE solvers
│   │   └── heat_equation.hpp   # 1D/2D/3D heat equation
│   └── errorfunc/              # Error metrics
│       ├── metric.hpp
│       └── funcs.hpp
├── src/                        # C++ implementation files
│   ├── main.cpp
│   ├── solver.cpp
│   ├── bindings.cpp            # Python bindings
│   ├── Steppers/               # Stepper implementations
│   └── pde/                    # PDE implementations
├── examples/                   # C++ example programs
│   ├── heat_1d.cpp
│   └── 2nd_order.cpp
├── python/                     # Python scripts and utilities
│   ├── plotter.py
│   ├── predator_prey.py
│   ├── plot_animate.py
│   ├── compare_steppers.py
│   └── ThermalLab/             # Jupyter notebook examples
└── build/                      # Build output (created by CMake)
```

## Examples

### ODE: Exponential Decay
Solve $\frac{dy}{dt} = -y$ with $y(0) = 1$:

**C++:**
```cpp
RHS f = [](double t, const Vec& y) {
    Vec dydt;
    dydt.resize(1);
    dydt[0] = -y[0];
    return dydt;
};

Vec y0 = {1.0};
IVPProblem prob(f, y0, 0.0);

RK4Stepper rk4;
Solver solver(rk4);
auto result = solver.integrateFixedSteps(prob, 1.0, 0.01);
```

**Python:**
```python
def f(t, y):
    return [-y[0]]

y0 = [1.0]
prob = diffeqpy.IVPProblem(f, y0, 0.0)
stepper = diffeqpy.RK4Stepper()
solver = diffeqpy.Solver(stepper)
result = solver.integrateFixedSteps(prob, 1.0, 0.01)
```

### PDE: 1D Heat Equation
The 1D heat equation $u_t = \alpha u_{xx}$ on $[0,1]$ with Dirichlet boundary conditions is solved using the Method of Lines.

See [heat_1d.cpp](examples/heat_1d.cpp) for a complete implementation that solves the heat equation with initial condition $u(x,0) = \sin(\pi x)$ and zero boundary conditions.