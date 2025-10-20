# Differential Equation Solver

This project is a simple yet efficient differential equation solver designed to handle **Initial Value Problems (IVPs)** using various **Runge–Kutta methods** up to third order.

### Features
- **Supported Methods:**
  - `EulerStepper` (1st-order Runge–Kutta)
  - `MidpointStepper` (2nd-order Runge–Kutta)
  - `RK3Stepper` (3rd-order Runge–Kutta)
- **Convergence Testing:** Built-in tools to analyze step-size convergence for each method.
- **C++ Backend:** Core numerical routines implemented in C++ for performance.
- **Python Bindings:** Python interface for easier experimentation, data analysis, and plotting.
- **Basic Plotting Utilities:** Includes simple visualizations to compare steppers and their convergence behavior.

### Future Improvements
- Support for **Boundary Value Problems (BVPs)**
- Additional **Runge–Kutta methods** (e.g., RK4)
- Expanded **plotting and data visualization** tools
- Potential **finite element method (FEM)** extensions for PDEs

---