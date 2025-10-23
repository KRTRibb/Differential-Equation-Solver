import numpy as np
import diffeqpy
from compare_steppers import StepperComparison


"First test for Circular Motion"
# def f(t, y):
#     return [y[1], -y[0]]

# def exact(t):
#     return [np.cos(t), -np.sin(t)]

# y0 = np.array([1.0, 0.0])
# t0 = 0.0
# t_end = 20
# h = 0.01
# error_metric = diffeqpy.L2

"Second test for Coupled Exponential Decay"
def f(t, y):
    return [-y[0] + 0.5*y[1], -2*y[1]]

def exact(t):
    return [np.exp(-t) + 0.5 * 2 * (np.exp(-t) - np.exp(-2 * t)), 2 * np.exp(-2 * t)]


y0 = np.array([1.0, 2.0])
t0 = 0.0
t_end = 1
h = 0.1
error_metric = diffeqpy.L2


comp = StepperComparison(f, exact, y0, t0, t_end, h, error_metric)

comp.add_stepper(diffeqpy.EulerStepper())
comp.add_stepper(diffeqpy.MidpointStepper())
comp.add_stepper(diffeqpy.RK3Stepper())
comp.add_stepper(diffeqpy.RK4Stepper())
comp.add_stepper(diffeqpy.RK5Stepper())

comp.plot_results()