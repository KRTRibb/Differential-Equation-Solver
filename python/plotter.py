import numpy as np
import diffeqpy
from compare_steppers import StepperComparison

def f(t, y):
    return [y[1], -y[0]]

def exact(t):
    return [np.cos(t), -np.sin(t)]

y0 = np.array([1.0, 0.0])
t0 = 0.0
t_end = 2 * np.pi
h = 0.01
error_metric = diffeqpy.L2

comp = StepperComparison(f, exact, y0, t0, t_end, h, error_metric)

comp.add_stepper(diffeqpy.EulerStepper())
comp.add_stepper(diffeqpy.MidpointStepper())
comp.add_stepper(diffeqpy.RK3Stepper())

comp.plot_results()