import sympy
from sympy import symbols, exp, sin, sqrt, integrate, Piecewise
from numpy import vectorize, linspace, pi
import matplotlib.pyplot as plt

print("Calculating symbolic integrals, please wait...")

m = 3
w_n = 16
zeta = 0.06

w_d = w_n * sqrt(1 - zeta**2)

t, tau = symbols('t tau', real=True, positive=True)

g = 1 / (m * w_d) * exp(-zeta * w_n * t) * sin(w_d * t)

g_tau = (1 / (m * w_d)) * exp(-zeta * w_n * (t - tau)) * sin(w_d * (t - tau))

F_1 = 5
F_2 = 11
w_1 = 13
w_2 = 20

T_1 = 2 * pi / w_1
T_2 = 2 * pi / w_2

t_vals = linspace(0, 5 * T_1, 2000)

f_1 = -F_1 * sin(w_1 * tau)
f_2 = F_2 * sin(w_2 * (tau - T_1))

def f_t(t):
    if 0 <= t <= T_1/2:
        return -F_1 * sin(w_1 * t)
    elif T_1 <= t <= T_1 + T_2/2:
        return F_2 * sin(w_2 * (t - T_1))
    else:
        return 0

f_vec = vectorize(f_t)

def solve_convolution(f_func, lower_limit, upper_limit):
    integrand = f_func * g_tau
    result = integrate(integrand, (tau, lower_limit, upper_limit))
    return result

x_1 = solve_convolution(f_1, 0, t)

x_2 = solve_convolution(f_1, 0, T_1/2)

x_3 = solve_convolution(f_2, T_1, t)

x_4 = solve_convolution(f_2, T_1, T_1 + T_2/2)

x_total = Piecewise(
    (x_1, (t >= 0) & (t <= T_1/2)),
    (x_2, (t > T_1/2) & (t < T_1)),
    (x_2 + x_3, (t >= T_1) & (t <= T_1 + T_2/2)),
    (x_2 + x_4, True)
)

print("Converting symbolic function to numerical function...")

x_total_func = sympy.lambdify(t, x_total, 'numpy')
x_vals = x_total_func(t_vals)

plt.figure(figsize=(10, 4))
plt.plot(t_vals, f_vec(t_vals))
plt.title("Excitation Force f(t)")
plt.xlabel("Time (s)")
plt.ylabel("Force (N)")
plt.grid()

plt.figure(figsize=(10, 4))
plt.plot(t_vals, x_vals)
plt.title("System Response x(t)")
plt.xlabel("Time (s)")
plt.ylabel("Displacement (m)")
plt.grid()

plt.axvline(T_1/2, color='k', linestyle='--', alpha=0.3)
plt.axvline(T_1, color='k', linestyle='--', alpha=0.3)
plt.axvline(T_1 + T_2/2, color='k', linestyle='--', alpha=0.3)

plt.show()