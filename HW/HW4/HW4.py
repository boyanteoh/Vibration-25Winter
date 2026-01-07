import matplotlib.pyplot as plt
import sympy as sp
import numpy as np

print("This code requires a few minutes to run! Please be patient...")

t, tau = sp.symbols('t, tau', real=True)

m = 3.0
wn = 16.0
zeta = 0.06
F1 = 5.0
w1 = 13.0
F3 = 11.0
w2 = 20.0

wd = wn * np.sqrt(1 - zeta ** 2)
T1 = 2 * np.pi / w1
T2 = 2 * np.pi / w2

F2 = F3 * np.exp(-1j * w2 * T1)

A1 = -1j / (2 * m * wd)
s1 = -zeta * wn + 1j * wd
s2 = -zeta * wn - 1j * wd

function_g = A1 * sp.exp(s1 * (t - tau)) + sp.conjugate(A1) * sp.exp(s2 * (t - tau))

f1 = (1j * F1 / 2) * sp.exp(1j * w1 * tau) - (1j * F1 / 2) * sp.exp(-1j * w1 * tau)
f2 = (-1j * F2 / 2) * sp.exp(1j * w2 * tau) + (1j * np.conjugate(F2) / 2) * sp.exp(-1j * w2 * tau)

print("Integrating x0 (Phase 1)...")
x0_sym = sp.integrate(f1 * function_g, (tau, 0, t))

print("Integrating x1 (Pulse 1 Memory)...")
x1_sym = sp.integrate(f1 * function_g, (tau, 0, T1 / 2))

print("Integrating x2 (Active Pulse 2)...")
x2_sym = sp.integrate(f2 * function_g, (tau, T1, t))

print("Integrating x3 (Pulse 2 Memory)...")
x3_sym = sp.integrate(f2 * function_g, (tau, T1, T1 + T2 / 2))
'''
x_total = sp.Piecewise(
    (x0_sym, (t >= 0) & (t <= T1/2)),
    (x1_sym, (t > T1/2) & (t < T1)),
    (x1_sym + x2_sym, (t >= T1) & (t <= T1 + T2/2)),
    (x1_sym + x3_sym, True)
)
print(x_total)
'''
modules = [{'exp': np.exp}, 'numpy']
x0_func = sp.lambdify(t, x0_sym, modules=modules)
x1_func = sp.lambdify(t, x1_sym, modules=modules)
x2_func = sp.lambdify(t, x2_sym, modules=modules)
x3_func = sp.lambdify(t, x3_sym, modules=modules)

t1 = np.linspace(0, T1/2, 500)
t2 = np.linspace(T1/2, T1, 500)
t3 = np.linspace(T1, T1 + T2/2, 500)
t4 = np.linspace(T1 + T2/2, 5 * T1, 500)

y1 = np.real(x0_func(t1))
y2 = np.real(x1_func(t2))
y3 = np.real(x2_func(t3) + x1_func(t3))
y4 = np.real(x3_func(t4) + x1_func(t4))

print("Plotting results...")

plt.figure(figsize=(10, 6))
plt.plot(t1, y1, 'b', label='Response')
plt.plot(t2, y2, 'b')
plt.plot(t3, y3, 'b')
plt.plot(t4, y4, 'b')

plt.title('Displacement Response x(t)')
plt.xlabel('Time (s)')
plt.ylabel('x(t) [m]')
plt.grid(True)
plt.axvline(T1/2, color='k', linestyle='--', alpha=0.3)
plt.axvline(T1, color='k', linestyle='--', alpha=0.3)
plt.axvline(T1 + T2/2, color='k', linestyle='--', alpha=0.3)
plt.show()