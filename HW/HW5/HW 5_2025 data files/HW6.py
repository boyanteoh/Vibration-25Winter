import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.linalg import solve, inv

# ==========================================
# 1. LOAD DATA
# ==========================================
# NOTE: You need to uncomment the lines below and load your actual .xlsx files
M_df = pd.read_excel('mass_matrix.xlsx', header=None)
C_df = pd.read_excel('damping_matrix.xlsx', header=None)
K_df = pd.read_excel('stiffness_matrix.xlsx', header=None)

# 打印一下形状，看看是不是 (8, 8)
print("Raw M shape:", M_df.shape)

# 自动修正：如果读成了 (8, 9)，说明第一列可能是序号，我们把它删掉
if M_df.shape[1] == 9:
    print("Detected extra column (index), removing it...")
    M = M_df.iloc[:, 1:].values
    C = C_df.iloc[:, 1:].values
    K = K_df.iloc[:, 1:].values
# 如果读成了 (9, 8)，说明第一行可能是表头
elif M_df.shape[0] == 9:
    print("Detected extra row (header), removing it...")
    M = M_df.iloc[1:, :].values
    C = C_df.iloc[1:, :].values
    K = K_df.iloc[1:, :].values
else:
    # 形状正常，直接取值
    M = M_df.values
    C = C_df.values
    K = K_df.values

print("Final M shape:", M.shape) # 这里必须显示 (8, 8)
# 确保 C 和 K 也是 (8, 8)
assert M.shape == (8, 8), "Error: Mass matrix is not 8x8"


# Helper function to get FRF matrix at a specific frequency (in Hz)
def get_frf_matrix(f_hz, M, C, K):
    omega = 2 * np.pi * f_hz
    # Dynamic Stiffness Matrix Z = K - w^2*M + i*w*C
    Z = K - (omega**2 * M) + (1j * omega * C)
    # The FRF is the inverse of Z
    return inv(Z)

# ==========================================
# 2. CALCULATION CASE 1
# ==========================================
# Specs: 8 N p-p force at 1.0 Hz at DOF 7. Find disp at DOF 3.
f1 = 1.0
F_pp_1 = 8.0
dof_force_1 = 6  # DOF 7 is index 6
dof_resp_1 = 2   # DOF 3 is index 2

# Force Amplitude
F_amp_1 = F_pp_1 / 2

# Get FRF matrix at 1.0 Hz
H_1 = get_frf_matrix(f1, M, C, K)

# Displacement Amplitude X = H * F
# Since F has only one non-zero component at dof_force_1:
# X_3 = H_37 * F_7
X_amp_1 = H_1[dof_resp_1, dof_force_1] * F_amp_1

# Convert back to peak-to-peak
disp_pp_1 = 2 * np.abs(X_amp_1)

print(f"--- Case 1 Results ---")
print(f"Excitation: {F_pp_1}N p-p at {f1} Hz (DOF 7)")
print(f"Displacement at DOF 3: {disp_pp_1:.6f} m (peak-to-peak)")
print(f"Displacement at DOF 3: {disp_pp_1*1000:.6f} mm (peak-to-peak)")

# ==========================================
# 3. CALCULATION CASE 2 & PHASE
# ==========================================
# Specs: 16 N p-p force at 2.85 Hz at DOF 3. Find disp at DOF 3.
f2 = 2.85
F_pp_2 = 16.0
dof_force_2 = 2 # DOF 3
dof_resp_2 = 2  # DOF 3

F_amp_2 = F_pp_2 / 2
H_2 = get_frf_matrix(f2, M, C, K)

# Complex displacement
X_complex_2 = H_2[dof_resp_2, dof_force_2] * F_amp_2

# Peak-to-peak
disp_pp_2 = 2 * np.abs(X_complex_2)

# Phase Angle (Lag)
# Angle in radians. If result is negative (e.g. -0.5), it lags by 0.5.
phase_rad = np.angle(X_complex_2)
# Usually 'lag' implies the positive delay.
# If phase is -1.2 rad, lag is 1.2 rad.
phase_lag_deg = np.degrees(-phase_rad) 

print(f"\n--- Case 2 Results ---")
print(f"Excitation: {F_pp_2}N p-p at {f2} Hz (DOF 3)")
print(f"Displacement at DOF 3: {disp_pp_2:.6f} m (peak-to-peak)")
print(f"Displacement at DOF 3: {disp_pp_2*1000:.6f} mm (peak-to-peak)")
print(f"Phase Lag: {-phase_rad:.4f} radians ({phase_lag_deg:.2f} degrees)")


# ==========================================
freqs = np.linspace(0, 10, 400) # 0 to 10 Hz, 400 points
mag_33 = [] # Resp 3 due to Force 3 (Blue)
mag_37 = [] # Resp 3 due to Force 7 (Orange)
mag_77 = [] # Resp 7 due to Force 7 (Green)

for f in freqs:
    H = get_frf_matrix(f, M, C, K)
    
    # We want the magnitude (absolute value) of the FRF elements
    # Note: We plot the FRF magnitude itself (response due to 1N excitation)
    mag_33.append(np.abs(H[2, 2])) 
    mag_37.append(np.abs(H[2, 6]))
    mag_77.append(np.abs(H[6, 6]))

plt.figure(figsize=(10, 6))
plt.semilogy(freqs, mag_33, 'b-', label='Resp at 3 / Force at 3')
plt.semilogy(freqs, mag_37, 'orange', label='Resp at 3 / Force at 7')
plt.semilogy(freqs, mag_77, 'g-', label='Resp at 7 / Force at 7')

plt.title('Frequency Response Functions')
plt.xlabel('Frequency (Hz)')
plt.ylabel('Magnitude (m/N) [Log Scale]')
plt.xlim([0, 10])
plt.ylim([4.0e-6, 0.4]) # As requested
plt.grid(True, which="both", ls="-", alpha=0.5)
plt.legend()
plt.tight_layout()

print("\nPlot generated. Please save or view the figure.")
plt.show()