import numpy as np
import pandas as pd
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter1d

# ------------------------------------------------------------------------------------
# 1. FERMI-DIRAC DISTRIBUTION WHEN IRRADIATED
## FOR THE CASE OF A NON ZERO TEMPERATURE

# Energy and temperature range
E = np.linspace(-1, 1, 500)
Ts = [0.1, 300, 600]  # K

# Fermi-Dirac distribution at T > 0 K
def fermi_T(E, Ef, T):
    kB = 8.617e-5  # eV/K
    return 1 / (1 + np.exp((E - Ef) / (kB * T)))

# State during irradiation: smooth transition between EFe and EFh
def irradiated_occupancy_T_temp(E, Efe, Efh, T):
    """
    Constructs a piecewise occupancy function with smoothing based on temperature.
    Sharp at low T, smooth at high T. Matches visual expectation.
    """
    f = np.zeros_like(E)
    mask = (E >= Efh) & (E <= Efe)
    f[E < Efh] = 1
    f[mask] = 1 - (E[mask] - Efh) / (Efe - Efh)  # Linear drop from 1 to 0

    # Smooth with Gaussian kernel whose width depends on T
    sigma = T / 10  # or adjust for visual consistency
    return gaussian_filter1d(f, sigma=sigma)


# Level definitions
Ef = 0.0                 # original Fermi level
Efe_initial = 0.5
Efh_initial = -0.5
Efe_mid = 0.2
Efh_mid = -0.2

# Temperatures
T = 300  # K
f_before = np.zeros((len(Ts), len(E)))
f_after = np.zeros((len(Ts), len(E)))
f_during = np.zeros((len(Ts), len(E)))
f_final = np.zeros((len(Ts), len(E)))

# Functions for plotting
for t in range(len(Ts)):
    f_before[t] = fermi_T(E, Ef, Ts[t])
    f_after[t] = irradiated_occupancy_T_temp(E, Efe_initial, Efh_initial, Ts[t])
    f_during[t] = irradiated_occupancy_T_temp(E, Efe_mid, Efh_mid, Ts[t])
    f_final[t] = fermi_T(E, Ef, Ts[t])

    # Plot the curves for each temperature in the for loop

# Plot
fig, axs = plt.subplots(1, 4, figsize=(20, 5), sharey=True)
titles = ['Before irradiation', 'After irradiation', 'During stimulation', 'Final equilibrium']
tags = ['(a)', '(b)', '(c)', '(d)']
occupancies = [f_before, f_after, f_during, f_final]
labels = ['EF', (Efe_initial, Efh_initial), (Efe_mid, Efh_mid), 'EF']

for i in range(4):
    for t in range(len(Ts)):
        axs[i].plot(E, occupancies[i][t], label=f'T = {Ts[t]} K')
    axs[i].set_title(titles[i], fontsize=14)
    axs[i].set_xlabel('Energy', fontsize=12)
    axs[i].set_ylim(-0.1, 1.1)
    axs[i].set_xticks([])
    axs[i].text(0.5, -0.13, tags[i], transform=axs[i].transAxes, fontsize=12, va='top', ha='center')

    # Vertical lines for Fermi levels
    if i == 1 or i == 2:
        Efe, Efh = labels[i]
        axs[i].axvline(Efe, color='green', linestyle='--', label='EFe')
        axs[i].axvline(Efh, color='purple', linestyle='--', label='EFh')
        axs[i].axvline(Ef, color='black', linestyle=':', label='EF')
    else:
        axs[i].axvline(Ef, color='black', linestyle=':', label='EF')

    if i == 0:
        axs[i].set_ylabel('Occupancy f(E)', fontsize=12)
    axs[i].legend(loc='lower left')
plt.tight_layout()

fig.savefig('LaTeX/Images/FD_irradiation.png', dpi=600, bbox_inches='tight')
plt.show()



# ------------------------------------------------------------------------------------
# 2. PHOTON ENERGY CALCULATION
nu_m = 90000                # from B.H. Bransden book (m-1)
nu_s = 2.9979e8 * nu_m      # (s-1)
h = 4.135667696e-15         # Planck's constant (eV·s)
S_hr = 10                   # Huang-Rhys factor

E_ph = S_hr * h * nu_s      # Photon energy (eV)
#print(f'Photon energy: {E_ph:.2f} eV')

# ------------------------------------------------------------------------------------
## 3. EXPERIMENTAL GLOW CURVE
'''
# Importing data
PathExcel = 'ExperimentalData/DatEx_TLD_100_Beta_1.xlsx'
ExperimentalData = pd.read_excel(PathExcel, sheet_name='GlowCurve', header=0, usecols='A:C')

T = ExperimentalData['T'].values
I = ExperimentalData['Glow'].values

# Find peaks
from scipy.signal import find_peaks
peaks, _ = find_peaks(I, prominence=500)
t_peaks = T[peaks]
i_peaks = I[peaks]

## Plotting
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple', 'tab:red', 'tab:brown']
label = ['Peak 1: 73ºC', 'Peak 2: 121ºC', 'Peak 3: 156ºC', 'Peak 5: 202ºC ', 'Peak 4: 182ºC']

plt.figure(figsize=(10, 6))
plt.plot(T, I, 'o-', label='Experimental Data', markersize=3, color='forestgreen')
for t_peak, i_peak, i in zip(t_peaks, i_peaks, range(len(t_peaks))):
    plt.plot(t_peak, i_peak, '*', markersize=10, color=colors[i], label=label[i])
    plt.plot([t_peak, t_peak], [0, i_peak], linestyle='dotted', color=colors[i])
    plt.plot(t_peak, 0, 'o', color=colors[i])
    plt.text(t_peak, -0.02 * max(I), f'{t_peak:.0f}°C', ha='center', va='top', fontsize=9, color=colors[i])

manual_T = 182.547
manual_I = 6812.25
plt.plot(manual_T, manual_I, '*', markersize=10, color=colors[4], label=label[4])
plt.plot([manual_T, manual_T], [0, manual_I], linestyle='dotted', color=colors[4])
plt.plot(manual_T, 0, 'o', color=colors[4])
plt.text(manual_T, -0.02 * max(I), '182ºC', ha='center', va='top', fontsize=9, color=colors[4])


plt.title('Experimental Glow Curve')
plt.xlabel('Temperature (ºC)')
plt.ylabel('Intensity (a.u.)')
plt.legend()
plt.savefig('LaTeX/Images/Experimental_Glow_Curve.png', dpi=300)
plt.show()'''
