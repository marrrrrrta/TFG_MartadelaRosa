import numpy as np
import pandas as pd
import os
import glob
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from types import SimpleNamespace
from scipy.signal import find_peaks

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

## FUNCTIONS
import sys
sys.path.append('Functions')
from Functions.diff_eqs_freqfactor import diff_eqs_freqfactor
from Functions.diff_eqs_growth import diff_eqs_growth
from Functions.diff_eqs_notemp import diff_eqs_notemp
from Functions.plotting import plot_results
from Functions.plotting import plot_results_flexible
from Functions.plotting import plot_glowcurve_flexible
from Functions.plotting import plot_column1
from Functions.plotting import plot_column2
from Functions.plotting import plot_column3

#-------------------------------------------------------------------------
#---------------------------- EDITABLE CODE ------------------------------

# Change this to the function you want to use
FunctionUsed = diff_eqs_freqfactor
GraphTitle = r' for a G = 10000 cm$^{-3}$ s$^{-1}$'

# Coment this line if you want to keep the results from last run
'''
results_folder = 'Results'
for file_path in glob.glob(os.path.join(results_folder, '*.png')):
    try:
        os.remove(file_path)
    except Exception as e:
        print(f"Could not delete {file_path}: {e}")
'''
#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

## IMPORTING DATA
# Structural parameters
PathStructural = 'ExperimentalData/ParametrosEstructurales.xlsx'
StructuralData = pd.read_excel(PathStructural, sheet_name='Hoja1', header=0, usecols=None, nrows=None)
# Total density of available positions for traps and Radiative/Non radiative recombination centers (cm-3)
N_I, N_II, N_III, N_IV, N_V, N_s = StructuralData.iloc[[0],[0]].values[0][0],StructuralData.iloc[[1],[0]].values[0][0],StructuralData.iloc[[2],[0]].values[0][0],StructuralData.iloc[[3],[0]].values[0][0],StructuralData.iloc[[4],[0]].values[0][0],StructuralData.iloc[[5],[0]].values[0][0]                    
M_R, M_NR=StructuralData.iloc[[0],[1]].values[0][0],StructuralData.iloc[[0],[2]].values[0][0]
# Electron trapping probability factor for traps and Radiative/Non radiative recombination centers (cm3/s)
A_I, A_II, A_III, A_IV, A_V, A_s = StructuralData.iloc[[0],[3]].values[0][0],StructuralData.iloc[[1],[3]].values[0][0],StructuralData.iloc[[2],[3]].values[0][0],StructuralData.iloc[[3],[3]].values[0][0],StructuralData.iloc[[4],[3]].values[0][0],StructuralData.iloc[[5],[3]].values[0][0]
# Hole trapping probability factor for Radiative/Non Radiative recomb. centers (cm3/s)
A_NR, A_R = StructuralData.iloc[[0],[6]].values[0][0], StructuralData.iloc[[0],[7]].values[0][0]
# Electron-hole trapping probability factor for radiative/non radiative recomb. centers (cm3/s)
A_mn_NR,A_mn_R=StructuralData.iloc[[0],[4]].values[0][0],StructuralData.iloc[[0],[5]].values[0][0]

# Cinetic parameters
PathCinetics = 'ExperimentalData/ParametrosCineticos.xlsx'
CineticsData = pd.read_excel(PathCinetics, sheet_name='Hoja1', header=0, usecols=None, nrows=None)
# Threshold energy for traps and Radiative/Non radiative recombination centers (eV)
E_I,E_II,E_III,E_IV,E_V,E_s=CineticsData.iloc[[0],[0]].values[0][0],CineticsData.iloc[[1],[0]].values[0][0],CineticsData.iloc[[2],[0]].values[0][0],CineticsData.iloc[[3],[0]].values[0][0],CineticsData.iloc[[4],[0]].values[0][0],CineticsData.iloc[[5],[0]].values[0][0]                            
E_R_h, E_NR_h = CineticsData.iloc[[0], [2]].values[0][0], CineticsData.iloc[[0], [4]].values[0][0]           
# Frequency factor for trap i (s-1)
S_I,S_II,S_III,S_IV,S_V,S_s=CineticsData.iloc[[0],[1]].values[0][0],CineticsData.iloc[[1],[1]].values[0][0],CineticsData.iloc[[2],[1]].values[0][0],CineticsData.iloc[[3],[1]].values[0][0],CineticsData.iloc[[4],[1]].values[0][0],CineticsData.iloc[[5],[1]].values[0][0]
S_R_h, S_NR_h = CineticsData.iloc[[0], [3]].values[0][0], CineticsData.iloc[[0], [5]].values[0][0]

#-------------------------------------------------------------------------
#-------------------------------------------------------------------------

## 2. DIFFERENTIAL EQUATIONS SYSTEM
# For the global variables, we need to define a SimpleNamespace
value = SimpleNamespace(
    growth = 1E-6,
    kB=0,
    G=0, hr=0, 
    N_I=N_I, N_II=N_II, N_III=N_III, N_IV=N_IV, N_V=N_V, N_s=N_s,
    M_R=M_R, M_NR=M_NR,
    A_I=A_I, A_II=A_II, A_III=A_III, A_IV=A_IV, A_V=A_V, A_s=A_s, A_NR=A_NR, A_R=A_R,
    E_I=E_I, E_II=E_II, E_III=E_III, E_IV=E_IV, E_V=E_V, E_s=E_s, E_R_h=E_R_h, E_NR_h=E_NR_h,
    S_I=S_I, S_II=S_II, S_III=S_III, S_IV=S_IV, S_V=S_V, S_s=S_s, S_R_h=S_R_h, S_NR_h=S_NR_h,
    A_mn_NR=A_mn_NR, A_mn_R=A_mn_R
)
# Define the columns for the dataframe
column_names = [
    'n_I', 'n_II', 'n_III', 'n_IV', 'n_V',
    'n_s', 'm_R', 'm_NR', 'n_c', 'n_v',
    'dm_R', 'dm_NR'
]

#-------------------------------------------------------------------------

## 2.1 IRRADIATION
# Parameters for IRRADIATION
value.kB = 8.617e-5        # Boltzmann constant (eV/K)
value.T_C = 25             # Temperature (ºC)
value.hr = 0               # Heating rate (ºC/s)
value.G = 1000             # Electron-hole pair generation (cm-3 s-1)

# Time vector (s)
npoints = 3600
t = np.linspace(0, npoints-1, npoints)

# Initial conditions vector
n_I_01, n_II_01, n_III_01, n_IV_01, n_V_01, n_s_01, m_NR_01, m_R_01, n_c_01, n_v_01 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
y01 = [n_I_01, n_II_01, n_III_01, n_IV_01, n_V_01, n_s_01, m_NR_01, m_R_01, n_c_01, n_v_01]

n_I_02, n_II_02, n_III_02, n_IV_02, n_V_02, n_s_02, m_NR_02, m_R_02, n_c_02, n_v_02 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
y02 = [n_I_02, n_II_02, n_III_02, n_IV_02, n_V_02, n_s_02, m_NR_02, m_R_02, n_c_02, n_v_02]

# Solving the differential equations system for TEMPERATURE DEPENDENCY
value.E_I, value.E_II, value.E_III, value.E_IV, value.E_V, value.E_s = E_I, E_II, E_III, E_IV, E_V, E_s
irradiation1 = odeint(FunctionUsed, y01, t, args=(value,))
n_I1, n_II1, n_III1 ,n_IV1 ,n_V1 ,n_s1 ,m_R1 ,m_NR1 ,n_c1 , n_v1 = irradiation1.T

value.E_I, value.E_II, value.E_III, value.E_IV, value.E_V, value.E_s = 1.19, 1.38, 1.68, 1.78, 2.12, 3.00  # Original values
irradiation2 = odeint(diff_eqs_notemp, y02, t, args=(value,))
n_I2, n_II2, n_III2 ,n_IV2 ,n_V2 ,n_s2 ,m_R2 ,m_NR2 ,n_c2 , n_v2 = irradiation2.T

'''# Plotting the results
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
plot_results_flexible(irradiation1, t, value, axes[0], 'upper left')
labels = ['(a)', '(b)', '(c)']
for ax, label in zip(axes[0], labels):
    ax.text(-0.1, 1.07, label, transform=ax.transAxes, fontsize=12, va='top')
plot_results_flexible(irradiation2, t, value, axes[1], 'upper left')
labels = ['(d)', '(e)', '(f)']
for ax, label in zip(axes[1], labels):
    ax.text(-0.1, 1.07, label, transform=ax.transAxes, fontsize=12, va='top')
fig.text(0.5, 1.01, 'Irradiation dependent (top) and independent (bottom) on temperature', ha='center', fontsize=13, weight='bold')
plt.tight_layout()
#plt.savefig('LaTeX/Images/Irradiation_Comparison.png', dpi=300, bbox_inches='tight')
#plt.show()'''

#plot_column1(irradiation1, irradiation2, 'LaTeX/Images/Irradiation ', t, 'upper left', 'Time (s)', 'Irradiation phase')
#plot_column2(irradiation1, irradiation2, 'LaTeX/Images/Irradiation ', t, value, 'upper left', 'Time (s)', 'Irradiation phase')
#plot_column3(irradiation1, irradiation2, 'LaTeX/Images/Irradiation ', t, value, 'Time (s)', 'Irradiation phase')

#-------------------------------------------------------------------------

## 2.2 RELAXATION
# Parameters for RELAXATION
value.T_C = 25             # Temperature (ºC)
value.hr = 0               # Heating rate (ºC/s)
value.G = 0                # Electron-hole pair generation (cm-3 s-1)

# Time vector (s)
npoints = 3600 * 24 * 7
t = np.linspace(0, npoints-1, npoints)

# Initial conditions vector
n_I_01, n_II_01, n_III_01, n_IV_01, n_V_01, n_s_01, m_NR_01, m_R_01, n_c_01, n_v_01 = n_I1[-1], n_II1[-1], n_III1[-1], n_IV1[-1], n_V1[-1], n_s1[-1], m_R1[-1], m_NR1[-1], n_c1[-1], n_v1[-1]
y01 = [n_I_01, n_II_01, n_III_01, n_IV_01, n_V_01, n_s_01, m_NR_01, m_R_01, n_c_01, n_v_01]

n_I_02, n_II_02, n_III_02, n_IV_02, n_V_02, n_s_02, m_NR_02, m_R_02, n_c_02, n_v_02 = n_I2[-1], n_II2[-1], n_III2[-1], n_IV2[-1], n_V2[-1], n_s2[-1], m_R2[-1], m_NR2[-1], n_c2[-1], n_v2[-1]
y02 = [n_I_02, n_II_02, n_III_02, n_IV_02, n_V_02, n_s_02, m_NR_02, m_R_02, n_c_02, n_v_02]

# Solving the differential equations system
value.E_I, value.E_II, value.E_III, value.E_IV, value.E_V, value.E_s = E_I, E_II, E_III, E_IV, E_V, E_s
relaxation1 = odeint(FunctionUsed, y01, t, args=(value,))
n_I1, n_II1, n_III1 ,n_IV1 ,n_V1 ,n_s1 ,m_R1 ,m_NR1 ,n_c1 , n_v1 = relaxation1.T

value.E_I, value.E_II, value.E_III, value.E_IV, value.E_V, value.E_s = 1.19, 1.38, 1.68, 1.78, 2.12, 3.00  # Original values
relaxation2 = odeint(diff_eqs_notemp, y02, t, args=(value,))
n_I2, n_II2, n_III2 ,n_IV2 ,n_V2 ,n_s2 ,m_R2 ,m_NR2 ,n_c2 , n_v2 = relaxation2.T

'''# Plotting the results
fig, axes = plt.subplots(2, 3, figsize=(15, 9))
plot_results_flexible(relaxation1, t, value, axes[0], 'upper right')
labels = ['(a)', '(b)', '(c)']
for ax, label in zip(axes[0], labels):
    ax.text(-0.1, 1.07, label, transform=ax.transAxes, fontsize=12, va='top')
plot_results_flexible(relaxation2, t, value, axes[1], 'upper right')
labels = ['(d)', '(e)', '(f)']
for ax, label in zip(axes[1], labels):
    ax.text(-0.1, 1.07, label, transform=ax.transAxes, fontsize=12, va='top')
fig.text(0.5, 1.01, 'Relaxation dependent (top) and independent (bottom) on temperature', ha='center', fontsize=13, weight='bold')
plt.tight_layout()
#plt.savefig('LaTeX/Images/Relaxation_Comparison.png', dpi=300, bbox_inches='tight')
#plt.show()'''

#plot_column1(relaxation1, relaxation2, 'LaTeX/Images/Relaxation ', t, 'center right', 'Time (s)', 'Relaxation phase')
#plot_column2(relaxation1, relaxation2, 'LaTeX/Images/Relaxation ', t, value, 'upper right', 'Time (s)', 'Relaxation phase')
#plot_column3(relaxation1, relaxation2, 'LaTeX/Images/Relaxation ', t, value, 'Time (s)', 'Relaxation phase')

#-------------------------------------------------------------------------

## 2.3 HEATING
# Parameters for HEATING
value.T_C = 0             # Temperature (ºC)
value.hr = 1.0             # Heating rate (ºC/s)
value.G = 0                # Electron-hole pair generation (cm-3 s-1)

# Time vector (s)
npoints = 400
t = np.linspace(0, npoints-1, npoints)

# Initial conditions vector
n_I_01, n_II_01, n_III_01, n_IV_01, n_V_01, n_s_01, m_NR_01, m_R_01, n_c_01, n_v_01 = n_I1[-1], n_II1[-1], n_III1[-1], n_IV1[-1], n_V1[-1], n_s1[-1], m_R1[-1], m_NR1[-1], n_c1[-1], n_v1[-1]
y01 = [n_I_01, n_II_01, n_III_01, n_IV_01, n_V_01, n_s_01, m_NR_01, m_R_01, n_c_01, n_v_01]

n_I_02, n_II_02, n_III_02, n_IV_02, n_V_02, n_s_02, m_NR_02, m_R_02, n_c_02, n_v_02 = n_I2[-1], n_II2[-1], n_III2[-1], n_IV2[-1], n_V2[-1], n_s2[-1], m_R2[-1], m_NR2[-1], n_c2[-1], n_v2[-1]
y02 = [n_I_02, n_II_02, n_III_02, n_IV_02, n_V_02, n_s_02, m_NR_02, m_R_02, n_c_02, n_v_02]

# Solving the differential equations system
value.E_I, value.E_II, value.E_III, value.E_IV, value.E_V, value.E_s = E_I, E_II, E_III, E_IV, E_V, E_s
heating1 = odeint(FunctionUsed, y01, t, args=(value,))
n_I1, n_II1, n_III1 ,n_IV1 ,n_V1 ,n_s1 ,m_R1 ,m_NR1 ,n_c1 , n_v1 = heating1.T
'''dm_R = m_R1 * value.A_mn_R * n_c1
dm_NR = m_NR1 * value.A_mn_NR * n_c1
heating1 = np.column_stack((heating1, dm_R, dm_NR))'''

value.E_I, value.E_II, value.E_III, value.E_IV, value.E_V, value.E_s = 1.19, 1.38, 1.68, 1.78, 2.12, 3.00  # Original values
heating2 = odeint(diff_eqs_notemp, y02, t, args=(value,))
n_I2, n_II2, n_III2 ,n_IV2 ,n_V2 ,n_s2 ,m_R2 ,m_NR2 ,n_c2 , n_v2 = heating2.T
'''dm_R = m_R2 * value.A_mn_R * n_c2
dm_NR = m_NR2 * value.A_mn_NR * n_c2
heating2 = np.column_stack((heating2, dm_R, dm_NR))'''

#plot_column1(heating1, heating2, 'LaTeX/Images/Heating ', t, 'upper left', 'Temperature (ºC)', 'Heating phase')
#plot_column2(heating1, heating2, 'LaTeX/Images/Heating ', t, value, 'upper right', 'Temperature (ºC)', 'Heating phase')
plot_column3(heating1, heating2, 'LaTeX/Images/Heating ', t, value, 'Temperature (ºC)', 'Heating phase')

# Saving and exporting to excel
results1 = pd.DataFrame(heating1, columns=column_names)
results2 = pd.DataFrame(heating2, columns=column_names)
results1.to_excel('Results/Heating_Temperature_Dependent.xlsx', index=False)
results2.to_excel('Results/Heating_Temperature_Independent.xlsx', index=False)


#-------------------------------------------------------------------------
#-------------------------------------------------------------------------
## 3. GLOW CURVE ANALYSIS

ni_curves1 = np.array([n_I1, n_II1, n_III1, n_IV1, n_V1, n_s1]).T 
ni_curves2 = np.array([n_I2, n_II2, n_III2, n_IV2, n_V2, n_s2]).T

#-------------------------------------------------------------------------
## 3.1  ACTIVATION TEMPERATURES

'''activation_temperatures = []
trap_labels = ['n$_I(t)$', 'n$_{II}(t)$', 'n$_{III}(t)$', 'n$_{IV}(t)$', 'n$_V(t)$']
colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']

plt.figure(figsize=(15, 8))
plt.subplot(1, 2, 1)
plt.gca().text(0.5, -0.13, '(a)', transform=plt.gca().transAxes, fontsize=12, va='top', ha='center')
plt.plot(t, ni_curves1[:, 5], label='n$_s(t)$', color=colors[5])
for i in reversed(range(5)):  # exclude n_s
    n = ni_curves1[:, i]
    dn_dT = -np.gradient(n, t)  # negative gradient

    # 1. Glow peak (max detrapping rate)
    T_peak = t[np.argmax(dn_dT)]
    plt.plot(T_peak, n[np.argmax(dn_dT)], '*', color=colors[i], markersize=10, label=f'Peak {trap_labels[i]}')
    plt.plot([T_peak, T_peak], [0, n[np.argmax(dn_dT)]], linestyle='dotted', color=colors[i])
    plt.text(T_peak, -0.05e6, f'{T_peak:.0f}°C', ha='center', va='top', fontsize=9, color=colors[i])
    
    # 2. Activation temperature (start of decrease)
    start_idx = np.argmax(dn_dT > 1e-2)  # threshold to detect slope change
    T_start = t[start_idx]
    plt.plot(T_start, n[start_idx], 'o', color=colors[i], markersize=6, label=f'Activation {trap_labels[i]}')
    
    # 3. Plot n_i(T)
    plt.plot(t, n, label=trap_labels[i], color=colors[i])

plt.xlabel('Temperature (°C)', fontsize=14)
plt.ylabel('Trap concentration (cm$^{-3}$)', fontsize=14)
plt.title('Temperature dependent model', fontsize=16)
plt.legend(reverse=True, loc='upper left')

plt.subplot(1, 2, 2)
plt.gca().text(0.5, -0.13, '(b)', transform=plt.gca().transAxes, fontsize=12, va='top', ha='center')
plt.plot(t, ni_curves2[:, 5], label='n$_s(t)$', color=colors[5])
for i in reversed(range(5)):  # exclude n_s
    n = ni_curves2[:, i]
    dn_dT = -np.gradient(n, t)  # negative gradient
    
    # 1. Glow peak (max detrapping rate)
    T_peak = t[np.argmax(dn_dT)]
    plt.plot(T_peak, n[np.argmax(dn_dT)], '*', color=colors[i], markersize=10, label=f'Peak {trap_labels[i]}')
    plt.plot([T_peak, T_peak], [0, n[np.argmax(dn_dT)]], linestyle='dotted', color=colors[i])
    plt.text(T_peak, -0.05e6, f'{T_peak:.0f}°C', ha='center', va='top', fontsize=9, color=colors[i])
    
    # 2. Activation temperature (start of decrease)
    start_idx = np.argmax(dn_dT > 1e-2)  # threshold to detect slope change
    T_start = t[start_idx]
    plt.plot(T_start, n[start_idx], 'o', color=colors[i], markersize=6, label=f'Activation {trap_labels[i]}')
    
    # 3. Plot n_i(T)
    plt.plot(t, n, label=trap_labels[i], color=colors[i])
    
plt.xlabel('Temperature (°C)', fontsize=14)
plt.ylabel('Trap concentration (cm$^{-3}$)', fontsize=14)
plt.title('Temperature independent model', fontsize=16)
plt.legend(reverse=True, loc='upper left')
plt.suptitle(r'n$_{i}$ evolution for the Heating phase', fontsize=18)
plt.tight_layout()
plt.savefig('LaTeX/Images/GC_ActivationAndPeakTemperatures.png', dpi=600, bbox_inches='tight', transparent=True)
plt.show()

#-------------------------------------------------------------------------
## 3.2  GLOW CURVE

dm_R1 = m_R1 * value.A_mn_R * n_c1
dm_NR1 = m_NR1 * value.A_mn_NR * n_c1

dm_R2 = m_R2 * value.A_mn_R * n_c2
dm_NR2 = m_NR2 * value.A_mn_NR * n_c2

# Find all maximums
peaks1, _ = find_peaks(dm_R1, prominence=50)
peaks2, _ = find_peaks(dm_R2, prominence=50)
T_peaks1 = t[peaks1]
T_peaks2 = t[peaks2]
I_peaks1 = dm_R1[peaks1]
I_peaks2 = dm_R2[peaks2]

colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:purple', 'tab:red', 'tab:brown']
# Plotting the glow curves
plt.figure(figsize=(15, 8))
plt.subplot(1, 2, 1)
plt.gca().text(0.5, -0.13, '(a)', transform=plt.gca().transAxes, fontsize=12, va='top', ha='center')
plt.plot(t, dm_R1, label=r'dm$_R$(t)')
plt.plot(t, dm_NR1, label=r'dm$_{NR}$(t)')
for t_peak, i_peak, i in zip(T_peaks1, I_peaks1, range(len(T_peaks1))):
    plt.plot(t_peak, i_peak, '*', markersize=10, color=colors[i])
    plt.plot([t_peak, t_peak], [0, i_peak], linestyle='dotted', color=colors[i])
    plt.plot(t_peak, 0, 'o', color=colors[i])
    plt.text(t_peak, -0.02 * max(dm_R1), f'{t_peak:.0f}°C', ha='center', va='top', fontsize=9, color=colors[i])
    
manual_T = 256
manual_I = dm_R1[manual_T]
plt.plot(manual_T, manual_I, '*', markersize=10, color=colors[4])
plt.plot([manual_T, manual_T], [0, manual_I], linestyle='dotted', color=colors[4])
plt.plot(manual_T, 0, 'o', color=colors[4])
plt.text(manual_T, -0.02 * max(dm_R1), '256°C', ha='center', va='top', fontsize=9, color=colors[4])

plt.legend()
plt.xlabel('Temperature (°C)', fontsize=14)
plt.ylabel('Intensity [u.a.]', fontsize=14)
plt.title('Temperature dependent model', fontsize=16)

plt.subplot(1, 2, 2)
plt.gca().text(0.5, -0.13, '(b)', transform=plt.gca().transAxes, fontsize=12, va='top', ha='center')
plt.plot(t, dm_R2, label=r'dm$_R$(T)')
plt.plot(t, dm_NR2, label=r'dm$_{NR}$(T)')
for t_peak, i_peak, i in zip(T_peaks2, I_peaks2, range(len(T_peaks2))):
    plt.plot(t_peak, i_peak, '*', markersize=10, color=colors[i])
    plt.plot([t_peak, t_peak], [0, i_peak], linestyle='dotted', color=colors[i])
    plt.plot(t_peak, 0, 'o', color=colors[i])
    plt.text(t_peak, -0.02 * max(dm_R2), f'{t_peak:.0f}°C', ha='center', va='top', fontsize=9, color=colors[i])

manual_T = 269
manual_I = dm_R2[manual_T]
plt.plot(manual_T, manual_I, '*', markersize=10, color=colors[4])
plt.plot([manual_T, manual_T], [0, manual_I], linestyle='dotted', color=colors[4])
plt.plot(manual_T, 0, 'o', color=colors[4])
plt.text(manual_T, -0.02 * max(dm_R2), '269°C', ha='center', va='top', fontsize=9, color=colors[4])
plt.xlabel('Temperature (°C)', fontsize=14)
plt.ylabel('Intensity [u.a.]', fontsize=14)
plt.title('Temperature independent model', fontsize=16)

plt.suptitle('Recombination for the Heating phase', fontsize=18)
plt.legend()
plt.tight_layout()
plt.savefig('LaTeX/Images/GC_GlowCurve.png', dpi=600, bbox_inches='tight')
plt.show()'''


