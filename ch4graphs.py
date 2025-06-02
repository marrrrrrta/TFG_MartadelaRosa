import numpy as np
import pandas as pd
import os
import glob
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from types import SimpleNamespace

import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

## FUNCTIONS
import sys
sys.path.append('Functions')
from Functions.diff_eqs_freqfactor import diff_eqs_freqfactor
from Functions.diff_eqs_growth import diff_eqs_growth
from Functions.diff_eqs_notemp import diff_eqs_notemp
from Functions.plotting import plot_results

#-------------------------------------------------------------------------
#---------------------------- EDITABLE CODE ------------------------------

# Change this to the function you want to use
FunctionUsed = diff_eqs_freqfactor
GraphTitle = r' for a G = 10000 cm$^{-3}$ s$^{-1}$'
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


# Time vector (s)
npoints = 3600
t = np.linspace(0, npoints-1, npoints)

# Initial conditions vector
n_I_01, n_II_01, n_III_01, n_IV_01, n_V_01, n_s_01, m_NR_01, m_R_01, n_c_01, n_v_01 = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0
y01 = [n_I_01, n_II_01, n_III_01, n_IV_01, n_V_01, n_s_01, m_NR_01, m_R_01, n_c_01, n_v_01]

# Solving the differential equations system for TEMPERATURE DEPENDENCY
value.G = 10000            # Electron-hole pair generation (cm-3 s-1)
irradiation1 = odeint(FunctionUsed, y01, t, args=(value,))
n_I1, n_II1, n_III1 ,n_IV1 ,n_V1 ,n_s1 ,m_R1 ,m_NR1 ,n_c1 , n_v1 = irradiation1.T

value.G = 1000
irradiation2 = odeint(FunctionUsed, y01, t, args=(value,))
n_I2, n_II2, n_III2 ,n_IV2 ,n_V2 ,n_s2 ,m_R2 ,m_NR2 ,n_c2 , n_v2 = irradiation2.T

# Plotting the results
# Frequency factor INDEPENDENT of temperature
plt.figure(figsize=(10, 4))

plt.subplot(1, 3, 1)
plt.plot(t, n_I1, label='n$_{I}$(t)', color='blue')
plt.plot(t, n_I2, label='n$_{I}$(t) 2', color='orange')
plt.plot(t, n_II1, label='n$_{II}$(t)', color='blue')
plt.plot(t, n_II2, label='n$_{II}$(t) 2', color='orange')
plt.plot(t, n_III1, label='n$_{III}$(t)', color='blue')
plt.plot(t, n_III2, label='n$_{III}$(t) 2', color='orange')
plt.plot(t, n_IV1, label='n$_{IV}$(t)', color='blue')
plt.plot(t, n_IV2, label='n$_{IV}$(t) 2', color='orange')
plt.plot(t, n_V1, label='n$_{V}$(t)', color='blue')
plt.plot(t, n_V2, label='n$_{V}$(t) 2', color='orange')
plt.plot(t, n_s1, label='n$_{s}$(t)', color='blue')
plt.plot(t, n_s2, label='n$_{s}$(t) 2', color='orange')

plt.xlabel('Time (s)')
plt.ylabel('Trap concentration (cm$^{-3}$)')
plt.title('n$_{i}$ evolution')
plt.legend(loc = 'upper left')

plt.subplot(1, 3, 2)
dm_R1 = m_R1 * A_mn_R * n_c1
dm_R2 = m_R2 * A_mn_R * n_c2
dm_NR1 = m_NR1 * A_mn_NR * n_c1
dm_NR2 = m_NR2 * A_mn_NR * n_c2
plt.plot(t, dm_R1, label='dm$_{R}$(t)', color='blue')
plt.plot(t, dm_R2, label='dm$_{R}$(t) 2', color='orange')
plt.plot(t, dm_NR1, label='dm$_{NR}$(t)', color='blue')
plt.plot(t, dm_NR2, label='dm$_{NR}$(t) 2', color='orange')

plt.xlabel('Time (s)')
plt.ylabel(r'$\frac{dm_{i=R,NR}}{dt}$ [u.a.]')
plt.title('Recombination')
plt.legend(loc = 'upper left')

plt.subplot(1, 3, 3)
plt.plot(t, (n_c1 + n_I1 + n_II1 + n_III1 + n_IV1 + n_V1 + n_s1)/(m_R1 + m_NR1), color='blue')
plt.plot(t, (n_c2 + n_I2 + n_II2 + n_III2 + n_IV2 + n_V2 + n_s2)/(m_R2 + m_NR2), color='orange')
plt.xlabel('Time (s)')
plt.title('Charge neutrality')

plt.suptitle('Frequency factor INDEPENDENT of temperature')

plt.tight_layout()
plt.show()




'''
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
n_I_0,n_II_0,n_III_0,n_IV_0,n_V_0,n_s_0,m_NR_0,m_R_0,n_c_0,n_v_0=n_I[-1],n_II[-1],n_III[-1],n_IV[-1],n_V[-1],n_s[-1],m_R[-1],m_NR[-1],n_c[-1], n_v[-1]
y0 = [n_I_0,n_II_0,n_III_0,n_IV_0,n_V_0,n_s_0,m_NR_0,m_R_0,n_c_0,n_v_0]

# Solving the differential equations system
relaxation = odeint(FunctionUsed, y0, t, args=(value,))
n_I, n_II, n_III ,n_IV ,n_V ,n_s ,m_R ,m_NR ,n_c , n_v = relaxation.T

# Plotting the results
plot_results(relaxation, 'Results/', 'RELAXATION' + GraphTitle, t, value)

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
n_I_0,n_II_0,n_III_0,n_IV_0,n_V_0,n_s_0,m_NR_0,m_R_0,n_c_0,n_v_0=n_I[1],n_II[1],n_III[1],n_IV[1],n_V[1],n_s[1],m_R[1],m_NR[1],n_c[1], n_v[1]
y0 = [n_I_0,n_II_0,n_III_0,n_IV_0,n_V_0,n_s_0,m_NR_0,m_R_0,n_c_0,n_v_0]

# Solving the differential equations system
heating = odeint(FunctionUsed, y0, t, args=(value,))
n_I, n_II, n_III ,n_IV ,n_V ,n_s ,m_R ,m_NR ,n_c , n_v = heating.T

# Plotting the results
temp_plot = value.T_C + value.hr * t
plot_results(heating, 'Results/', 'HEATING' + GraphTitle, temp_plot, value)
'''