# Plots and saves the results to the 'Results' folder
from matplotlib import gridspec
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import find_peaks

def plot_results(odeint_solution, save_path, title, x_axis, value):
    """
    Plots the three main graphs of a certain stage. The results are saved in the 'Results' folder with the given title. 
    
    The plots are:
    - n_i evolution
    - Recombination
    - Charge neutrality ratio
    """    
    # Unpack the variables
    n_I, n_II, n_III ,n_IV ,n_V ,n_s ,m_R ,m_NR ,n_c , n_v = odeint_solution.T
    
    # Intensity of the glow curve (eq. 3.7)
    dm_R = m_R * value.A_mn_R * n_c
    dm_NR = m_NR * value.A_mn_NR * n_c
    odeint_solution = np.column_stack((odeint_solution, dm_R, dm_NR))
    
    # Plotting
    plt.figure(figsize=(15, 8))
    
    plt.subplot(1, 3, 1)
    plt.plot(x_axis, n_I, label='n$_{I}$(t)')
    plt.plot(x_axis, n_II, label='n$_{II}$(t)')
    plt.plot(x_axis, n_III, label='n$_{III}$(t)')
    plt.plot(x_axis, n_IV, label='n$_{IV}$(t)')
    plt.plot(x_axis, n_V, label='n$_{V}$(t)')
    plt.plot(x_axis, n_s, label='n$_{s}$(t)')
    plt.xlabel('Time (s)', fontsize=14)
    plt.ylabel('Trap concentration (cm$^{-3}$)', fontsize=14)
    plt.title('n$_{i}$ evolution', fontsize=16)
    plt.legend(loc = 'upper left')

    plt.subplot(1, 3, 2)
    plt.plot(x_axis, dm_R, label='dm$_{R}$(t)')
    plt.plot(x_axis, dm_NR, label='dm$_{NR}$(t)')
    plt.xlabel('Time (s)', fontsize=14)
    plt.ylabel(r'$\frac{dm_{i=R,NR}}{dt}$ [u.a.]', fontsize=14)
    plt.title('Recombination', fontsize=16)
    plt.legend(loc = 'upper left')

    plt.subplot(1, 3, 3)
    plt.plot(x_axis, (n_c + n_I + n_II + n_III + n_IV + n_V + n_s)/(m_R + m_NR + n_v))
    plt.xlabel('Time (s)', fontsize=14)
    plt.title('Charge neutrality', fontsize=16)
    plt.suptitle(title)
    plt.tight_layout()
    
    # Saving the results
    plt.savefig(save_path + title + '.png', dpi=600, bbox_inches='tight')
    plt.show()
    
def plot_glowcurve(odeint_solution, save_path, title, x_axis, value):
    """
    Plots the three main graphs for the heating stage. The results are saved in the 'Results' folder with the given title. 
    
    The plots are:
    - n_i evolution
    - Recombination
    - Charge neutrality ratio
    """
    # Unpack the variables
    n_I, n_II, n_III ,n_IV ,n_V ,n_s ,m_R ,m_NR ,n_c , n_v = odeint_solution.T
    
    # Intensity of the glow curve (eq. 3.7)
    dm_R = m_R * value.A_mn_R * n_c
    dm_NR = m_NR * value.A_mn_NR * n_c
    odeint_solution = np.column_stack((odeint_solution, dm_R, dm_NR))
    
    # Plotting
    plt.figure(figsize=(15, 8))
    
    plt.subplot(1, 3, 1)
    plt.plot(x_axis, n_I, label='n$_{I}$(t)')
    plt.plot(x_axis, n_II, label='n$_{II}$(t)')
    plt.plot(x_axis, n_III, label='n$_{III}$(t)')
    plt.plot(x_axis, n_IV, label='n$_{IV}$(t)')
    plt.plot(x_axis, n_V, label='n$_{V}$(t)')
    plt.plot(x_axis, n_s, label='n$_{s}$(t)')
    plt.xlabel('Temperature (ºC)', fontsize=14)
    plt.ylabel('Trap concentration (cm$^{-3}$)', fontsize=14)
    plt.title('n$_{i}$ evolution', fontsize=16)
    plt.legend(loc = 'upper left')

    plt.subplot(1, 3, 2)
    plt.plot(x_axis, dm_R, label='dm$_{R}$(t)')
    plt.plot(x_axis, dm_NR, label='dm$_{NR}$(t)')
    plt.xlabel('Temperature (ºC)', fontsize=14)
    plt.ylabel('I [u.a.]', fontsize=14)
    plt.title('Recombination', fontsize=16)
    plt.legend(loc = 'upper left')

    plt.subplot(1, 3, 3)
    plt.plot(x_axis, (n_c + n_I + n_II + n_III + n_IV + n_V + n_s)/(m_R + m_NR + n_v))
    plt.xlabel('Temperature (ºC)', fontsize=14)
    plt.title('Charge neutrality', fontsize=16)
    plt.suptitle(title)
    plt.tight_layout()
    
    # Saving the results
    plt.savefig(save_path + title + '.png', dpi=600, bbox_inches='tight')
    plt.show()

def plot_results_flexible(odeint_solution, x_axis, value, axes_row, legendloc):
    """
    Plots the three main graphs of a certain stage. The results are not saved, but displayed in the given axes_row. Prepared for the comparison of more than one model in the same figure.
    
    The plots are:
    - n_i evolution
    - Recombination
    - Charge neutrality ratio
    """
    # Unpack the variables
    n_I, n_II, n_III ,n_IV ,n_V ,n_s ,m_R ,m_NR ,n_c , n_v = odeint_solution.T
    
    # Intensity of the glow curve (eq. 3.7)
    dm_R = m_R * value.A_mn_R * n_c
    dm_NR = m_NR * value.A_mn_NR * n_c
    odeint_solution = np.column_stack((odeint_solution, dm_R, dm_NR))
    
    # Axes row
    ax1, ax2, ax3 = axes_row
    
    # Plotting
    ax1.plot(x_axis, n_I, label='n$_{I}$(t)')
    ax1.plot(x_axis, n_II, label='n$_{II}$(t)')
    ax1.plot(x_axis, n_III, label='n$_{III}$(t)')
    ax1.plot(x_axis, n_IV, label='n$_{IV}$(t)')
    ax1.plot(x_axis, n_V, label='n$_{V}$(t)')
    ax1.plot(x_axis, n_s, label='n$_{s}$(t)')
    ax1.set_xlabel('Time (s)', fontsize=14)
    ax1.set_ylabel('Trap concentration (cm$^{-3}$)', fontsize=14)
    ax1.set_title('n$_{i}$ evolution', fontsize=16)
    ax1.legend(loc = legendloc)

    ax2.plot(x_axis, dm_R, label='dm$_{R}$(t)')
    ax2.plot(x_axis, dm_NR, label='dm$_{NR}$(t)')
    ax2.set_xlabel('Time (s)', fontsize=14)
    ax2.set_ylabel('Intensity [u.a.]', fontsize=14)
    ax2.set_title('Recombination', fontsize=16)
    ax2.legend(loc = legendloc)

    # test
    denominator = m_R + m_NR
    c_neutrality = (n_c + n_I + n_II + n_III + n_IV + n_V + n_s) / denominator
    ax3.plot(x_axis, c_neutrality)
    ax3.set_ylim(0.95, 1.05)
    ax3.set_xlabel('Time (s)', fontsize=14)
    ax3.set_title('Charge neutrality ratio', fontsize=16)

def plot_glowcurve_flexible(odeint_solution, x_axis, value, axes_row):
    """
    Plots the three main graphs for the heating stage. The results are not saved, but displayed in the given axes_row. Prepared for the comparison of more than one model in the same figure.
    The plots are:
    - n_i evolution
    - Recombination
    - Charge neutrality ratio
    """
    
    # Unpack the variables
    n_I, n_II, n_III ,n_IV ,n_V ,n_s ,m_R ,m_NR ,n_c , n_v = odeint_solution.T
    ni_curves = np.array([n_I, n_II, n_III, n_IV, n_V, n_s]).T
    
    # Intensity of the glow curve (eq. 3.7)
    dm_R = m_R * value.A_mn_R * n_c
    dm_NR = m_NR * value.A_mn_NR * n_c
    odeint_solution = np.column_stack((odeint_solution, dm_R, dm_NR))
    
    # Axes row
    ax1, ax2, ax3 = axes_row
    
    ## Plotting
    # n_i evolution
    trap_labels = ['n$_I(t)$', 'n$_{II}(t)$', 'n$_{III}(t)$', 'n$_{IV}(t)$', 'n$_V(t)$']
    colors = ['tab:blue', 'tab:orange', 'tab:green', 'tab:red', 'tab:purple', 'tab:brown']
    ax1.plot(x_axis, ni_curves[:, 5], label='n$_{s}$(t)', color=colors[5])  # n_s curve
    for i in reversed(range(5)):
        n = ni_curves[:, i]
        dn_dT = -np.gradient(n, x_axis)  # negative gradient

        # 1. Glow peak (max detrapping rate)
        T_peak = x_axis[np.argmax(dn_dT)]
        ax1.plot(T_peak, n[np.argmax(dn_dT)], '*', color=colors[i], markersize=10)
        ax1.plot([T_peak, T_peak], [0, n[np.argmax(dn_dT)]], linestyle='dotted', color=colors[i])
        ax1.text(T_peak, -0.05e6, f'{T_peak:.0f}°C', ha='center', va='top', fontsize=9, color=colors[i])
        
        # 2. Activation temperature (start of decrease)
        start_idx = np.argmax(dn_dT > 1e-2)  # threshold to detect slope change
        T_start = x_axis[start_idx]
        ax1.plot(T_start, n[start_idx], 'o', color=colors[i], markersize=6)
        
        # 3. Plot n_i(T)
        ax1.plot(x_axis, n, label=trap_labels[i], color=colors[i])


    ax1.set_xlabel('Temperature (ºC)', fontsize=14)
    ax1.set_ylabel('Trap concentration (cm$^{-3}$)', fontsize=14)
    ax1.set_title('n$_{i}$ evolution', fontsize=16)
    ax1.legend(reverse=True, loc = 'upper left', fontsize=10)

    # Recombination
    peaks, _ = find_peaks(dm_R, prominence=50)
    T_peaks = x_axis[peaks]
    I_peaks = dm_R[peaks]
    
    ax2.plot(x_axis, dm_R, label='dm$_{R}$(t)')
    ax2.plot(x_axis, dm_NR, label='dm$_{NR}$(t)')
    for t_peak, i_peak, i in zip(T_peaks, I_peaks, range(len(T_peaks))):
        ax2.plot(t_peak, i_peak, '*', markersize=10, color=colors[i])
        ax2.plot([t_peak, t_peak], [0, i_peak], linestyle='dotted', color=colors[i])
        ax2.plot(t_peak, 0, 'o', color=colors[i])
        ax2.text(t_peak, -0.02 * max(dm_R), f'{t_peak:.0f}°C', ha='center', va='top', fontsize=9, color=colors[i])
    manual_T = 356
    manual_I = dm_R[manual_T]
    ax2.plot(manual_T, manual_I, '*', markersize=10, color=colors[4])
    ax2.plot([manual_T, manual_T], [0, manual_I], linestyle='dotted', color=colors[4])
    ax2.plot(manual_T, 0, 'o', color=colors[4])
    ax2.set_xlabel('Temperature (ºC)', fontsize=14)
    ax2.set_ylabel('Intensity [u.a.]', fontsize=14)
    ax2.set_title('Recombination', fontsize=16)
    ax2.legend(loc = 'upper left')

    # test
    denominator = m_R + m_NR + n_v
    c_neutrality = (n_c + n_I + n_II + n_III + n_IV + n_V + n_s) / denominator
    ax3.plot(x_axis, c_neutrality)
    ax3.set_ylim(0.95, 1.05)
    ax3.set_xlabel('Temperature (ºC)', fontsize=14)
    ax3.set_title('Charge neutrality', fontsize=16)

def plot_column1(odeint_solution1, odeint_solution2, save_path, x_axis, legendloc, xaxis_label, phase):
    """
    Plots the comparison of the n_i evolution for two different models (temperature dependent and temperature independent) in a single figure.The results are saved in the 'Results' folder with the given title.
    """
    # Unpack the variables
    n_I1, n_II1, n_III1 ,n_IV1 ,n_V1 ,n_s1 ,m_R1 ,m_NR1 ,n_c1 , n_v1 = odeint_solution1.T
    n_I2, n_II2, n_III2 ,n_IV2 ,n_V2 ,n_s2 ,m_R2 ,m_NR2 ,n_c2 , n_v2 = odeint_solution2.T
    
    # Plotting
    plt.figure(figsize=(15, 8))
    plt.subplot(1, 2, 1)
    plt.gca().text(0.5, -0.13, '(a)', transform=plt.gca().transAxes, fontsize=12, va='top', ha='center')
    plt.plot(x_axis, n_I1, label='n$_{I}$(t)')
    plt.plot(x_axis, n_II1, label='n$_{II}$(t)')
    plt.plot(x_axis, n_III1, label='n$_{III}$(t)')
    plt.plot(x_axis, n_IV1, label='n$_{IV}$(t)')
    plt.plot(x_axis, n_V1, label='n$_{V}$(t)')
    plt.plot(x_axis, n_s1, label='n$_{s}$(t)')
    plt.xlabel(xaxis_label, fontsize=14)
    plt.ylabel('Trap concentration (cm$^{-3}$)', fontsize=14)
    plt.title('Temperature dependent model', fontsize=16)
    plt.legend(loc = legendloc)
    
    plt.subplot(1, 2, 2)
    plt.gca().text(0.5, -0.13, '(b)', transform=plt.gca().transAxes, fontsize=12, va='top', ha='center')
    plt.plot(x_axis, n_I2, label='n$_{I}$(t)')
    plt.plot(x_axis, n_II2, label='n$_{II}$(t)')
    plt.plot(x_axis, n_III2, label='n$_{III}$(t)')
    plt.plot(x_axis, n_IV2, label='n$_{IV}$(t)')
    plt.plot(x_axis, n_V2, label='n$_{V}$(t)')
    plt.plot(x_axis, n_s2, label='n$_{s}$(t)')
    
    plt.xlabel(xaxis_label, fontsize=14)
    plt.ylabel('Trap concentration (cm$^{-3}$)', fontsize=14)
    plt.title('Temperature independent model', fontsize=16)
    plt.legend(loc = legendloc)
    plt.suptitle(r'n$_{i}$ evolution for the ' + phase, fontsize=18)
    plt.tight_layout()
    plt.savefig(save_path + 'n_i evolution' + '.png', dpi=600, bbox_inches='tight')
    plt.show()
    
def plot_column2(odeint_solution1, odeint_solution2, save_path, x_axis, value, legendloc, xaxis_label, phase):
    """
    Plots the comparison of the recombination for two different models (temperature dependent and temperature independent) in a single figure. The results are saved in the 'Results' folder with the given title.
    """
    # Unpack the variables
    n_I1, n_II1, n_III1 ,n_IV1 ,n_V1 ,n_s1 ,m_R1 ,m_NR1 ,n_c1 , n_v1 = odeint_solution1.T
    n_I2, n_II2, n_III2 ,n_IV2 ,n_V2 ,n_s2 ,m_R2 ,m_NR2 ,n_c2 , n_v2 = odeint_solution2.T
    
    # Intensity of the glow curve (eq. 3.7)
    dm_R1 = m_R1 * value.A_mn_R * n_c1
    dm_NR1 = m_NR1 * value.A_mn_NR * n_c1
    odeint_solution1 = np.column_stack((odeint_solution1, dm_R1, dm_NR1))
    
    dm_R2 = m_R2 * value.A_mn_R * n_c2
    dm_NR2 = m_NR2 * value.A_mn_NR * n_c2
    odeint_solution2 = np.column_stack((odeint_solution2, dm_R2, dm_NR2))
    
    # Plotting
    plt.figure(figsize=(15, 8))
    
    plt.subplot(1, 2, 1)
    plt.gca().text(0.5, -0.13, '(a)', transform=plt.gca().transAxes, fontsize=12, va='top', ha='center')
    plt.plot(x_axis, dm_R1, label='dm$_{R}$(t)')
    plt.plot(x_axis, dm_NR1, label='dm$_{NR}$(t)')
    plt.xlabel(xaxis_label, fontsize=14)
    plt.ylabel('Intensity [u.a.]', fontsize=14)
    plt.title('Temperature dependent model', fontsize=16)
    plt.legend(loc = legendloc)

    plt.subplot(1, 2, 2)
    plt.gca().text(0.5, -0.13, '(b)', transform=plt.gca().transAxes, fontsize=12, va='top', ha='center')
    plt.plot(x_axis, dm_R2, label='dm$_{R}$(t)')
    plt.plot(x_axis, dm_NR2, label='dm$_{NR}$(t)')
    plt.xlabel(xaxis_label, fontsize=14)
    plt.ylabel('Intensity [u.a.]', fontsize=14)
    plt.title('Temperature independent model', fontsize=16)
    plt.legend(loc = legendloc)
    
    plt.suptitle('Recombination for the ' + phase, fontsize=18)
    plt.tight_layout()
    
    # Saving the results
    plt.savefig(save_path + 'Recombination' + '.png', dpi=600, bbox_inches='tight')
    plt.show()

def plot_column3(odeint_solution1, odeint_solution2, save_path, x_axis, value, xaxis_label, phase):
    """
    Plots the comparison of the charge neutrality ratio for two different models (temperature dependent and temperature independent) in a single figure. The results are saved in the 'Results' folder with the given title.
    """
    # Unpack the variables
    n_I1, n_II1, n_III1 ,n_IV1 ,n_V1 ,n_s1 ,m_R1 ,m_NR1 ,n_c1 , n_v1 = odeint_solution1.T
    n_I2, n_II2, n_III2 ,n_IV2 ,n_V2 ,n_s2 ,m_R2 ,m_NR2 ,n_c2 , n_v2 = odeint_solution2.T
    
    denominator = m_R1 + m_NR1
    c_neutrality1 = (n_c1 + n_I1 + n_II1 + n_III1 + n_IV1 + n_V1 + n_s1)/denominator

    denominator = m_R2 + m_NR2
    c_neutrality2 = (n_c2 + n_I2 + n_II2 + n_III2 + n_IV2 + n_V2 + n_s2)/denominator
    
    # Plotting
    plt.figure(figsize=(15, 6))
    
    plt.subplot(1, 2, 1)
    plt.gca().text(0.5, -0.13, '(a)', transform=plt.gca().transAxes, fontsize=12, va='top', ha='center')
    plt.plot(x_axis, c_neutrality1)
    plt.xlabel(xaxis_label, fontsize=14)
    plt.ylabel('Charge neutrality ratio', fontsize=14)
    plt.ylim(0.95, 1.05)

    plt.title('Temperature dependent model', fontsize=16)

    plt.subplot(1, 2, 2)
    plt.gca().text(0.5, -0.13, '(b)', transform=plt.gca().transAxes, fontsize=12, va='top', ha='center')
    plt.plot(x_axis, c_neutrality2)
    plt.xlabel(xaxis_label, fontsize=14)
    plt.ylabel('Charge neutrality ratio', fontsize=14)
    plt.title('Temperature independent model', fontsize=16)
    plt.ylim(0.95, 1.05)
    plt.suptitle('Charge neutrality for the ' + phase, fontsize=18)
    
    # Saving the results
    plt.savefig(save_path + 'Charge neutrality' + '.png', dpi=600, bbox_inches='tight')
    plt.show()
