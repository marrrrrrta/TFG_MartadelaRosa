# Plots and saves the results to the 'Results' folder
import numpy as np
import matplotlib.pyplot as plt

def plot_results(odeint_solution, save_path, title, x_axis, value):
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
    ax2.set_ylabel(r'$\frac{dm_{i=R,NR}}{dt}$ [u.a.]', fontsize=14)
    ax2.set_title('Recombination', fontsize=16)
    ax2.legend(loc = legendloc)

    # test
    denominator = m_R + m_NR
    denominator[denominator == 0] = np.nan
    ax3.plot(x_axis, (n_c + n_I + n_II + n_III + n_IV + n_V + n_s)/denominator)
    ax3.set_xlabel('Time (s)', fontsize=14)
    ax3.set_title('Charge neutrality', fontsize=16)

def plot_glowcurve_flexible(odeint_solution, x_axis, value, axes_row):
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
    ax1.set_xlabel('Temperature (ºC)', fontsize=14)
    ax1.set_ylabel('Trap concentration (cm$^{-3}$)', fontsize=14)
    ax1.set_title('n$_{i}$ evolution', fontsize=16)
    ax1.legend(loc = 'upper left')

    ax2.plot(x_axis, dm_R, label='dm$_{R}$(t)')
    ax2.plot(x_axis, dm_NR, label='dm$_{NR}$(t)')
    ax2.set_xlabel('Temperature (ºC)', fontsize=14)
    ax2.set_ylabel(r'$\frac{dm_{i=R,NR}}{dt}$ [u.a.]', fontsize=14)
    ax2.set_title('Recombination', fontsize=16)
    ax2.legend(loc = 'upper left')

    # test
    denominator = m_R + m_NR + n_v
    numerator = n_c + n_I + n_II + n_III + n_IV + n_V + n_s
    charge_ratio = numerator / denominator
    charge_ratio[np.isnan(charge_ratio) | np.isinf(charge_ratio)] = np.nan  # hide bad points
    ax3.plot(x_axis, charge_ratio)
    ax3.set_xlabel('Temperature (ºC)', fontsize=14)
    ax3.set_title('Charge neutrality', fontsize=16)

def plot_column1(odeint_solution1, odeint_solution2, save_path, x_axis, legendloc, xaxis_label, phase):
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
    # Unpack the variables
    n_I1, n_II1, n_III1 ,n_IV1 ,n_V1 ,n_s1 ,m_R1 ,m_NR1 ,n_c1 , n_v1 = odeint_solution1.T
    n_I2, n_II2, n_III2 ,n_IV2 ,n_V2 ,n_s2 ,m_R2 ,m_NR2 ,n_c2 , n_v2 = odeint_solution2.T
    
    denominator = m_R1 + m_NR1 + n_v1
    denominator[denominator == 0] = np.nan
    c_neutrality1 = (n_c1 + n_I1 + n_II1 + n_III1 + n_IV1 + n_V1 + n_s1)/denominator
    c_neutrality1 = c_neutrality1[:-50]
    x_axis = x_axis[:-50]
    
    denominator = m_R2 + m_NR2 + n_v2
    denominator[denominator == 0] = np.nan
    c_neutrality2 = (n_c2 + n_I2 + n_II2 + n_III2 + n_IV2 + n_V2 + n_s2)/denominator
    c_neutrality2 = c_neutrality2[:-50]
    
    # Plotting
    plt.figure(figsize=(15, 6))
    
    plt.subplot(1, 2, 1)
    plt.gca().text(0.5, -0.13, '(a)', transform=plt.gca().transAxes, fontsize=12, va='top', ha='center')
    plt.plot(x_axis, c_neutrality1)
    plt.xlabel(xaxis_label, fontsize=14)
    plt.ylabel('Charge neutrality ratio', fontsize=14)
    plt.ylim(0.5, 10)
    plt.title('Temperature dependent model', fontsize=16)

    plt.subplot(1, 2, 2)
    plt.gca().text(0.5, -0.13, '(b)', transform=plt.gca().transAxes, fontsize=12, va='top', ha='center')
    plt.plot(x_axis, c_neutrality2)
    plt.xlabel(xaxis_label, fontsize=14)
    plt.ylabel('Charge neutrality ratio', fontsize=14)
    plt.ylim(0.5, 10)
    plt.title('Temperature independent model', fontsize=16)
    plt.suptitle('Charge neutrality for the ' + phase, fontsize=18)
    
    # Saving the results
    plt.savefig(save_path + 'Charge neutrality' + '.png', dpi=600, bbox_inches='tight')
    plt.show()
    