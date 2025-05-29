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
    plt.figure(figsize=(15, 6))
    
    plt.subplot(1, 3, 1)
    plt.plot(x_axis, n_I, label='n$_{I}$(t)')
    plt.plot(x_axis, n_II, label='n$_{II}$(t)')
    plt.plot(x_axis, n_III, label='n$_{III}$(t)')
    plt.plot(x_axis, n_IV, label='n$_{IV}$(t)')
    plt.plot(x_axis, n_V, label='n$_{V}$(t)')
    plt.plot(x_axis, n_s, label='n$_{s}$(t)')
    plt.xlabel('Time (s)')
    plt.ylabel('Trap concentration (cm$^{-3}$)')
    plt.title('n$_{i}$ evolution')
    plt.legend(loc = 'upper left')

    plt.subplot(1, 3, 2)
    plt.plot(x_axis, dm_R, label='dm$_{R}$(t)')
    plt.plot(x_axis, dm_NR, label='dm$_{NR}$(t)')
    plt.xlabel('Time (s)')
    plt.ylabel(r'$\frac{dm_{i=R,NR}}{dt}$ [u.a.]')
    plt.title('Recombination')
    plt.legend(loc = 'upper left')

    plt.subplot(1, 3, 3)
    plt.plot(x_axis, (n_c + n_I + n_II + n_III + n_IV + n_V + n_s)/(m_R + m_NR))
    plt.xlabel('Time (s)')
    plt.title('Charge neutrality')
    plt.suptitle(title)
    plt.tight_layout()
    
    # Saving the results
    plt.savefig(save_path + title + '.png', dpi=300, bbox_inches='tight')
    #plt.show()
    
def plot_glowcurve(odeint_solution, save_path, title, x_axis, value):
    # Unpack the variables
    n_I, n_II, n_III ,n_IV ,n_V ,n_s ,m_R ,m_NR ,n_c , n_v = odeint_solution.T
    
    # Intensity of the glow curve (eq. 3.7)
    dm_R = m_R * value.A_mn_R * n_c
    dm_NR = m_NR * value.A_mn_NR * n_c
    odeint_solution = np.column_stack((odeint_solution, dm_R, dm_NR))
    
    # Plotting
    plt.figure(figsize=(15, 6))
    
    plt.subplot(1, 3, 1)
    plt.plot(x_axis, n_I, label='n$_{I}$(t)')
    plt.plot(x_axis, n_II, label='n$_{II}$(t)')
    plt.plot(x_axis, n_III, label='n$_{III}$(t)')
    plt.plot(x_axis, n_IV, label='n$_{IV}$(t)')
    plt.plot(x_axis, n_V, label='n$_{V}$(t)')
    plt.plot(x_axis, n_s, label='n$_{s}$(t)')
    plt.xlabel('Temperature (ºC)')
    plt.ylabel('Trap concentration (cm$^{-3}$)')
    plt.title('n$_{i}$ evolution')
    plt.legend(loc = 'upper left')

    plt.subplot(1, 3, 2)
    plt.plot(x_axis, dm_R, label='dm$_{R}$(t)')
    plt.plot(x_axis, dm_NR, label='dm$_{NR}$(t)')
    plt.xlabel('Temperature (ºC)')
    plt.ylabel('I [u.a.]')
    plt.title('Recombination')
    plt.legend(loc = 'upper left')

    plt.subplot(1, 3, 3)
    plt.plot(x_axis, (n_c + n_I + n_II + n_III + n_IV + n_V + n_s)/(m_R + m_NR + n_v))
    plt.xlabel('Temperature (ºC)')
    plt.title('Charge neutrality')
    plt.suptitle(title)
    plt.tight_layout()
    
    # Saving the results
    plt.savefig(save_path + title + '.png', dpi=300, bbox_inches='tight')
    plt.show()