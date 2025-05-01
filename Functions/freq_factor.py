def freq_factor(deltaS):
    """
    Computes the frequency factor that follows the expression:
                 S = nu * K * exp(DeltaS / kB)
    
    * nu: lattive phonon vibration frequency (s-1)  
    * K: transition probability constant
    * kB: Boltzmann constant (eV/K)

    Args:
        deltaS (_type_): entropy change
    """
    kB = 8.6173335e-5          # Boltzmann constant (eV/K)
    K = 9.39e+02                # from 'other_plots.ipynb'
    nu_m = 90000                # from B.H. Bransden book (m-1)
    nu_s = 2.9979e8 * nu_m      # (s-1)
    
    freq_factor = nu_s * K * np.exp(deltaS / kB)
    return freq_factor