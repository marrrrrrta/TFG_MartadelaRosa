# Luminescent response study of ionic crystals used in dosimetry

This repository contains code and data for simulating and analyzing thermoluminescence (TL) processes in LiF:Mg,Ti, including mathematical modeling, data visualization, and parameter extraction.

## Quick Start
1.  **Install dependencies**   
    Make sure you have Python 3.12+ and the required packages.

2. **Run a simulation**   
    Edit and execute `simulations.py` to perform irradiation, relaxation, and heating simulations.
    Results will be saved in the Results folder.

3. **Explore and plot**   
    Use the provided functions in the Functions folder to analyze the results using the proposed mathematical models. Further analysis can be done in the Jupyter notebook `other_plots.ipynb`.

## Project Structure
This project is organized to allow flexible simulation and analysis of TL processes.
You can easily switch between different mathematical models by changing the `FunctionUsed` variable in simulations.py. Both Cinetic and Structural parameters can be modified in the `ParametrosCineticos.xlsx` and `ParametrosEstructurales.xlsx` files, respectively.

## Contents

- ExperimentalData/
  - `DatEx_TLD_100_Beta_1.xlsx`   
    Experimental TL data.
  - `ParametrosCineticos.xlsx`  
    Kinetic parameters for traps and recombination centers.
  - `ParametrosEstructurales.xlsx`   
    Structural parameters for the simulation.

- Functions/
  - `diff_eqs_freqfactor.py`   
    Differential equations for frequency factor-dependent kinetics.
  - `diff_eqs_growth.py`   
    Differential equations for growth kinetics.
  - `diff_eqs_notemp.py`   
    Differential equations for temperature-independent kinetics.
  - `freq_factor.py`   
    Frequency factor calculations.
  - `plotting.py`   
    Plotting and visualization utilities.

- Results/  
    Output folder for generated plots and figures.

- Notebooks and Scripts
  - `other_plots.ipynb`   
    Jupyter notebook for additional plots and analysis.
  - `simulations.py`   
    Main script for running TL simulations (irradiation, relaxation, heating).
  - `TFG Marta de la Rosa.pdf`   
    Thesis document containing detailed explanations and results.

## About
This project is a comprehensive study of thermoluminescence (TL) processes in ionic crystals, with a focus on LiF:Mg,Ti. The goal is to develop a deeper understanding of the underlying mechanisms and to provide a valuable tool for dosimetry applications.