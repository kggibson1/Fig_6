# -*- coding: utf-8 -*-
"""
Generate data for figure 6
"""

import os
from BacStroke import main as bs
import BacPlot as bp
import Functions as f
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl
from pathlib import Path

def run_sims_6(config_folder, ic_folder, output_folder, no_runs):
    '''
    Generate data for figure 6. 
    
    Runs simulations in combination of all configs within the config folder
    and all inital condition files inside the ic_folder, repeating the simulations
    no_run times.
    
    The starting position of the bacteria are randomised at the beginning of each simulation

    Parameters
    ----------
    config_folder : STRING
        Path to the config folder containing each config file.
    output_folder : STRING
        Path to folder that subdata folder will be created and data written into.

    Returns
    -------
    None.

    '''
    
    
    # first create folder with same name as config folder, inside the output folder, it it doesnt exist
    
    # get all of the config files inside of the config folder
    configs = glob.glob(f'{config_folder}/*.txt')
    print(configs)
    
    # get all of the initial condition files
    ics = glob.glob(f'{ic_folder}/*.txt')
    print(ics)
    
    # loop over configuration files inside config_folder
    for i in range(len(configs)):
        
        # read the config file
        vals = f.read_config(configs[i])
        print(vals)
        
        rpm = float(vals[3]) # rotation rate in rpm
        omega = rpm*2*np.pi/60 # rotation rate in rad/s
        
        # clinostat measurements for coordinate randomisation
        R = float(vals[4]) # outer clinostat radius [m]
        r = float(vals[5]) # inner radius [m]
        H = float(vals[6]) # length/ height [m]
        
        g = float(vals[8]) # gravitational acceleration [m/s^2]
        
        # create an output folder inside output_folder for specific rotation rate
        rounded_omega = f.round_to_first_non_zero(omega) # for decoration purposes
        output_path = f'{output_folder}/omega={rounded_omega}'
        
        # create output folder if it doesnt exist already
        if not os.path.isdir(output_path): # check if folder exists
                os.makedirs(output_path) 
         
        # loop over initial condition files inside ic_folder        
        for j in range(len(ics)):
            
            # read the swimming speed
            ic_vals = f.read_ic(ics[j])
            vs = ic_vals[-2] # swimming speed [m/s]
            
            # create output folder for the g, vs combination inside output_path, if it doesnt already exist
            g_vs_output_path = f'{output_path}/g={g}/vs={vs}'
            
            if not os.path.isdir(g_vs_output_path): # check if folder exists
                os.makedirs(g_vs_output_path) 
            
            # change the initial conditions path in config to current ic
            f.edit_config(configs[i], 0, ics[j])
            
            for k in range(no_runs):
                
                # check if simulation has already been carried out
                run_output = f'{g_vs_output_path}/run_{k+1}.csv'
                
                if not os.path.isfile(run_output):
                
                    # randomise simulation starting coordinate
                    f.change_starting_coords(ics[j], f.sample_from_hollow_cylinder(r, R, H))  # change this to the dimensions gotten from the config file
                    
                    # begin simulations using the edited config file
                    bs(configs[i], run_output)


parent_directory = Path(__file__).resolve().parent

config_files = f'{parent_directory}/run_files/configs'
ic_files = f'{parent_directory}/run_files/initial_conditions'
data_files = f'{parent_directory}/data'

# run sims for omega = 0.0 rad/s
omega1_config_path = f'{config_files}/omega=0'
run_sims_6(omega1_config_path, ic_files, data_files, 20)
    
# run sims for omega = 0.0001 rad/s
omega1_config_path = f'{config_files}/omega=0.0001'
run_sims_6(omega1_config_path, ic_files, data_files, 20)    

# run sims for omega = 0.0000001 rad/s
omega2_config_path = f'{config_files}/omega=0.0000001'
run_sims_6(omega2_config_path, ic_files, data_files, 20)  

# run sims for omega = 32 rad/s
omega3_config_path = f'{config_files}/omega=32'
#run_sims_6(omega3_config_path, ic_files, data_files, 20) 


