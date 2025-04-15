# -*- coding: utf-8 -*-
"""
Launch single simulation for figure 6.

To run launch this script in a python editor or the command line. You will be 
prompted for two arguments

1. Absolute path to config file
2. Output folder for sim data
    - This folder is the master data folder and you do not need to create 
      subdirectories for each parameter configuration, this is done automatically.
      
Command line input:

Absolute_path_to_config_file Absolute_path_to_master_data_folder      


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



def get_ic_path(vs = '0'):
    '''
    Change the initial conditions path in the config file

    Returns
    -------
    None.

    '''
    
    parent_directory = Path(__file__).resolve().parent
    
    # get path to initial conditions file
    run_path = os.path.join(parent_directory, 'run_files', 'initial_conditions', f'vs={vs}.txt')
    
    return run_path  


def find_between(s, start, end):
    start_idx = s.find(start) + len(start)
    end_idx = s.find(end, start_idx)
    return s[start_idx:end_idx]

                    
def run_sim(config, output_folder):
    '''
    Run a single simulation of BacStroke. 
    
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
    
    # read config vals 
    vals = f.read_config(config)
    
    ic_path = vals[0]
    print(ic_path)
    
    # nums = [n for n in ic_path if n.isdigit()]
    # print(nums)
    
    swimming_speed = find_between(ic_path, 'vs=', '.txt')   
    print(swimming_speed)
    
    
    # swimming_speed = f'{nums[-2]}{nums[-1]}' # microns per second
    
    # replace path to initial conditions file in config with system specific path
    updated_ic_path = get_ic_path(vs = swimming_speed)
    f.edit_config(config, 0, updated_ic_path)
    
    # get rotation rate and clinostat size from config file
    rpm = float(vals[3]) # rotation rate in rpm
    omega = rpm*2*np.pi/60 # rotation rate in rad/s
    
    # clinostat measurements for coordinate randomisation
    R = float(vals[4]) # outer clinostat radius [m]
    r = float(vals[5]) # inner radius [m]
    H = float(vals[6]) # length/ height [m]
    
    g = float(vals[8]) # gravitational acceleration [m/s^2]
    
    # create an output folder inside output_folder for specific rotation rate
    rounded_omega = f.round_to_first_non_zero(omega) # for decoration purposes
    output_path = os.path.join(output_folder, f'omega={rounded_omega}', f'vs={swimming_speed}', f'g={g}')
    
    # create output folder if it doesnt exist already
    if not os.path.isdir(output_path): # check if folder exists
            os.makedirs(output_path)     

    # count how many sims have been carried out already
    num_sims = len(os.listdir(output_path))
    # print(num_sims)
    
    # establish file name for data output
    output_file = os.path.join(output_path, f'run_{num_sims+1}')
    # print(output_file)
    
    # randomise the starting coordinate
    f.change_starting_coords(updated_ic_path, f.sample_from_hollow_cylinder(r, R, H))
    
    # run simulation
    print('Simulation params')
    print(f'vs = {swimming_speed} µm/s')
    print(f'g = {g} m/s²')
    print(f'ω = {rounded_omega} rad/s')
    bs(config, output_file)
                
    # if not os.path.isfile(run_output):
                
    #     # randomise simulation starting coordinate
    #     f.change_starting_coords(ic, f.sample_from_hollow_cylinder(r, R, H))  # change this to the dimensions gotten from the config file
                        
    #     # begin simulations using the edited config file
    #     bs(config, run_output)
    
config, output_folder = input("Abs path to config folder & output folder: ").split()
    
run_sim(config, output_folder)

