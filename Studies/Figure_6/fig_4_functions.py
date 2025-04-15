# -*- coding: utf-8 -*-
"""
Run bacstroke and generate multiple data sets from two given folders containing
initial condition files and configuration files. 
"""
# imports 

import os
from BacStroke import main as bs
import BacPlot as bp
import Functions as f
import glob
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import matplotlib as mpl


def rotated_plot(config_path, data, output_path, save = False):
    '''
    Plots raw coordinates side by side with the rotated version of the coordinates.
    '''
    
    # get constants related to the data file
    constants = f.read_config(config_path)
    
    # rotation rate
    RPM = float(constants[3]) # rotation rate in RPM
    omega = RPM*2*np.pi/60# rotation rate in rad/s
    print(omega)
    
    # clinostat dimensions
    R = float(constants[4])
    r = float(constants[5])
    H = float(constants[6])
    
    # read in output from bacstroke
    df = np.array(pd.read_csv(data, sep=',', header = None))
    
    x_coords = df[:, 0]
    y_coords = df[:, 1]
    z_coords = df[:, 2]
    time = df[:, 3]
    
    # angle to rotate positions by
    theta = -1*omega*time # -1 for a negative rotation
    print(theta)
    
    costheta = np.cos(theta)
    sintheta = np.sin(theta)
    
    # rotated x and y positions, z is not effected by the rotation.
    rotx = (x_coords*costheta) - (y_coords*sintheta)
    roty = (x_coords*sintheta) + (y_coords*costheta)  
    
    fig, ax = plt.subplots(1, 2, figsize = (20, 9), gridspec_kw={'width_ratios': [1, 1]})
    
    ax[0].scatter(x_coords, y_coords, zorder = 20, c=time, cmap = 'coolwarm')
    
    # outer circle patch
    cir = plt.Circle((0, 0), R, facecolor='#c7c7c7', alpha=1, linewidth=3, linestyle='--', edgecolor='black')#color='darkorange',fill=False)
    ax[0].add_patch(cir)
    
    # inner circle patch
    cir2 = plt.Circle((0, 0), r, facecolor='white', alpha=1, linewidth=3, linestyle='--', edgecolor='black')#color='darkorange',fill=False)
    ax[0].add_patch(cir2)
      
    scatter0 = ax[1].scatter(rotx, roty, zorder = 20, c=time, cmap = 'coolwarm')
    ax[1].plot(rotx[0], roty[0], zorder = 100, color = 'pink')
    
    # outer circle patch
    cir = plt.Circle((0, 0), R, facecolor='#c7c7c7', alpha=1, linewidth=3, linestyle='--', edgecolor='black')#color='darkorange',fill=False)
    ax[1].add_patch(cir)
    
    # inner circle patch
    cir2 = plt.Circle((0, 0), r, facecolor='white', alpha=1, linewidth=3, linestyle='--', edgecolor='black')#color='darkorange',fill=False)
    ax[1].add_patch(cir2)
    
    # # ax[1].scatter(rotx, roty, zorder = 20)
    # # ax[1].plot(rotx[0], roty[0], zorder = 100, color = 'pink')
    
    # # Plot the second subplot with a color bar
    # scatter = ax[1].scatter(rotx, roty, zorder=20, c=time, cmap = 'coolwarm')
    # ax[1].plot(rotx[0], roty[0], 'rX', zorder = 1000) # highlight first point, to see direction of movement
    
    # # Add a color bar to the second plot
    cbar = plt.colorbar(scatter0, ax=ax[1], fraction = 0.046, pad = 0.04)
    cbar.set_label('Time, t (s)')
    
    ax[0].set_title('Lab Frame', fontsize = 40)
    ax[1].set_title('Rotated Frame', fontsize = 40)    
    
    plt.show()

def run_sims(config_folder, no_runs, output_folder):
    '''
    Runs simulations for multiple configs

    Parameters
    ----------
    config_folder : STRING
        Absolute path to folder containing config files needed for BacStroke
    no_runs : INTEGER
        Number of runs to be executed for each config file.

    Returns
    -------
    None.

    '''
    
    # grab config file names from config folder
    config_list = os.listdir(config_folder)
    
    # create output folder if it doesnt exist already
    if not os.path.isdir(output_folder): # check if folder exists
            os.makedirs(output_folder) 
    
    for i in range(len(config_list)):# for each config file
    
        config_path = f'{config_folder}/{config_list[i]}'
        config_output_folder = f'{output_folder}/{config_list[i]}' # subfolder within output_folder
        
        # create output subfolder for config within output folder
        if not os.path.isdir(config_output_folder): # check if folder exists
                os.makedirs(config_output_folder) 
                
        # get path to initial conditions file
        config_vals = f.read_config(config_path) # first parameter in config is the ic path
        
        ic_path = config_vals[0] # first parameter in config is the ic path
        
        R = float(config_vals[4]) # outer radius
        r = float(config_vals[5]) # inner radius
        H = float(config_vals[6]) # height
    
        for j in range(no_runs): # run BacStroke for specified number of runs
        
            # randomise starting position of the bacteria
            f.change_starting_coords(ic_path, f.sample_from_hollow_cylinder(r, R, H))  # change this to the dimensions gotten from the config file
        
            run_output = f'{config_output_folder}/run_{j + 1}'
        
            # run BacStroke simulation with randomised starting position
            if not os.path.isfile(run_output): # checks if run has already occured
                bs(config_path, run_output)
            
            
def plot_multi_files(data_file_path, config_file, number_to_plot): 
    '''
    Plot multiple data sets from the same configuration file. (May adapt to 
                                                               save figure)

    Parameters
    ----------
    data_file_path : STRING
        absolute path to folder containing all of the datasets that are to be
        plotted.
    config_file : STRING
        absolute path to the configuration .txt file used to create all of the
        data files inside the data_file_path folder.

    Returns
    -------
    None.

    '''
    
    data = np.array(glob.glob(f'{data_file_path}/*')) # get all paths inside this folder
    
    # read config file
    config_vals = f.read_config(config_file) # first parameter in config is the ic path
    
    R = float(config_vals[4]) # outer radius
    r = float(config_vals[5]) # inner radius
    H = float(config_vals[6]) # height
    
    # set up plot
    fig, ax = plt.subplots(1, 1, figsize = (10, 9)) # set up plots
    
    ax.set_ylabel('y')
    ax.set_xlabel('x')
    
    # outer circle patch
    cir = plt.Circle((0, 0), R, facecolor='#c7c7c7', alpha=1, linewidth=3, linestyle='--', edgecolor='black')
    ax.add_patch(cir)
    
    # inner circle patch
    cir2 = plt.Circle((0, 0), r, facecolor='white', alpha=1, linewidth=3, linestyle='--', edgecolor='black')
    ax.add_patch(cir2)
       
    # loop over the data and plot to the subplot
    for i in range(number_to_plot):
        
        # read data and plot onto the plot
        df = np.array(pd.read_csv(data[i], sep=',', header=None))
        
        x_coords = df[:, 0]
        y_coords = df[:, 1]
        
        ax.scatter(x_coords, y_coords, s = 5)
    
    plt.show()
    

def get_mpl_colours(colour_map, no_colours):
    '''
    Generate no_colours from a matplotlib colour map

    Parameters
    ----------
    cmap : TYPE
        DESCRIPTION.

    Returns
    -------
    None.

    '''
    
    # grab the colour map from matplotlib
    cmap = mpl.colormaps[colour_map]
    
    # get the specified number of colours from matplotlib
    colours = cmap(np.linspace(0, 1, no_colours))
    
    return colours
    
    
    
def plot_multi_files_rotated(data_file_path, config_file, number_to_plot, save_to, cmap): 
    '''
    Plot multiple data sets from the same configuration file. (May adapt to 
                                                               save figure)

    Parameters
    ----------
    data_file_path : STRING
        absolute path to folder containing all of the datasets that are to be
        plotted.
    config_file : STRING
        absolute path to the configuration .txt file used to create all of the
        data files inside the data_file_path folder.

    Returns
    -------
    None.

    '''
    
    data = np.array(glob.glob(f'{data_file_path}/*')) # get all paths inside this folder
    print(data)
    # read config file
    config_vals = f.read_config(config_file) # first parameter in config is the ic path
    
    R = float(config_vals[4]) # outer radius
    r = float(config_vals[5]) # inner radius
    H = float(config_vals[6]) # height
    
    RPM = float(config_vals[3])
    omega = RPM*2*np.pi/60
    print(omega)
    
    # set up plot
    fig, ax = plt.subplots(1, 1, figsize = (9, 9)) # set up plots
    
    # remove ticks
    ax.set_xticks([])
    ax.set_yticks([])
    
    ax.set_facecolor('white')
    
    ax.spines['bottom'].set_color('white')
    ax.spines['top'].set_color('white')
    ax.spines['left'].set_color('white')
    ax.spines['right'].set_color('white')
    
    # ax.set_ylabel('y')
    # ax.set_xlabel('x')
    
    # outer circle patch
    cir = plt.Circle((0, 0), R, facecolor='#c7c7c7', alpha=1, linewidth=3, linestyle='--', edgecolor='black')
    ax.add_patch(cir)
    
    # inner circle patch
    cir2 = plt.Circle((0, 0), r, facecolor='white', alpha=1, linewidth=3, linestyle='--', edgecolor='black')
    ax.add_patch(cir2)
    
    # generate enough colours for plotting
    colours = get_mpl_colours(cmap, number_to_plot)
       
    # loop over the data and plot to the subplot
    for i in range(number_to_plot):
        
        # read data and plot onto the plot
        df = np.array(pd.read_csv(data[i], sep=',', header=None))
        
        x_coords = df[:, 0]
        y_coords = df[:, 1]
        time = df[:, 3]
        
        # angle to rotate positions by
        theta = -1*omega*time # -1 for a negative rotation
        
        costheta = np.cos(theta)
        sintheta = np.sin(theta)
        
        # rotated x and y positions, z is not effected by the rotation.
        rotx = (x_coords*costheta) - (y_coords*sintheta)
        roty = (x_coords*sintheta) + (y_coords*costheta)
        
        ax.scatter(rotx,roty, s = 5, c = colours[i])
    
    plt.savefig(save_to, dpi = 2000)
    plt.show()
    

def plot_frac_at_wall_vs_time(data_file_path, config_file):
    '''
    Generate dataset of fraction of bacteria at the clinostat wall vs time
    
    Parameters
    ----------
    data_file_path: STRING
        Absolute path to folder containing all datasets generated using a 
        single configuration file. These all must contain the same number of
        datapoints and be consistent in time.
        
    config_file: STRING
        Absolute path to the config file that was used to generate the data 
        files within data_file_path
        
    Returns
    -------
    frac_at_wall : NUMPY ARRAY (LENGTH OF A DATA FILE WITHIN DATA_FILE_PATH)
        contains the fraction of bacteria at a wall as a function of time.

    '''
    
    # grab all the data files inside the master data folder
    data = np.array(glob.glob(f'{data_file_path}/*'))
    no_bacteria = len(data)
    
    # read config file to determine how many timesteps are in each file
    config_vals = f.read_config(config_file)
    print(config_vals)
    
    # grab relevent parameters from the config file
    dt = float(config_vals[1])
    total_time = float(config_vals[2])
    output_interval = int(config_vals[-1])
    print(dt, total_time, output_interval)
    
    no_timesteps = int((total_time/dt)/output_interval)
    
    # store count for number of bacteria at wall for each timestep
    bac_at_wall = np.zeros(no_timesteps - 1) # -1 makes sure the shape aligns with the data array shapes
    print(bac_at_wall)
    
    # get size of clinostat from the config file, used to define if a bacteria
    # is at the wall in the for loop
    R = float(config_vals[4])
    R_min = float(config_vals[5]) # inner radius
    print(R, R_min)
    
    RPM = float(config_vals[3])
    omega = RPM*2*np.pi/60
    print(omega)
    
    # read in bacteria size from initial conditions file 
    ic_file = config_vals[0]
    print(ic_file)
    
    # read ic file
    ic_params = f.read_ic(ic_file)
    a = float(ic_params[1]) # bacteria radius

    # loop through all data files
    for i in range(len(data)):
        
        # read datafile
        df = np.array(pd.read_csv(data[i]))
        #print(df.shape)
        
        # rotate coordinates
        # pull coords to determine radial position of bacteria
        x = df[:, 0]
        y = df[:, 1]
        time = df[:, 3]
        # plt.scatter(x, y)
        # plt.show()
        
        # angle to rotate positions by
        theta = -1*omega*time # -1 for a negative rotation
        
        costheta = np.cos(theta)
        sintheta = np.sin(theta)
        
        # rotated x and y positions, z is not effected by the rotation.
        rotx = (x*costheta) - (y*sintheta)
        roty = (x*sintheta) + (y*costheta)
        
        r = np.sqrt(rotx**2 + roty**2) # radial position of bacteria
        
        # checks if bacteria at the wall, 1 = yes, 0 = no 
        at_wall_outer = np.where(r > R - (2*a), 1, 0) # !!! not sure about this 2*a condition, not convinced should discuss with tyler once code finished
        at_wall_inner = np.where(r < R_min + (2*a), 1, 0)
        
        # add the wall count to the master array
        bac_at_wall += (at_wall_outer + at_wall_inner)
        #print(np.sum(bac_at_wall))
     
    # # simulation time to match bac_at_wall
    # time = df[:, 3]
    
    # convert from number to fraction of bacteria at wall
    frac_at_wall = bac_at_wall/no_bacteria # convert from number to fraction of bacteria at wall
        
    return frac_at_wall, time


# generate data for plots
#run_sims('C:/Users/Kenzie/Documents/Education/Edinburgh_University/Summer_Masters/BacStroke/Studies/Figure_4/run_files/configs', 100, 'C:/Users/Kenzie/Documents/Education/Edinburgh_University/Summer_Masters/BacStroke/Studies/Figure_4/data')


a, b = plot_frac_at_wall_vs_time('C:/Users/Kenzie/Documents/Education/Edinburgh_University/Summer_Masters/BacStroke/Studies/Figure_4/data/omega=0.0001_old', 'C:/Users/Kenzie/Documents/Education/Edinburgh_University/Summer_Masters/BacStroke/Studies/Figure_4/run_files/configs/omega=0.0001_old')

#c, d = plot_frac_at_wall_vs_time('C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/data_output/omega=0.0001', 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/run_files/configs/omega=0.0001')

#e, g = plot_frac_at_wall_vs_time('C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/data_output/omega=0.00001', 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/run_files/configs/omega=0.00001')

#print(a, b)   

# plt.scatter(b, a, label = '0.0001', zorder = 1)  
# # plt.scatter(g, e, label = '$10^{-4}$', zorder = 2)
# # plt.scatter(d, c, label = '$10^{-3}$', zorder = 3)

# plt.xlabel('Time, t (s)')
# plt.ylabel('Fraction at wall')
# plt.legend(title = '$\omega$ (rad/s)')

# t_g = (1e-2)/(5.88e-8)


# plt.axvline(x=t_g, color = 'black', zorder = 4)#, color='b', label='axvline - full height')
# plt.show()
    
     
    
    
    
# plot_multi_files('C:/Users/Kenzie/Documents/Education/Edinburgh_University/Summer_Masters/BacStroke/Studies/Figure_4/data/omega=0.0001_old', 'C:/Users/Kenzie/Documents/Education/Edinburgh_University/Summer_Masters/BacStroke/Studies/Figure_4/run_files/configs/omega=0.0001_old', 20)
# plot_multi_files_rotated('C:/Users/Kenzie/Documents/Education/Edinburgh_University/Summer_Masters/BacStroke/Studies/Figure_4/data/omega=0.0001_old', 'C:/Users/Kenzie/Documents/Education/Edinburgh_University/Summer_Masters/BacStroke/Studies/Figure_4/run_files/configs/omega=0.0001_old', 20, 'test.png', 'Blues')
# # plot_multi_files_rotated('C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/data_output/omega=0.00001', 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/run_files/configs/omega=0.00001', 50, 'C:/Users/kenzi/Documents/Masters/Summer/Figures/Figure 4/omega=0.00001_rotated.png', 'Pinks')
# plot_multi_files_rotated('C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/data_output/omega=0', 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/run_files/configs/omega=0', 50, 'C:/Users/kenzi/Documents/Masters/Summer/Figures/Figure 4/omega=0.png', 'Blues')
# # rotated_plot('C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/run_files/configs/omega=0.00001', 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/data_output/omega=0.00001/run_1', '', save = False)

# # NO ROTATION
#run_sims('C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/run_files/configs', 100, 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/data_output')


#plot_frac_at_wall_vs_time('C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/data_output/omega=0', 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/run_files/configs/omega=0')

# configs = np.array(glob.glob("C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/run_files/configs/*"))
# data = np.array(glob.glob('C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/data_output/omega=0/*'))

# bp.plot_multi_trajectories(configs, data, 'C:/Users/kenzi/Documents/Masters/Summer/Studies/Figure_4/fig.png')


