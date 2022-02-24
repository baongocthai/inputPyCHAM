import os
import pandas as pd
import numpy as np
import input_processing as ip
from write_conc import write_conc
from write_particle import write_particle
from write_temp_rh import write_temp_rh
from write_init import write_init

def main(folder_name, toggle_particle, species_list, inorganic_list, start_time, end_time, toggle_log, 
        lower_part_size, upper_part_size, number_size_bins, file_name, update_step, umansysprop_update, 
        p_init, lat, lon, DayOfYear, daytime_start, light_time, light_status, coag_on, index, const_comp,
        tracked_comp):
    
    # Read file inputs accordingly

    dir_list = os.listdir(folder_name)

    for i in dir_list:
        if ('CO' in i): # CO
            wb_co = pd.read_csv(str(folder_name+'/'+i), index_col = 0,delimiter=",")
        if ('SO2' in i): # SO2
            wb_so2 = pd.read_csv(str(folder_name+'/'+i), index_col = 0,delimiter=",")
        if ('O3' in i): # O3
            wb_o3 = pd.read_csv(str(folder_name+'/'+i), index_col = 0,delimiter=",")
        if ('NO2' in i): # NO2
            wb_no2 = pd.read_csv(str(folder_name+'/'+i), index_col = 0, delimiter=",")
        if ('Meteorology' in i): # temperature and rh
            wb_temp_rh = pd.read_csv(str(folder_name+'/'+i), delimiter=',')
        if ('Minute data_BN3' in i): # gas components
            wb_species = pd.read_csv(str(folder_name+'/'+i))
        if ('Species as inputs' in i): # Species names
            wb_species_inputs = pd.read_csv(str(folder_name+'/'+i))
        if ('SMPS' in i): # particle input
            file_particle = str(folder_name+'/'+i)

    # Read csv file for the Particle Size
    if toggle_particle == True :
        pb = pd.read_csv(file_particle, delimiter = ",", dtype = object)

    # Processing for temperature and Rh
    time_and_temperature, time_and_Rh = ip.proc_temp_rh(wb_temp_rh)

    # Processing input names
    mcm_complete, mcm_data = ip.proc_species_inputs_names(wb_species_inputs, species_list, inorganic_list)

    # Processing concentration input for every minute
    wb_species = ip.proc_species_inputs(wb_species, mcm_data)

    # Processing hourly NO2, CO, O3, SO2 
    inorganic = ip.proc_inorganic(wb_no2, wb_co, wb_o3, wb_so2, start_time, end_time)

    # Combine inorganic and organic species
    result = wb_species[start_time:end_time].join(inorganic, how="outer")

    # Choose which input we are interested in
    result = result[species_list]

    # Particle processing

    particle = ip.proc_particle(pb, toggle_log, lower_part_size, upper_part_size, number_size_bins, start_time, end_time)

    ### Writing the .txt file starts here #####################################################################################################
    ###########################################################################################################################################

    # write items with less/no processing
    write_init(file_name, start_time, end_time, update_step, umansysprop_update, p_init, lat, lon, 
                    DayOfYear, daytime_start, light_time, light_status, coag_on)

    # Write temperature and relative humidity inputs
    write_temp_rh(file_name, start_time, end_time, time_and_temperature, time_and_Rh)

    # Write concentration inputs
    write_conc(file_name, result, start_time, inorganic_list, mcm_complete, index, species_list, const_comp, tracked_comp)

    # Write particle inputs
    write_particle(file_name, toggle_particle, particle, number_size_bins, lower_part_size, upper_part_size, toggle_log)
    
    return ()
