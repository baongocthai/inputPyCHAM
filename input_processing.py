import pandas as pd
import numpy as np
import datetime

### Day of Year and Day Time Start ################################################################################

def proc_date_time(start_time):
    # Convert local to (UTC+8) to GMT (UTC+0)
    actual_datetime = start_time - datetime.timedelta(hours=8, minutes=00, seconds=00)
    DayOfYear = int(actual_datetime.strftime('%j')) #no. of day in the year
    daytime_start = actual_datetime.hour*60*60 + actual_datetime.minute*60 #in seconds

    return DayOfYear, daytime_start

### Temperature and Rh ################################################################################

def proc_temp_rh(wb):

    # Setting date and time as index
    wb.index = pd.to_datetime(wb.pop("Date&Time"), format = '%d/%m/%Y %H:%M')

    # Processing for time and temperature
    time_and_temperature = wb[["Temperature"]]
    time_and_temperature = time_and_temperature.dropna()

    # Processing for time and Rh
    time_and_Rh = wb[["Humidity"]]
    time_and_Rh = time_and_Rh.dropna()

    return time_and_temperature, time_and_Rh

### Species Input Names ###############################################################################

def proc_species_inputs_names(wb, species_list, inorganic_list):

    # Processing input names
    wb = wb[wb["Required unit"] == "ppb"]
    mcm = dict(zip(wb["Inputs"], wb["Name in MCM"]))
    mcm_data = {x: mcm[x] for x in mcm if x not in inorganic_list}

    # Double check number of species
    mcm_keys = [x for x in enumerate(mcm) if x in species_list]
    if len(mcm_keys) == 0: print("All species are accounted for")

    return mcm, mcm_data

### Organic Species Concentraion Input ################################################################

def proc_species_inputs(wb, mcm_data):

    # Processing input for every minute
    wb.index = pd.to_datetime(wb.pop("DateTime"),format = '%d/%m/%Y %H:%M')
    # Replace all BD values to 0
    wb.replace({"BD":0},inplace=True)
    wb = wb.astype(float)

    # Data check
    print("All numbers are float" if sum([1 if x=="float64" else 0 for x in wb.dtypes]) == wb.dtypes.count() \
        else "Error: not all numbers are float")
    wb = wb[mcm_data.keys()]

    return wb

### Inorganic Species Concentration Input #############################################################

def proc_inorganic(wb_no2, wb_co, wb_o3, wb_so2, start_time, end_time):
    ### Processing hourly NO2, CO, O3, SO2 
    # Note that CO, O3 and SO2 has to be converted to hourly from NEA average data first
    wb_no2.index = pd.to_datetime(wb_no2.pop("Date&Time"), format='%d/%m/%Y %H:%M')
    wb_no2 = wb_no2.loc[~wb_no2.index.duplicated(keep='first')]

    # co = pd.read_csv("CO_5min.csv",index_col = 0,delimiter=",")
    wb_co.index = pd.to_datetime(wb_co.pop("Date&Time"), format='%d/%m/%Y %H:%M')
    wb_co = wb_co.loc[~wb_co.index.duplicated(keep='first')]

    # o3 = pd.read_csv("O3_5min.csv",index_col = 0,delimiter=",")
    wb_o3.index = pd.to_datetime(wb_o3.pop("Date&Time"), format='%d/%m/%Y %H:%M')
    wb_o3 = wb_o3.loc[~wb_o3.index.duplicated(keep='first')]

    # so2 = pd.read_csv("SO2_5min.csv",index_col = 0,delimiter=",")
    wb_so2.index = pd.to_datetime(wb_so2.pop("DateTime"), format='%d/%m/%Y %H:%M')
    wb_so2 = wb_so2.loc[~wb_so2.index.duplicated(keep='first')]

    ### Convert to ppb
    # Hourly temperature
    hourly_temperature = pd.read_csv("Hourly Meteorology.csv",index_col = 0,delimiter=",")
    hourly_temperature = hourly_temperature[["Date&Time", "Air Temp (oC)"]]
    hourly_temperature.index = pd.to_datetime(hourly_temperature.pop('Date&Time'),format = '%d/%m/%Y %H:%M')

    # Combine inputs
    inorganic = wb_no2.join([wb_co, wb_o3, wb_so2, hourly_temperature], how="outer")
    inorganic = inorganic.rename(columns = {'NO2 (ug/m3)': 'NO2', 'Hourly CO (mg/m3)': \
        'CO', 'Hourly O3 (ug/m3)': 'O3', 'SO2 (ug/m3)': 'SO2'}, inplace = False)
    inorganic = inorganic[start_time:end_time]

    # Convert to ppb
    inorganic["NO2"] = inorganic["NO2"].div(12.187*46.0055/(273.15 + inorganic["Air Temp (oC)"]))
    inorganic["CO"] = inorganic["CO"].div(12.187*28.01/(273.15 + inorganic["Air Temp (oC)"])/1000)
    inorganic["O3"] = inorganic["O3"].div(12.187*48/(273.15 + inorganic["Air Temp (oC)"]))
    inorganic["SO2"] = inorganic["SO2"].div(12.187*64.07/(273.15 + inorganic["Air Temp (oC)"]))

    return inorganic

### Particle Concentration ###########################################################################

def proc_particle(pb, toggle_log, lower_part_size, upper_part_size, number_size_bins, start_time, end_time):

    # Particle processing
    particle = pb
    particle.index = pd.to_datetime(particle.pop("DateTime"), format='%Y-%m-%d %H:%M:%S')
    particle = particle.loc[~particle.index.duplicated(keep='first')]

    # Set the bin spacing
    if toggle_log == True:
        # logarithmic method
        rad_bounds = 10.0**(np.linspace(np.log10(lower_part_size), 
                                np.log10(upper_part_size), num=(number_size_bins+1)))
        rwid = (rad_bounds[1::]-rad_bounds[0:-1]) # width of size bins (um)
        x_output = rad_bounds[0:-1]+rwid/2.0 # particle radius (um)

    else:
        rad_bounds = np.linspace(lower_part_size, upper_part_size, (number_size_bins+1))
    	# width of size bins (um)
        rwid = np.array((rad_bounds[1]-rad_bounds[0])).reshape(1)
        x_output = rad_bounds[0:-1]+rwid/2.0 # particle radius (um)

    # enhance upper radius bound (um) to reduce possibility of particles growing beyond 
    upper_bin_rad_amp = 1.0e6
    rad_bounds[-1] = rad_bounds[-1]*upper_bin_rad_amp

    # Create a deep copy
    conv_particle = particle.copy()
    conv_particle = conv_particle.iloc[:,:0]

    # Classify according to size
    i=0
    while i < (len(rad_bounds)-1):
        temp = np.zeros(len(particle.index))
        lower_bound = rad_bounds[i]
        upper_bound = rad_bounds[i+1]
        for size in particle.columns:
            if (float(size) >= lower_bound and float(size) < upper_bound):
                temp = np.add(temp, particle[size].astype(float))
        conv_particle["%f" %(x_output[i])] = temp
        i=i+1

    # Set index to seconds
    conv_particle = conv_particle.loc[start_time:end_time]
    date_time_part = np.array(conv_particle.index)
    time_seconds_part = (np.divide((np.subtract(date_time_part, np.datetime64(start_time))), 1E9)).astype(int)
    conv_particle = conv_particle.reset_index(drop=True)
    conv_particle = conv_particle.set_index(time_seconds_part)

    return conv_particle

