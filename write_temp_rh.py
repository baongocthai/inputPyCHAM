import pandas as pd
import numpy as np

def write_temp_rh(file_name, start_time, end_time, time_and_temperature, time_and_Rh):

    f = open(file_name, "a")

    ### Write for Temperature
    # Convert to seconds
    time_and_temperature = time_and_temperature[start_time:end_time]
    date_time_temperature = np.array(time_and_temperature.index)
    time_seconds_temperature = np.divide((np.subtract(date_time_temperature, np.datetime64(start_time))), 1E9)
    air_temp = np.add(np.array(time_and_temperature["Temperature"]), 273.15)
    temperature_array = np.array([time_seconds_temperature,air_temp], dtype = object) 

    # Write
    f.write("temperature = ")
    for i in range(len(temperature_array[1]) - 1):
        f.write("%f, " %(temperature_array[1][i]))
    f.write("%f\n" %(temperature_array[1][len(temperature_array[1]) - 1]))

    f.write("tempt = ")
    for i in range(len(temperature_array[0]) - 1):
        f.write("%d, " %(temperature_array[0][i]))
    f.write("%d\n" %(temperature_array[0][len(temperature_array[0]) - 1]))

    ### Write for Rh
    # Convert to seconds
    time_and_Rh = time_and_Rh[start_time:end_time]
    date_time_Rh = np.array(time_and_Rh.index)
    time_seconds_Rh = np.divide((np.subtract(date_time_Rh, np.datetime64(start_time))), 1E9)
    air_Rh = np.divide(np.array(time_and_Rh["Humidity"]), 100)
    Rh_array = np.array([time_seconds_Rh,air_Rh], dtype = object)

    # Write
    f.write("rh = ")
    for i in range(len(Rh_array[1]) - 1):
        f.write("%f, " %(Rh_array[1][i]))
    f.write("%f\n" %(Rh_array[1][len(Rh_array[1]) - 1]))

    f.write("rht = ")
    for i in range(len(Rh_array[0]) - 1):
        f.write("%d, " %(Rh_array[0][i]))
    f.write("%d\n" %(Rh_array[0][len(Rh_array[0]) - 1]))

    f.close()