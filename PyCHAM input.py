import math
import datetime
import pandas as pd
import numpy as np

### Manual Input #############################################################################################
##############################################################################################################

# Set working directory
import os
path = r'C:\Users\janet\Documents\nBox\SOAFP\6. PyCHAM simulation\Input'
os.chdir(path)

# Input time manually
start_time = pd.to_datetime('2020-11-17 09:00:00')
end_time = pd.to_datetime('2020-11-17 15:00:00')

# Input file_name for model variable input
file_name = "run_particle.txt"

###
# EDIT: To import all data from individual csv files for smaller size file
# Enter file location for hourly data and particle data here
#file_location = r"C:\Users\janet\Desktop\working folder\19 Aug 2021 Inputs\PyCHAM inputs.xlsx"
###

# Insert the desired the species
# In total there are 40 overlapping species between MCM inventory and PTR-MS
# 40 VOC species include those with isomers
# Xylene = o-xylene, TMB = TM124B, Pinene = a-pinene, Butene = mepropene, Pentene = trans-2-pentene
# Butanol was not measured in 2011-2012
# 4 inorganic species
species_list = ["Chloromethane", "Dimethyl.Sulfide", "Isoprene", "Dichloro.methane", "Acetone", #5
                "Toluene", "Styrene", "Acetaldehyde","Acetic.acid","Acetic.anhydride", #10
                "Acrolein","Acrylic.acid","Benzene","Butanone","Catechol", #15
                "Cresol","Ethanol","Ethyl.nitrate","Formaldehyde","Formic.acid", #20
                "Hexanal","Methanesulfonic.acid","Methyl.nitrate","MVK","p.Benzoquinone", #25
                "Pentanal","Peroxyacetyl.Nitrate","Phenol","Pinic.acid","Pinonic.acid", #30
                "Propene","Propyl.acetate","Pyruvic.acid","Acetol", #34
                # 6 species with isomers
                "Xylene", "Trimethylbenzene", "Pinene", "Butene", "Pentene", "Butanol", #40
                "NO2", "O3", "SO2", "CO"] # Inorganic species
inorganic_list = ["NO2", "O3", "SO2", "CO"]

## Toggle particle consideration
toggle_particle = True

# Toggle logarithmic scale consideration
toggle_log = True

# Enter number_size_bins
number_size_bins = 5

# Set upper and lower bounds (um)
lower_part_size = 0.0136*1000
upper_part_size = 0.790*1000

#######################################################################################################################
#######################################################################################################################

###
# EDIT: To import all data from individual csv files for smaller size file
# Read Excel File For Concentration Data
#wb = pd.read_excel(file_location, sheet_name = None)
#for key in wb.keys():
#    wb[key].to_csv(('%s.csv' %key))
###

# Read csv file for the Particle Size
if toggle_particle == True :
    file_particle = r"SMPS.csv"
    pb = pd.read_csv(file_particle, delimiter = ",", dtype = object)

# Processing for time and temperature
#time_and_temperature = wb["Meteorology"]
#time_and_temperature = wb["Meteorology"][["Date&Time", "Air Temp (oC)"]]
Meteorology = pd.read_csv("Meteorology.csv",index_col = 0,delimiter=",")
time_and_temperature = Meteorology[["Date&Time", "Air Temp (oC)"]]
time_and_temperature = time_and_temperature.dropna()
time_and_temperature.index = pd.to_datetime(time_and_temperature.pop("Date&Time"))

# Processing for time and Rh
#time_and_Rh = wb["Meteorology"][["Date&Time", "Relative Humidity (%)"]]
time_and_Rh = Meteorology[["Date&Time", "Relative Humidity (%)"]]
time_and_Rh = time_and_Rh.dropna()
time_and_Rh.index = pd.to_datetime(time_and_Rh.pop("Date&Time"))

# Processing input names
#input_names = wb["Species as inputs"]
input_names = pd.read_csv("Species as inputs.csv",index_col = 0,delimiter=",")
input_names = input_names[input_names["Required unit"] == "ppb"]
mcm = dict(zip(input_names["Inputs"], input_names["Name in MCM"]))
mcm_hourly_data = {x: mcm[x] for x in mcm if x not in inorganic_list}

# Double check number of species
mcm_keys = [x for x in enumerate(mcm) if x in species_list]
if len(mcm_keys) == 0: print("All species are accounted for")

# Processing input hourly VOC data
#hourly_data = wb["Hourly data"]
hourly_data = pd.read_csv("Hourly data.csv",index_col = 0,delimiter=",")
hourly_data.index = pd.to_datetime(hourly_data.pop("DateTime"))
hourly_data = hourly_data[mcm_hourly_data.keys()]
hourly_data = hourly_data.loc[~hourly_data.index.duplicated(keep='first')]

# Processing hourly NO2, CO, O3, SO2 
# Note that CO, O3 and SO2 has to be converted to hourly from NEA average data first
#no2 = wb["NO2"]
no2 = pd.read_csv("NO2.csv",index_col = 0,delimiter=",")
no2.index = pd.to_datetime(no2.pop("Date").astype(str) + " " + no2.pop("Time").astype(str))
no2 = no2.loc[~no2.index.duplicated(keep='first')]

#co = wb["CO"]
co = pd.read_csv("CO.csv",index_col = 0,delimiter=",")
co.index = pd.to_datetime(co.pop("Date&Time"))
co = co.loc[~co.index.duplicated(keep='first')]

#o3 = wb["O3"]
o3 = pd.read_csv("O3.csv",index_col = 0,delimiter=",")
o3.index = pd.to_datetime(o3.pop("Date&Time"))
o3 = o3.loc[~o3.index.duplicated(keep='first')]

#so2 = wb["SO2"][["DateTime", "SO2"]]
so2 = pd.read_csv("SO2.csv",index_col = 0,delimiter=",")
so2.index = pd.to_datetime(so2.pop("DateTime"))
so2 = so2.loc[~so2.index.duplicated(keep='first')]

# # Select the Max_Min Time(Automatically identify start and end date where all data are present)
# start_time = max([np.amin(time_and_temperature.index), np.amin(time_and_Rh.index),\
#     np.amin(hourly_data.index),np.amin(no2.index),np.amin(co.index),np.amin(o3.index), np.min(so2.index)])

# end_time = min([np.amax(time_and_temperature.index), np.amax(time_and_Rh.index),\
#     np.amax(hourly_data.index),np.amax(no2.index),np.amax(co.index),np.amax(o3.index), np.max(so2.index)])


# Combine inputs
result = pd.concat([hourly_data, no2, co, o3, so2, time_and_temperature], axis=1)
result = result.rename(columns = {'NO2 (ug/m3)': 'NO2', 'Hourly CO (mg/m3)': \
    'CO', 'Hourly O3 (ug/m3)': 'O3', 'SO2 (ug/m3)': 'SO2'}, inplace = False)
result = result[start_time:end_time].fillna(0)

# Convert to ppb
result["NO2"] = result["NO2"].div(12.187*46.0055/(273.15 + result["Air Temp (oC)"]))
result["CO"] = result["CO"].div(12.187*28.01/(273.15 + result["Air Temp (oC)"])/1000)
result["O3"] = result["O3"].div(12.187*48/(273.15 + result["Air Temp (oC)"]))
result["SO2"] = result["SO2"].div(12.187*64.07/(273.15 + result["Air Temp (oC)"]))

# Choose which input we are interested in
result = result[species_list]

# Particle processing
particle = pb
#particle.index = pd.to_datetime(particle.pop("Date").astype(str) + " " + particle.pop("Start Time").astype(str), errors='ignore')
particle.index = pd.to_datetime(particle.pop("DateTime"))
particle = particle.loc[~particle.index.duplicated(keep='first')]
#particle = particle.iloc[:, 0:]

# Set the bin spacing
if toggle_log == True:
    # logarithmic method
	rad_bounds = 10.0**(np.linspace(np.log10(lower_part_size), 
						np.log10(upper_part_size), num=(number_size_bins+1)))
	rwid = (rad_bounds[1::]-rad_bounds[0:-1]) # width of size bins (um)
	x_output = rad_bounds[0:-1]+rwid/2.0 # particle radius (um)
        
    #interval = np.logspace(math.log10(lower_part_size), math.log10(upper_part_size), num=(number_size_bins + 1))
else:
    #interval = np.linspace(particle.columns[0], particle.columns[-1], num=number_size_bins)
    rad_bounds = np.linspace(lower_part_size, upper_part_size, (number_size_bins+1))
	# width of size bins (um)
    rwid = np.array((rad_bounds[1]-rad_bounds[0])).reshape(1)
    x_output = rad_bounds[0:-1]+rwid/2.0 # particle radius (um)
# ---------------------------------------
# enhance upper radius bound (um) to reduce possibility of particles growing beyond 
# this (reversed in saving.py)
upper_bin_rad_amp = 1.0e6
rad_bounds[-1] = rad_bounds[-1]*upper_bin_rad_amp
    
# Create a deep copy
conv_particle = particle.copy()
conv_particle = conv_particle.iloc[:,:0]

'''
# if number concentration (#/cc (air)) explicitly stated in inputs
pmode = 1
if (pmode == 1):
    Nperbin = np.array((pconc))
    Nperbin = Nperbin.reshape(-1, 1) # ensure correct shape
	# volume of single particles per size bin (um3) - use with lognormal method
    Varr = ((4.0*np.pi)/3.0)*(x_output**3.0)
	# volume bounds (um3)
    V_bounds = ((4.0*np.pi)/3.0)*(rad_bounds**3.0)
'''
'''
# Classify according to size
for bounds in interval:
    temp = np.zeros(len(particle.index))
    for size in particle.columns:
        if (abs(bounds - float(size)) > 1e-8):
            temp = np.add(temp, (particle.pop(size).astype(float)))
            break
    conv_particle["%f" %(bounds)] = temp
'''
# Classify according to size
i=0
while i < (len(rad_bounds)-1):
    temp = np.zeros(len(particle.index))
    lower_bound = rad_bounds[i]
    upper_bound = rad_bounds[i+1]
    #print(lower_bound, upper_bound)
    for size in particle.columns:
        #size = particle.columns[0]
#if (abs(bounds - float(size)) > 1e-8):
        if (float(size) >= lower_bound and float(size) < upper_bound):
            temp = np.add(temp, (particle.pop(size).astype(float)))
            break
    conv_particle["%f" %(x_output[i])] = temp
    i=i+1
    #print (i)


# Set index to seconds
conv_particle = conv_particle.sort_index().loc[start_time:end_time]
date_time_part = np.array(conv_particle.index)
time_seconds_part = (np.divide((np.subtract(date_time_part, np.datetime64(start_time))), 1E9)).astype(int)
conv_particle = conv_particle.reset_index(drop=True)
conv_particle = conv_particle.set_index(time_seconds_part)


### Writing the .txt file starts here #####################################################################################################
###########################################################################################################################################

f = open(file_name, "w+")

# Total model time
f.write("total_model_time = %d\n" %((end_time - start_time).total_seconds()))

# Write for Temperature

# Convert to seconds
time_and_temperature = time_and_temperature[start_time:end_time]
date_time_temperature = np.array(time_and_temperature.index)
time_seconds_temperature = np.divide((np.subtract(date_time_temperature, np.datetime64(start_time))), 1E9)
air_temp = np.add(np.array(time_and_temperature["Air Temp (oC)"]), 273.15)

temperature_array = np.array([time_seconds_temperature,air_temp]) 

# Write
f.write("temperature = ")
for i in range(len(temperature_array[1]) - 1):
    f.write("%f, " %(temperature_array[1][i]))
f.write("%f\n" %(temperature_array[1][len(temperature_array[1]) - 1]))

f.write("tempt = ")
for i in range(len(temperature_array[0]) - 1):
    f.write("%d, " %(temperature_array[0][i]))
f.write("%d\n" %(temperature_array[0][len(temperature_array[0]) - 1]))

# Write for Rh

# Convert to seconds
time_and_Rh = time_and_Rh[start_time:end_time]
date_time_Rh = np.array(time_and_Rh.index)
time_seconds_Rh = np.divide((np.subtract(date_time_Rh, np.datetime64(start_time))), 1E9)
air_Rh = np.divide(np.array(time_and_Rh["Relative Humidity (%)"]), 100)

Rh_array = np.array([time_seconds_Rh,air_Rh])

# Write
f.write("rh = ")
for i in range(len(Rh_array[1]) - 1):
    f.write("%f, " %(Rh_array[1][i]))
f.write("%f\n" %(Rh_array[1][len(Rh_array[1]) - 1]))

f.write("rht = ")
for i in range(len(Rh_array[0]) - 1):
    f.write("%d, " %(Rh_array[0][i]))
f.write("%d\n" %(Rh_array[0][len(Rh_array[0]) - 1]))

# Write for C0 and Comp0

c0 = result.iloc[0, :len(result.columns)]

f.write("C0 = ")
for i in range(len(c0) -1):
    f.write("%f, " %(c0[i]))
f.write("%f\n" %(c0[len(c0) - 1]))

f.write("Comp0 = ")
for i in range(len(c0.index) -1):
    f.write("%s, " %(mcm[c0.index[i]]))
f.write("%s\n" %(mcm[c0.index[len(c0.index) - 1]]))

# Processing for Ct, Compt, injectt

# Set new index to seconds
ct = result
date_time_ct = np.array(ct.index)
time_seconds_ct = (np.divide((np.subtract(date_time_ct, np.datetime64(start_time))), 1E9)).astype(int)
ct = ct.reset_index(drop=True)
ct= ct.set_index(time_seconds_ct)

# Clear rows with empty values
ct = ct.loc[(ct != 0).any(axis=1)]
ct = ct.iloc[1: , :]

# Write Ct
f.write("Ct = ")
last_col = ct.columns[-1]
last_row = ct.index[-1]
for col in ct.columns:
    for row in ct.index:
        if (col == last_col and row == last_row):
            f.write("%f\n" %(ct[col][row]))
            break
        if row == last_row:
            f.write("%f; " %(ct[col][row]))
            break
        f.write("%f, " %(ct[col][row]))

# Write Compt
f.write("Compt = ")
last = ct.columns[-1]
for items in ct.columns:
    if items == last:
        f.write("%s\n" %(mcm[items]))
        break
    f.write("%s, " %(mcm[items]))

# Write injectt
f.write("injectt = ")
for i in range(len(ct.index) - 1):
    f.write("%d, " %(ct.index[i]))
f.write("%d\n" %(ct.index[len(ct.index) - 1]))

# Write number_size_bins
if toggle_particle == False:
    f.write("number_size_bins = 0\n")
else:
    f.write("number_size_bins = %d\n" %(number_size_bins))

# Write upper and lower bin size bounds
if toggle_particle == False:
    pass
else:
    f.write("lower_part_size = %f\n" %(lower_part_size))
    f.write("upper_part_size = %f\n" %(upper_part_size))

# Write space_mode
if toggle_log == False:
    f.write("space_mode = lin\n")
else:
    f.write("space_mode = log\n")

# Write pcont
f.write("pcont = ")
for i in range(len(conv_particle.index) - 1):
    f.write("%d; " %(0))
f.write("%d\n" %(0))

# Write pconc and pconct
f.write("pconc = ")
last_col = conv_particle.columns[-1]
last_row = conv_particle.index[-1]
for row in conv_particle.index:
    for col in conv_particle.columns:
        if (col == last_col and row == last_row):
            f.write("%f\n" %(conv_particle[col][row]))
            break
        if col == last_col:
            f.write("%f; " %(conv_particle[col][row]))
            break
        f.write("%f, " %(conv_particle[col][row]))

f.write("pconct = ")
for i in range(len(conv_particle.index) - 1):
    f.write("%d; " %(conv_particle.index[i]))
f.write("%d\n" %(conv_particle.index[len(conv_particle.index) - 1]))

f.close()
