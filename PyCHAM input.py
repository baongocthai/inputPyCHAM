import datetime
import pandas as pd
import numpy as np

### Manual Input #############################################################################################
##############################################################################################################

# Set working directory
import os
folder_name = r'C:\Users\janet\Documents\nBox\SOAFP\6. PyCHAM simulation\Input\Input'
os.chdir(folder_name)

# Input time manually: Format: 'year-month-day hour:minute:second')
start_time = pd.to_datetime('2020-11-05 09:00:00')
end_time = pd.to_datetime('2020-11-05 17:59:00')

# Input file_name for model variable input
file_name = "05Nov_0900_1759_20bins_89VOCs.txt"

# Desired VOC data resolution input. For example, input of VOC every 5 min
resolution = 300 #every 5min
index = list(range(0,(end_time.hour-start_time.hour+1)*60*60, resolution))

# Enter number_size_bins
number_size_bins = 164

# Set upper and lower bounds (um)
lower_part_size = 0.002*1000
upper_part_size = 0.760*1000

### Other parameters #########################################################################################
##############################################################################################################

# Input time interval for each time step (s)
update_step = 60

# Toggle unmansysprop update: 1=True, 0=False
umansysprop_update = 1

# Constant pressure
p_init = 100000

# Co-ordinates of measurement location
lat = 1.3 
lon = -103.8 

### Day of Year & day time start
# Convert local to (UTC+8) to GMT (UTC+0)
actual_datetime = start_time - datetime.timedelta(hours=8, minutes=00, seconds=00)
DayOfYear = int(actual_datetime.strftime('%j')) #no. of day in the year
daytime_start = actual_datetime.hour*60*60 + actual_datetime.minute*60 #in seconds

# Light status & light time
light_time = 0
light_status = 1

# Toggle for coagulation: 1=True, 0=False
coag_on = 1

### Input name of species to track for rate of loss/production
# Must include coresponding MCM name instead of normal name; i.e., Acrolein as ACR
tracked_comp = "OH, O3, HONO"

# species with constant concentration
# these 25 species are not measured in 2011-2012, hence not estimated for
# these 25 species are measured in TD-GC/MS
# Therefore, the median concentration of these 25 species for 9 days in Sep.-Oct. 2021 are used and set as constant
const_comp = ["BENZAL", "BUTACID", "C5H11CHO", "C6H13CHO", "CYHEXONE", "DCBENE", "DIETETHER", "DIME35EB", "DMC", "ETHACET",  #10
            "HEX1ENE", "HEX3ONE", "IPROPOL", "LIMONENE", "M3F", "MACR", "METHACET", "MIBK", "MPRK", "NBUTACET", #20
            "NC11H24", "NC12H26", "PENTACID", "SBUTACET", "THEX2ENE"] #25

# Insert the desired the species
# Overlapping species between PTRMS & TDGCMS, and required to use TDGCMS data has been scaled by a constant ratio
# The ratio = median concentration (Sep.-Oct.2021, PTRMS) / median concentration (Sep.-Oct.2021, TDGCMS)
# These 13 species are: total pentene, total trimethylbenzene, total xylene, total pinene, toluene, benzene, total butene, 
# dichloromethane, isoprene, MVK, butanone (MEK),  nonanal, decanal
# median dichloromethane (Sep.-Oct.2021, PTRMS) = 0 ppb => replace all values with median concentration (Sep.-Oct.2021, TDGCMS) 
# nonanal and decanal arenot available in MCM inventory
species_list = ["Chloromethane", "Dichloro.methane", "Ethyl.nitrate", "Pinic.acid", "Pinonic.acid",#5
                "Dimethyl.Sulfide", 
                "Isoprene", 
                "Acetone", 
                "Toluene", 
                "Styrene", #10
                "Acetic.anhydride", "Acetaldehyde", "Acetic.acid",
                "Acrolein", "Acrylic.acid", "Benzene", "Butanone", "Catechol", "Cresol", "Ethanol", #20
                "Formaldehyde", "Formic.acid",
                "Hexanal", "Methanesulfonic.acid", "Methyl.nitrate",#25
                "MVK", "p.Benzoquinone", 
                "Pentanal", "Peroxyacetyl.Nitrate", "Phenol",#30
                "Propene","Propyl.acetate","Pyruvic.acid","Acetol", #34
                # 6 species with isomers
                "OXYL", "MXYL", "PXYL", # Xylene
                "TM123B", "TM124B","TM135B", #TMB
                "APINENE", "BPINENE", #Pinene
                "BUT1ENE", "CBUT2ENE", "MEPROPENE", "TBUT2ENE", #Butene
                "PENT1ENE", "ME2BUT1ENE", "ME2BUT2ENE", "ME3BUT1ENE", "CPENT2ENE","TPENT2ENE", #Pentene
                "Butanol", #53
                # Other 2011-2012 species
                "CDICLETH", "C4H6", "M22C4", "M23C4", "NC4H9NO3", "M2HEX", "M2PE", #60
                "PEANO3", "M3HEX", "M3PE", "TCE", "TRICLETH", #65
                "CH3BR", "CH3CCL3", "CH4", "CHEX", "C2H6", #70
                "C2H4", "EBENZ", "C2H2","IC4H10", "IC5H12", #75
                "IC3H7NO3", "NC4H10", "NC10H22","NC7H16", "PBENZ",#80
                "NC6H14", "NC9H20", "NC8H18", "NC5H12", "NC3H7NO3", #85
                "C3H8","OETHTOL", "METHTOL", "PETHTOL", #89
                # "CHCL3", #90 can be estimated from 2011-2012
                "NO2", "O3", "SO2", "CO"] # Inorganic species
'''
species_list = ['APINENE','Toluene', "C2H4", "Propene", "Formaldehyde", "OXYL", "Acetaldehyde", "C4H6", "Isoprene",
                "NO2", "O3", "SO2", "CO"]
'''

inorganic_list = ["NO2", "O3", "SO2", "CO"]

## Toggle particle consideration
toggle_particle = True

# Toggle logarithmic scale consideration
toggle_log = True

#######################################################################################################################
#######################################################################################################################

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
    file_particle = r"SMPS.csv"
    pb = pd.read_csv(file_particle, delimiter = ",", dtype = object)

# Processing for time and temperature
#Meteorology = pd.read_csv("Meteorology.csv",index_col = 0,delimiter=",")
wb_temp_rh.index = pd.to_datetime(wb_temp_rh.pop("Date&Time"), format = '%d/%m/%Y %H:%M')
time_and_temperature = wb_temp_rh[["Temperature"]]
time_and_temperature = time_and_temperature.dropna()

# Processing for time and Rh
time_and_Rh = wb_temp_rh[["Humidity"]]
time_and_Rh = time_and_Rh.dropna()

# Processing input names
input_names = pd.read_csv("Species as inputs_BN2.csv",index_col = 0,delimiter=",")
input_names = input_names[input_names["Required unit"] == "ppb"]
mcm = dict(zip(input_names["Inputs"], input_names["Name in MCM"]))
mcm_data = {x: mcm[x] for x in mcm if x not in inorganic_list}

# Double check number of species
mcm_keys = [x for x in enumerate(mcm) if x in species_list]
if len(mcm_keys) == 0: print("All species are accounted for")

# Processing input for every minute
#VOC_data = pd.read_csv("Minute data_BN2.csv",delimiter=",")
# VOC_data = pd.read_csv("Minute data_BN2.csv",delimiter=",")
wb_species.index = pd.to_datetime(wb_species.pop("DateTime"),format = '%d/%m/%Y %H:%M')
# Replace all BD values to 0
wb_species.replace({"BD":0},inplace=True)
wb_species = wb_species.astype(float)

# Data check
print("All numbers are float" if sum([1 if x=="float64" else 0 for x in wb_species.dtypes]) == wb_species.dtypes.count() else "Error: not all numbers are float")
wb_species = wb_species[mcm_data.keys()]

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

# Hourly temperature
hourly_temperature = pd.read_csv("Hourly Meteorology.csv",index_col = 0,delimiter=",")
hourly_temperature = hourly_temperature[["Date&Time", "Air Temp (oC)"]]
hourly_temperature.index = pd.to_datetime(hourly_temperature.pop('Date&Time'),format = '%d/%m/%Y %H:%M')


# Combine inputs
result = wb_species.join([wb_no2, wb_co, wb_o3, wb_so2, hourly_temperature], how="outer")
result = result.rename(columns = {'NO2 (ug/m3)': 'NO2', 'Hourly CO (mg/m3)': \
    'CO', 'Hourly O3 (ug/m3)': 'O3', 'SO2 (ug/m3)': 'SO2'}, inplace = False)
result = result[start_time:end_time]


# Convert to ppb
result["NO2"] = result["NO2"].div(12.187*46.0055/(273.15 + result["Air Temp (oC)"]))
result["CO"] = result["CO"].div(12.187*28.01/(273.15 + result["Air Temp (oC)"])/1000)
result["O3"] = result["O3"].div(12.187*48/(273.15 + result["Air Temp (oC)"]))
result["SO2"] = result["SO2"].div(12.187*64.07/(273.15 + result["Air Temp (oC)"]))


# Choose which input we are interested in
result = result[species_list]

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

### Writing the .txt file starts here #####################################################################################################
###########################################################################################################################################

f = open(file_name, "w+")
# resutl file name
f.write("res_file_name = %s\n" %(file_name[:len(file_name)-4]))

# Total model time
f.write("total_model_time = %d\n" %((end_time - start_time).total_seconds()+60))

# Time step
f.write("update_step = %d\n" %(update_step))
f.write("recording_time_step = %d\n" %(update_step))

# umansysprop_update
f.write("umansysprop_update = %d\n" %(umansysprop_update))

# Constant pressure
f.write("p_init = %f\n" %(p_init))

# Co-ordiante of SG
f.write("lat = %f\n" %(lat))
f.write("lon = %f\n" %(lon))

# Date time of year
f.write("DayOfYear = %d\n" %(DayOfYear))
f.write("daytime_start = %d\n" %(daytime_start))

# Natural light status
f.write("light_time = %d\n" %(light_time))
f.write("light_status = %d\n" %(light_status))

# Toggle for coagulation
f.write("coag_on = %d\n" %(coag_on))

# Track rate of loss/production of species
f.write("tracked_comp = %s\n" %(tracked_comp))

# species with constant concentration
f.write("const_comp = ")
for i in range(len(const_comp) -1):
    f.write("%s, " %(const_comp[i]))
f.write("%s\n" %(const_comp[len(const_comp) - 1]))

# Turn off wall-partitioning
f.write("wall_on = 0\n")

# Size structure = moving-center
f.write("size_structure = 0\n")

# Input for nucleation
f.write("nucv1 = 30000\n")
f.write("nucv2 = -4\n")
f.write("nucv3 = 100\n")

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

# Write for C0 and Comp0
c0 = result.iloc[0, :len(result.columns)]
# Add constant concentration in C0
c0_additional_name = pd.Series(["BENZAL", "BUTACID", "C5H11CHO", "C6H13CHO", "CYHEXONE", "DCBENE", "DIETETHER", "DIME35EB", "DMC", "ETHACET",  #10
            "HEX1ENE", "HEX3ONE", "IPROPOL", "LIMONENE", "M3F", "MACR", "METHACET", "MIBK", "MPRK", "NBUTACET", #20
            "NC11H24", "NC12H26", "PENTACID", "SBUTACET", "THEX2ENE", 'OH'])
c0_additional_conc = pd.Series([0.09703243, 0, 0, 0.028099814, 0, 0.016587651, 0, 0, 0.01925449, 0, 0.266989242, 0, 0.081819296, 0.081030143, 
        0.004791663, 0.009176899, 0.071893593, 0.136275387, 0.030715693, 0.071789396, 0.031384498, 0.047690948, 0, 0, 0, 6.5e-5])
c0_combined = pd.concat([c0_additional_name, c0_additional_conc], axis=1)
c0_combined.index = c0_combined[0]
c0_combined = c0_combined[1]
c0 = pd.concat([c0, c0_combined], axis=0)

f.write("C0 = ")
for i in range(len(c0) -1):
    f.write("%f, " %(c0[i]))
f.write("%f\n" %(c0[len(c0) - 1]))

f.write("Comp0 = ")
for i in range(len(species_list) -1):
        f.write("%s, " %(mcm[c0.index[i]]))

# Add in species with constant concentration
for i in const_comp:
    f.write("%s, " %i)
# Add in for OH
f.write("%s\n" %(c0.index[len(c0.index) - 1]))

### Processing for Ct, Compt, injectt
# Set new index to seconds
ct = result
date_time_ct = np.array(ct.index)
time_seconds_ct = (np.divide((np.subtract(date_time_ct, np.datetime64(start_time))), 1E9)).astype(int)
ct = ct.reset_index(drop=True)
ct= ct.set_index(time_seconds_ct)

# Clear rows with empty values
ct = ct.loc[(ct != 0).any(axis=1)]
ct_new = pd.DataFrame([ct.iloc[i] for i in range(len(ct)) if ct.index[i] in index])
ct_new = ct_new.iloc[1: , :]

# Write Ct
f.write("Ct = ")
last_col = ct_new.columns[-1]
last_row = ct_new.index[-1]
for col in ct_new.columns:
    if col not in inorganic_list:
        for row in ct_new.index:
            if (col == last_col and row == last_row):
                f.write("%f\n" %(ct_new[col][row]))
                break
            if row == last_row:
                f.write("%f; " %(ct_new[col][row]))
                break
            f.write("%f, " %(ct_new[col][row]))

# Write Compt
f.write("Compt = ")
last = ct.columns[-1]
for items in ct.columns:
    if col not in inorganic_list:
        if items == last:
            f.write("%s\n" %(mcm[items]))
            break
        f.write("%s, " %(mcm[items]))

# Write injectt
f.write("injectt = ")
for i in range(len(ct_new.index) - 1):
    f.write("%d, " %(ct_new.index[i]))
f.write("%d\n" %(ct_new.index[len(ct_new.index) - 1]))

# Write number_size_bins
if toggle_particle == False:
    f.write("number_size_bins = 0\n")
else:
    f.write("number_size_bins = %d\n" %(number_size_bins))

# Write upper and lower bin size bounds
if toggle_particle == False:
    pass
else:
    f.write("lower_part_size = %f\n" %(lower_part_size/1000))
    f.write("upper_part_size = %f\n" %(upper_part_size/1000))

# Write space_mode
if toggle_log == False:
    f.write("space_mode = lin\n")
else:
    f.write("space_mode = log\n")
    
# Write mean_rad
if toggle_particle == True:
    f.write("mean_rad = ")
    for i in range(len(conv_particle.index) - 1):
        f.write("%f; " %((upper_part_size + lower_part_size)/2))
    f.write("%f\n" %((upper_part_size + lower_part_size)/2))

# Write std
if toggle_particle == True:
    f.write("std = ")
    for i in range(len(conv_particle.index) - 1):
        f.write("%f; " %(1.2))
    f.write("%f\n" %(1.2))
    
# Write pcont
f.write("pcont = ")
for i in range(len(conv_particle.index) - 1):
    f.write("%d; " %(0))
f.write("%d\n" %(0))

# Write pconc and pconct
if toggle_particle == True:
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
    f.write("%d" %(conv_particle.index[len(conv_particle.index) - 1]))


f.close()
