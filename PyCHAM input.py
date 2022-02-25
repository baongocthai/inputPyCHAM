import pandas as pd
import numpy as np
import input_processing as ip
import input_main

### Manual Input #############################################################################################
##############################################################################################################

# Set working directory
import os
folder_name = r'C:\Users\janet\Documents\nBox\SOAFP\6. PyCHAM simulation\Input\Input'
os.chdir(folder_name)

# Input time manually: Format: 'year-month-day hour:minute:second')
start_time = pd.to_datetime('2020-11-20 15:00:00')
end_time = pd.to_datetime('2020-11-20 15:59:00')

# Input file_name for model variable input
file_name = "20Nov_1500_1559_20bins.txt"

# Desired VOC data resolution input. For example, input of VOC every 5 min
resolution = 300 #every 5min
index = list(range(0,(end_time.hour-start_time.hour+1)*60*60, resolution))

# Enter number_size_bins
number_size_bins = 20

# Set upper and lower bounds (um)
lower_part_size = 0.001*1000
upper_part_size = 0.380*1000

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
DayOfYear, daytime_start = ip.proc_date_time(start_time)

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
# DCBENE does not have any reaction in MCM inventory -> to be removed. In total, there are 25-1=24 species from TDGCMS
const_comp = ["BENZAL", "BUTACID", "C5H11CHO", "C6H13CHO", "CYHEXONE", #"DCBENE", 
              "DIETETHER", "DIME35EB", "DMC", "ETHACET",  #10
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
                "CHCL3", #90 can be estimated from 2011-2012
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

# Input processing and writing

input_main.main(folder_name, toggle_particle, species_list, inorganic_list, start_time, end_time, toggle_log, 
        lower_part_size, upper_part_size, number_size_bins, file_name, update_step, umansysprop_update, 
        p_init, lat, lon, DayOfYear, daytime_start, light_time, light_status, coag_on, index, const_comp,
        tracked_comp)
