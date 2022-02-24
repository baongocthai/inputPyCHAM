'''module to write gas concentration inputs'''

import pandas as pd
import numpy as np

def write_conc(file_name, result, start_time, inorganic_list, mcm, index, species_list, const_comp,
                tracked_comp):

    f = open(file_name, "a")

    # Track rate of loss/production of species
    f.write("tracked_comp = %s\n" %(tracked_comp))
    
    # species with constant concentration
    f.write("const_comp = ")
    for i in range(len(inorganic_list)):
        if inorganic_list[i] != "O3":
            f.write("%s, " %inorganic_list[i])
    
    for i in range(len(const_comp) -1):
        f.write("%s, " %(const_comp[i]))
    f.write("%s\n" %(const_comp[len(const_comp) - 1]))

    # Write for C0 and Comp0
    c0 = result.iloc[0, :len(result.columns)]
    # Add constant concentration in C0
    c0_additional_name = pd.Series(["BENZAL", "BUTACID", "C5H11CHO", "C6H13CHO", "CYHEXONE", #"DCBENE", 
                                    "DIETETHER", "DIME35EB", "DMC", "ETHACET",  #10
                                    "HEX1ENE", "HEX3ONE", "IPROPOL", "LIMONENE", "M3F", "MACR", "METHACET", "MIBK", "MPRK", "NBUTACET", #20
                                    "NC11H24", "NC12H26", "PENTACID", "SBUTACET", "THEX2ENE", 'OH'])
    c0_additional_conc = pd.Series([0.09703243, 0, 0, 0.028099814, 0, #0.016587651, 
                                    0, 0, 0.01925449, 0, #10
                                    0.266989242, 0, 0.081819296, 0.081030143, 0.004791663, 0.009176899, 0.071893593, 0.136275387, 0.030715693, 0.071789396, #20
                                    0.031384498, 0.047690948, 0, 0, 0, 6.5e-5])
    c0_combined = pd.concat([c0_additional_name, c0_additional_conc], axis=1)
    c0_combined.index = c0_combined[0]
    c0_combined = c0_combined[1]
    c0 = pd.concat([c0, c0_combined], axis=0)

    f.write("C0 = ")
    for i in range(len(c0) -1):
        f.write("%f, " %(c0[i]))
    f.write("%f\n" %(c0[len(c0) - 1]))

    f.write("Comp0 = ")
    for i in range(len(species_list)):
            f.write("%s, " %(mcm[c0.index[i]]))

    # Add in species with constant concentration
    for species in const_comp:
        f.write("%s, " %species)
    f.write("%s\n" %(c0.index[len(c0.index) - 1]))

    ### Processing for Ct, Compt, injectt
    # Set new index to seconds
    ct = result
    date_time_ct = np.array(ct.index)
    time_seconds_ct = (np.divide((np.subtract(date_time_ct, np.datetime64(start_time))), 1E9)).astype(int)
    ct = ct.reset_index(drop=True)
    ct= ct.set_index(time_seconds_ct)

    # Clear rows with empty values
    ct = ct.loc[(ct != np.NaN).any(axis=1)]
    ct_new = pd.DataFrame([ct.iloc[i] for i in range(len(ct)) if ct.index[i] in index])
    ct_new = ct_new.iloc[1: , :]

    # Write Ct
    f.write("Ct = ")
    for col in ct_new.columns:
        if col not in inorganic_list:
            last_col = col
    last_row = ct_new.index[-1]
    for col in ct_new.columns:
        if col not in inorganic_list:
            for row in ct_new.index:
                f.write("%f, "%(ct_new[col][row]))
                if (col == last_col and row == last_row):
                    f.write("%f\n" %(ct_new[col][row]))
                    break        
                if row == last_row:
                    f.write("%f; " %(ct_new[col][row]))
                    break

    # Write Compt
    f.write("Compt = ")
    for i in range(len(ct_new.columns)-2):
        if ct_new.columns[i] not in inorganic_list:
            if ct_new.columns[i+1] not in inorganic_list: 
                f.write("%s, " %(mcm[ct_new.columns[i]]))
            else:
                f.write("%s\n" %(mcm[ct_new.columns[i]]))

    # Write injectt
    f.write("injectt = ")
    for i in range(len(ct_new.index) - 1):
        f.write("%d, " %(ct_new.index[i]))
    f.write("%d\n" %(ct_new.index[len(ct_new.index) - 1]))

    f.close()

