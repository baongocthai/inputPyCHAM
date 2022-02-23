import pandas as pd
import numpy as np

def write_particle(file_name, toggle_particle, conv_particle, number_size_bins, lower_part_size, upper_part_size, toggle_log):
    
    f = open(file_name, "a")

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
    if toggle_particle == True and len(conv_particle) > 1:
        f.write("mean_rad = ")
        for i in range(len(conv_particle.index) - 1):
            f.write("%f; " %((upper_part_size + lower_part_size)/2))
        f.write("%f\n" %((upper_part_size + lower_part_size)/2))

    # Write std
    if toggle_particle == True and len(conv_particle) > 1:
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
