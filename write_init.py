import pandas as pd
import numpy as np
# All initial & constant variables which does not require external inputs
def write_init(file_name, start_time, end_time, update_step, umsp_up, p_init, lat, lon, 
                DayOfYear, daytime_start, light_time, light_status, coag_on):

    f = open(file_name, "w+")

    # result file name
    f.write("res_file_name = %s\n" %(file_name[:len(file_name)-4]))

    # Total model time
    f.write("total_model_time = %d\n" %((end_time - start_time).total_seconds()+60))

    # Time step
    f.write("update_step = %d\n" %(update_step))
    f.write("recording_time_step = %d\n" %(update_step))

    # umansysprop_update
    f.write("umansysprop_update = %d\n" %(umsp_up))

    # Constant pressure
    f.write("p_init = %f\n" %(p_init))

    # Co-ordinate of SG
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

    # Turn off wall-partitioning
    f.write("wall_on = 0\n")

    # Size structure = moving-center
    f.write("size_structure = 0\n")

    # Input for nucleation
    f.write("nucv1 = 30000\n")
    f.write("nucv2 = -4\n")
    f.write("nucv3 = 100\n")

    f.close()

    return()
