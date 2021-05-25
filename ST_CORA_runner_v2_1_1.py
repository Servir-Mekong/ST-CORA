# -*- coding: utf-8 -*-
"""
###############################################################################
                        ST-CORA v2.2.1 master
###############################################################################
                        
@author: M. Laverde-Barajas
@email: m.laverde@un-ihe.org
        mlaverdeb@gmail.com

@company SERVIR- MEKONG program
        IHE Delft Institute for Water Education
        Delft University of Technology

@citation: Laverde-Barajas, M., et al(2019). 
        Spatiotemporal analysis of extreme rainfall 
        events using an object-based approach. 
        In Spatiotemporal Analysis of Extreme Hydrological Events 
        (pp. 95-112). Elsevier.
        ISBN 9780128116890,
        https://doi.org/10.1016/B978-0-12-811689-0.00005-7.
        (https://www.sciencedirect.com/science/article/pii/B9780128116890000057)
        
## Paramaters
T = 0.5 # wet values
T2 = 0.5  # delineation
Minsize = 15 # 64 km2
kernel =  0.25  # 1 = kernel segmentation 4D 0= 3D
Psize = 100   # Min Object size
pixel_value = 0.1  # resolution

StartTime = datetime(2015, 6, 1, 0, 0, 0) #  LOCAL TIME (GMT + 7)
EndTime = datetime(2015, 10, 31, 23, 0, 0) 

MATRIX = 3D rainfall matrix
boundary = [Xmin_lbm,Xmax_lbm,Ymin_lbm,Ymax_lbm] 

"""
import functions_stcora_v2_1_1 as cora
from datetime import datetime, timedelta
import pandas as pd
import time
import warnings
warnings.filterwarnings("ignore")
#from multiprocessing import Pool
import sys, getopt
import os

###############################################################################
#%%###########                    MASTER ST-CORA
###############################################################################
    
def main():
    
    start_time = time.time()
    
    P={}
    P['T'] = 3  # deliniation  mm/h
    P['T2'] = 20 # segmenation mm/h
    P['kernel'] = 1   # Kernel yes =1  
    P['Minsize'] = 10 # Min Object size to be considered as noise 
    P['Psize'] = 100   # pixel extension in km2
    P['pixel_value'] = 0.1  # resolution in degree
    # P['Mask_path'] = r'Path mask' (optional )
    P['boundary'] = [min_lon, maxlon,Minlat,Maxlat] #Boundary domain 
    P['Storm_dir'] = r' path' 
    P['CPU_cores'] = 1
    # Select initial and end date   
    StartTime ='2020-11-27'
    EndTime = '2021-01-01'
    
    Dates =  pd.date_range(StartTime,EndTime, freq='h')
    
    # RUNNER ST-CORA
    cora.runner(Dates,MATRIX,P)
        
if __name__ == '__main__':
    main()         