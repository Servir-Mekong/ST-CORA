# -*- coding: utf-8 -*-
"""
###############################################################################
                        ST-CORA v2.1 master
###############################################################################
                        
@author: M. Laverde-Barajas
@email: m.laverde@un-ihe.org

@company SERVIR- MEKONG program
        IHE Delft Institute for Water Education
        Delft University of Technology


## Paramaters
T = 0.5 # wet values
T2 = 0.5  # delineation
Minsize = 15 # 64 km2
Prob =  0.25  # kernel segmentation 
Psize = 100   # Min Object size
pixel_value = 0.1  # resolution
method = 4  # magnitude type
Inter_seach = 0.3 # minimum intersection rate betwwen storm objects


StartTime = datetime(2015, 6, 1, 0, 0, 0) #  LOCAL TIME (GMT + 7)
EndTime = datetime(2015, 10, 31, 23, 0, 0) 

[Xmin_lbm,Xmax_lbm,Ymin_lbm,Ymax_lbm] = [100,106,14,19] #Boundary THAILAND

## SPATIAL SEARCH
Inten_Seach = 10
Min_seachsize = 20

"""

#import pickle as pkl
#import ellipsoid as ET
import functions_stcora as cora_bias
from datetime import datetime
import pandas as pd
import time
import warnings
warnings.filterwarnings("ignore")
#from multiprocessing import Pool

###############################################################################
############                    MASTER ST-CORA
###############################################################################

class BICO(object):    
   
    def runner(Dates,P):
        #%   RUNNER 
        if yyy == 2018:
            ranst = 53 #75 in 2015 and 2016 ;  90 119 123 211 238 in 2017 ; 96 110 111 2018check
        else:
            ranst = 0
        Temporal_filter = cora_bias.spatial_ref(P,Dates)       
        Total_time = time.time() 
        for st in range(ranst,len(Temporal_filter)):
            print('Total_event ' + str(len(Temporal_filter)))
            start_time = time.time() 
            Group = Dates[Temporal_filter[st,0]:Temporal_filter[st,1]]
            try:
                GPM_cluster = cora_bias.main(P,Group,st)
            except:
                print('error in biascora no. {0} from {1}'.format(st,yyy))
            cora_bias.sim_time(start_time)    
        print('TOTAL TIME')
        cora_bias.sim_time(Total_time) 
        
        return GPM_cluster
    
def main(yyy):
    P={}
    P['T'] = 0.1 # wet values 0.5
    P['T2'] = 1 # delineation n0.5 
    P['Minsize'] = 15 # 64 km2
    P['Prob'] =  0.5  # kernel segmentation 0.25
    P['Psize'] = 100   # Min Object size
    P['pixel_value'] = 0.1  # resolution
    P['method'] = 4  # magnitude type
    P['Inter_seach'] = 0.25 # minimum intersection rate between storm objects 0.5
    P['HAII_dir']= r'D:\CORA_BIAS\netcdf_temp\Outputs\HAII_asc'
    P['GPM_dir'] = r'D:\CORA_BIAS\netcdf_temp\Outputs\GPM_asc'
    P['Mask_path'] = r'D:\CORA_BIAS\Files\Mask_rec.tif'    
    P['boundary'] = [100,106,14,19] #Boundary THAILAND
    P['Inten_Seach'] = 10 ## SPATIAL SEARCH
    P['Min_seachsize'] = 20 ## SPATIAL SEARCH
    P['Storm_dir'] = r'D:\CORA_BIAS\Corrected_Objects\{}/'.format(yyy)  
#    P['Storm_dir'] = r'D:\CORA_BIAS\Objects_Files\{}/'.format(yyy) 
    P['KNN'] = 0  # KNN function for displacment 0 no 1 yes
    
    # Select initial and end date
    StartTime = datetime(yyy, 6, 1, 0, 0, 0) #  LOCAL TIME (GMT + 7)
    EndTime = datetime(yyy, 10, 31, 23, 0, 0) 
    Dates =  pd.date_range(StartTime,EndTime, freq='h')
    
    BICO.runner(Dates,P)
        
if __name__ == '__main__':
    YEAR = [2018]
    for yyy in YEAR:
        main(yyy)         