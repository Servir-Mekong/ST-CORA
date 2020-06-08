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

import pickle as pkl
import numpy as np
#import scipy.io
#import ellipsoid as ET
import ST_MultiCORA_v1_0 as stcora
import Displacement_extraction_Multi as disp
import Magnitude_extract_Multi as mag
from datetime import  timedelta
import time
from osgeo import gdal
import os
import warnings
warnings.filterwarnings("ignore")
#from multiprocessing import Pool
import skimage as skm 


############################################################################
#%%                              FUNCTIONS
#############################################################################
def spatial_ref(P,Dates):
    
    [Xmin_lbm,Xmax_lbm,Ymin_lbm,Ymax_lbm] = P['boundary'] #Boundary THAILAND  
    lon = np.int((Xmax_lbm - Xmin_lbm) / P['pixel_value'])
    lat = np.int((Ymax_lbm - Ymin_lbm) / P['pixel_value'])
    
    try:
        Mask = gdal.Open(P['Mask_path']).ReadAsArray()
    except:
        Mask = np.empty([lat,lon]) 
        Mask.fill(1)
    
    DATA_OBS = np.empty([lon,len(Dates)])
    DATA_SAT = np.empty([lon,len(Dates)])
    t=0
    for Date in Dates:
        DateH= Date + timedelta(hours=14)         
        GPM_name = 'rain_{0}{1:02d}{2:02d}{3:02d}0000.asc'.format(Date.year,Date.month,Date.day,Date.hour)
        HAII_name = 'rain_{0}{1:02d}{2:02d}{3:02d}.asc'.format(DateH.year,DateH.month,DateH.day,DateH.hour) 
        DATA_SAT[:,t] = np.max(gdal.Open(os.path.join(P['GPM_dir'],GPM_name)).ReadAsArray() * Mask,0 )
        DATA_OBS[:,t] = np.max(gdal.Open(os.path.join(P['HAII_dir'],HAII_name)).ReadAsArray()[:-1,:-1] * Mask,0) # remove last  row and column
        t=t+1  
    SPACE_SEARCH = np.maximum(DATA_SAT.T,DATA_OBS.T )
    #SPACE_SEARCH[SPACE_SEARCH<10]=0
    SPACE_SEARCH = (SPACE_SEARCH > P['Inten_Seach']).astype(int)
    SPACE_label  = skm.measure.label(SPACE_SEARCH, background=0)
    SPACE_clean = skm.morphology.remove_small_objects(SPACE_label, P['Min_seachsize'])
    SPACE_regions = skm.measure.regionprops(SPACE_clean)
    SPACE_regions = np.array([list(props.bbox) for props in SPACE_regions] )
    
    Temporal_filter = np.array([SPACE_regions[:,0]-5,SPACE_regions[:,2]+5]).T
    Temporal_filter[Temporal_filter<0]=0
    #plt.figure()
    #plt.imshow(CC2D_clean)
    return Temporal_filter

def sim_time(start_time):
    elapsed_time = time.time() - start_time
    print(time.strftime("%H:%M:%S", time.gmtime(elapsed_time)))
def dist(CenSat,CenObs):
    D = np.sqrt( (CenSat[0]-CenObs[0])**2 + (CenSat[1]-CenObs[1])**2 + (CenSat[2]-CenObs[2])**2)
    return D

def bb_intersection_over_union(boxA, boxB):
	# determine the (x, y)-coordinates of the intersection rectangle
    xA = max(boxA[0], boxB[0])
    yA = max(boxA[1], boxB[1])
    zA = max(boxA[2], boxB[2])
    xB = min(boxA[3], boxB[3])
    yB = min(boxA[4], boxB[4]) 
    zB = min(boxA[5], boxB[5])
        
	# compute the area of intersection rectangle
    interArea = max(0, xB - xA + 1) * max(0, yB - yA + 1) * max(0, zB - zA + 1) 
	# compute the area of both the prediction and ground-truth
	# rectangles
    boxAArea = (boxA[3] - boxA[0] + 1) * (boxA[4] - boxA[1] + 1) * (boxA[5] - boxA[2] + 1)
    boxBArea = (boxB[3] - boxB[0] + 1) * (boxB[4] - boxB[1] + 1) * (boxA[5] - boxA[2] + 1) 
	# compute the intersection over union by taking the intersection
	# area and dividing it by the sum of prediction + ground-truth
	# areas - the interesection area
    iou = interArea / float(boxAArea + boxBArea - interArea) 
	# return the intersection over union value
    return iou
###############################################################################
#%%               #############  main process ##########
###############################################################################

#GPM_cluster ={}
#Total_time = time.time() 
#for Group in Date_group:
def main(P,Group,st):      
    
    [Xmin_lbm,Xmax_lbm,Ymin_lbm,Ymax_lbm] = P['boundary'] #Boundary THAILAND   

    lat = np.int((Ymax_lbm - Ymin_lbm) / P['pixel_value'])
    lon = np.int((Xmax_lbm - Xmin_lbm) / P['pixel_value'])

    ###### Start initial values #############
    try:
        Mask = gdal.Open(P['Mask_path']).ReadAsArray()
    except:
        Mask = np.empty([lat,lon]) 
        Mask.fill(1)

    # SRE Coordinates
    LAT= np.linspace(Ymin_lbm,Ymax_lbm, num = lat)  # Array latitude
    LON= np.linspace(Xmin_lbm,Xmax_lbm, num = lon)  # Array longitude
    LAT = LAT[::-1]
    
    
#    Group=Date_group[4]
    DATA_OBS = np.empty([lat,lon,len(Group)])
    DATA_SAT = np.empty([lat,lon,len(Group)])
    t=0
    for Date in Group:
        DateH= Date + timedelta(hours=14)         
        GPM_name = 'rain_{0}{1:02d}{2:02d}{3:02d}0000.asc'.format(Date.year,Date.month,Date.day,Date.hour)
        HAII_name = 'rain_{0}{1:02d}{2:02d}{3:02d}.asc'.format(DateH.year,DateH.month,DateH.day,DateH.hour) 
        DATA_SAT[:,:,t] = gdal.Open(os.path.join(P['GPM_dir'],GPM_name)).ReadAsArray()  * Mask  
        DATA_OBS[:,:,t] = gdal.Open(os.path.join(P['HAII_dir'],HAII_name)).ReadAsArray()[:-1,:-1]  * Mask # remove last  row and column
        t=t+1         
    ###############################################################################
    ############    1.              ST-CORA
    ###############################################################################
    print('Processing ST-CORA...{}'.format(st)) 
#    start_time = time.time()
    
    #ST_CORA
    Objects_OBS =  stcora.ST_CORA(DATA_OBS,P['T'],P['T2'],P['Prob'],LAT,LON,P['Minsize'],P['Psize'])
    Objects_SAT =  stcora.ST_CORA(DATA_SAT,P['T'],P['T2'],P['Prob'],LAT,LON,P['Minsize'],P['Psize'])    
    
    DATA_SAT,DATA_OBS = [None,None] # flush
    # math events
    WCentroides_OBS = np.array([np.array(Objects_OBS.get(item)['centroid'])
     for item in range(len(Objects_OBS))])
    WCentroides_SAT = np.array([np.array(Objects_SAT.get(item)['centroid']) 
    for item in range(len(Objects_SAT))])
    
    #Time in meters    
    #SL = 154285714.2
    #WCentroides_OBS[:,2] = WCentroides_OBS[:,2] *SL
    #WCentroides_SAT[:,2] = WCentroides_SAT[:,2] *SL   
    
    # Matching objects
    bbox_OBS = [list(Objects_OBS.get(item)['Region_box'])
     for item in range(len(Objects_OBS))]
    bbox_SAT = [list(Objects_SAT.get(item)['Region_box']) 
    for item in range(len(Objects_SAT))]
        
    #  bb_intersection_over_union
    matrix_R  = np.zeros([len(WCentroides_OBS),len(WCentroides_SAT)])   
    for cen_obs in range(len(WCentroides_OBS)):
        for cen_sat in range(len(WCentroides_SAT)):    
            matrix_R[cen_obs,cen_sat] = bb_intersection_over_union(bbox_OBS[cen_obs], bbox_SAT[cen_sat])        
    
    Min_cor = np.argmax(matrix_R,1)
    Object_matched = Min_cor[np.max(matrix_R,1) >= P['Inter_seach']]   # number ob satObject per satObj
    
    if len(Object_matched) == 0:
        print('not matched found')   


