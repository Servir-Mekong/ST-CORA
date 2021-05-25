# -*- coding: utf-8 -*-
"""
###############################################################################
                        ST-CORA  v2.2.1 functions
###############################################################################
                        
@author: M. Laverde-Barajas
@email: m.laverde@un-ihe.org

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

boundary = [Xmin_lbm,Xmax_lbm,Ymin_lbm,Ymax_lbm] 

"""

#import pickle as pkl
import numpy as np
import ST_MultiCORA_v2_1_1 as stcora 
import time
import os
import warnings
warnings.filterwarnings("ignore")
import skimage as skm 
from joblib import Parallel, delayed
import multiprocessing
from cc3d import connected_components

############################################################################
#%%                              FUNCTIONS
#############################################################################
                       
############ Clean maps (2D analysis)
def clean_SAT(DATA_SAT,T,Minsize):

    [lat,lon,t]= DATA_SAT.shape  
    DATA_SATClean = np.empty([lat,lon,t])
    for Step in range(t):        
    #    Step=113
        GPM_step = DATA_SAT[:,:,Step]        
        #  remove small dark spots
        SEE = GPM_step > T       
        SEE_clean = skm.morphology.closing(SEE, skm.morphology.square(3)) 
        CC2D  = skm.measure.label(SEE_clean, background=0)          
        # np.count_nonzero((CC2D == [4]))        
        CC2D_clean = skm.morphology.remove_small_objects(CC2D, Minsize)    
 
        Object = np.empty([lat,lon])
        labels = np.unique(CC2D_clean)
        for label in labels[1:]:
            temp = np.empty([lat,lon])
            temp[CC2D_clean == label] = GPM_step[CC2D_clean == label] 
            # remove intense dark spot
            # Noise = np.quantile(temp[temp > 0 ], [0.9999,1],interpolation = 'nearest')
            hist, bin_edges = np.histogram(temp[temp > 0 ])
            if any(hist[-3:] <3):
                _noise = bin_edges[-3:]
                noise = np.min(_noise[hist[-3:] <4])                
                temp[temp >= noise ] = 0   
                # print('clean_' + str(Step) + 'label_'+str(label))
            Object[CC2D_clean == label] = temp[CC2D_clean == label]  
        
        DATA_SATClean[:,:,Step] =  Object        

    DATA_SATClean[np.isnan(DATA_SATClean)] = 0 #remove isnan isinf
    DATA_SATClean[np.isinf(DATA_SATClean) ] = 0 
    return  DATA_SATClean

def sim_time(start_time):
    elapsed_time = time.time() - start_time
    Time = time.strftime("%H:%M:%S", time.gmtime(elapsed_time))
    print(Time)
    return Time
    

# Storm analysis # 1. Connected labelling component T1
###############################################################################      

def Storm_list(P,GPM_event,LAT,LON,connec=6): # connec=6 # only 26, 18, and 6 are allowed  6 is the strongest connection
    #######################################  
    # CC3D_regions = connected3D(GPM_event,T,connec)
    BN = np.array(GPM_event > P['T'], dtype=int) 
    CC3D = connected_components(BN,connectivity= connec) 
    CC3D_regions = skm.measure.regionprops(CC3D,GPM_event)
    #Statistics
    Areas,Labels,Region_box = [],[],[]
    for props in CC3D_regions:
        Areas.append(props.area)
        Labels.append(props.label)
        Region_box.append(props.bbox)
    
    # indetify storm within  the domain
    Region_box = np.array(Region_box) 
    minY,minX,minZ = [Region_box[:,0],Region_box[:,1],Region_box[:,2]]
    maxY,maxX,maxZ = [Region_box[:,3]-1,Region_box[:,4]-1,Region_box[:,5]-1]    
    
    minLat = np.array([LAT[i] for i in minY])
    minLon = np.array([LON[i] for i in minX])
    maxLat = np.array([LAT[i] for i in maxY])
    maxLon = np.array([LON[i] for i in maxX])   
        
    [Xmin_lbm,Xmax_lbm,Ymin_lbm,Ymax_lbm] = P['Area_basin']  
    
    Duration = Region_box[:,5] - Region_box[:,2]
    Extension = (Region_box[:,4]-Region_box[:,1]) * (Region_box[:,3]-Region_box[:,0])         
       
    NumCLust = np.array([Areas, Labels,range(len(Labels)),Duration,Extension]).T  
    NumCLust = NumCLust[(NumCLust[:,3]  > 3 ) & (NumCLust[:,0] > 50)] # Select storm within basin and bigger than 3  
    NumCLust = NumCLust[NumCLust[:,0].argsort()[::-1]]   

    return [CC3D, CC3D_regions, NumCLust]

###############################################################################
#%%               #############  main process ##########
###############################################################################

def runner(Dates,MATRIX,P):      
    
    ############    1.   Storm analysis
    ###############################################################################
    print('ST-CORA analysis '+ Dates[0].strftime('%m/%d/%Y') + ' - ' + Dates[-1].strftime('%m/%d/%Y'))
    Total_time = time.time()   # Total time
    [Xmin_lbm,Xmax_lbm,Ymin_lbm,Ymax_lbm] = P['boundary'] # Domain   
    lat = np.int((Ymax_lbm - Ymin_lbm) / P['pixel_value'])
    lon = np.int((Xmax_lbm - Xmin_lbm) / P['pixel_value'])
    
    # Parallel
    num_cores = P['CPU_cores']# multiprocessing.cpu_count()
    
    if not os.path.exists(P['Storm_dir']):
        os.mkdir(P['Storm_dir'])     
    if not os.path.exists(os.path.join(P['Storm_dir'],'pkl_files')):
        os.mkdir(os.path.join(P['Storm_dir'],'pkl_files'))    

    # SRE Coordinates
    LAT= np.linspace(Ymin_lbm,Ymax_lbm, num = lat)  # Array latitude
    LON= np.linspace(Xmin_lbm,Xmax_lbm, num = lon)  # Array longitude
    LAT = LAT[::-1]   
    
    DATA_SAT = MATRIX
    [lat,lon,t]= DATA_SAT.shape
    MATRIX_SAT = clean_SAT(DATA_SAT,0.1,P['Minsize']) # clean Data
    MATRIX_SAT[np.isnan(MATRIX_SAT)]=0  

    #   list of storm events using ST-CORA    
    CC3D, CC3D_regions, NumCLust = Storm_list(P,MATRIX_SAT,LAT,LON,connec=6)
    
    #1. Storm segmentation
    ###############################################################################    
    print('Processing KDE segmentation...')
    if P['kernel'] == 1: 
        start_time = time.time() 
        # [X,Y,Z,R]=[[],[],[],[]]  
        Inputs = range(len(NumCLust))
        # X,Y,Z,R = stcora.KED_segmentation(NumCLust,CC3D_regions,P['T2'])
        StomArrays = Parallel(n_jobs=num_cores, backend="multiprocessing")(delayed( stcora.KED_segmentation)(No_c,NumCLust,CC3D_regions,P['T2'])
                                                        for No_c in Inputs)

        [X,Y,Z,R] = [[],[],[],[]]
        for st in StomArrays:
            # i = results[0]
            X.append(st[0])
            Y.append(st[1])
            Z.append(st[2])
            R.append(st[3])   
            
        X,Y,Z,R= [np.concatenate(X),np.concatenate(Y), np.concatenate(Z),np.concatenate(R)]

        MATRIX_SAT = np.zeros([lat,lon,t])
        for i in range(len(Z)):
            MATRIX_SAT[Y[i],X[i],Z[i]] = R[i]
        
        # Update storm list
        CC3D, CC3D_regions, NumCLust = Storm_list(P,MATRIX_SAT,LAT,LON,connec=6)
        sim_time(start_time) 
    
    #2. Storm extraction
    ###############################################################################    
    start_time = time.time() 
    Inputs = range(len(NumCLust))    
    Parallel(n_jobs=num_cores, backend="multiprocessing")(delayed( stcora.ST_extraction)(No_c,P,Dates,LAT,LON,P['Psize'],MATRIX_SAT,CC3D, CC3D_regions, NumCLust)
                                                        for No_c in Inputs)
    sim_time(start_time)  
        
    print('TOTAL TIME')
    Total_time = sim_time(Total_time)  
    
    Perfom_file = open(os.path.join(P['Storm_dir'],'Performance.txt'),'a')  # error file        
    Perfom_file.write( P['Storm_dir'] + '  Time' + Total_time +' \n') 
    Perfom_file.close()

    # return Counter
