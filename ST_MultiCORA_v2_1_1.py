# -*- coding: utf-8 -*-
"""

######              ST-CORA             ############ 
######              version 2.2.1         ############
######              multi objects         ###########

@author: laverde-Barajas M.A.
@date: Aug 26 17:00:35 2019

input:
ST_CORA(No_c,P,Dates,LAT,LON,Psize,GPM_thresholded,CC3D, CC3D_regions, NumCLust)

Output:
EVENT = {MCtime, MCspace, Centroid, weighted_centroid, GPM_storm, Kernel, Event}


GPM_event= DATA_SATclean
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
import numpy as np
import skimage as skm 
from skimage import morphology
# from statsmodels.nonparametric.kernel_density import KDEMultivariate
from cc3d import connected_components
import pickle as pkl
import os
#from compress_pickle import dump, load
from scipy import stats

def histc(X, bins):
       map_to_bins = np.digitize(X,bins)
       r = np.zeros(bins.shape)
       for i in map_to_bins:
           r[i-1] += 1
       return [r, map_to_bins]  
    
def condi_area(P,Lat,Lon):
    [Xmin_lbm,Xmax_lbm,Ymin_lbm,Ymax_lbm] = P['Area_basin']
    Condi = np.asanyarray(np.where( (Lon >= Xmin_lbm ) & (Lon <= Xmax_lbm) &
                                   (Lat >= Ymin_lbm ) & (Lat <= Ymax_lbm) ) )
    Rate = len(Condi.T) /len(Lon)
    return Rate       

def connected3D(Object,T,connec):
    BN = np.array(Object  > T, dtype=int)
    CC3D = connected_components(BN,connectivity= connec) 
    CC3D_regions = skm.measure.regionprops(CC3D,Object )    
    return CC3D_regions  
    
def check_sub(Object,T2): 
    # Region_box_sub = Cluster.bbox             
    BN = np.array(Object  > T2, dtype=int) 
    CC3D = connected_components(BN,connectivity= 26) 
    CC3D_regionsSub = skm.measure.regionprops(CC3D,Object) 
    if len(CC3D_regionsSub)< 1:
        subNumber = 0
    else: 
         # dimension of the substorms                
        Dim = np.array([props.intensity_image.shape  for props in CC3D_regionsSub] )        
        Dim = np.array([Dim[:,0]*Dim[:,1],Dim[:,2]]).T
        # Dim = Dim[Dim[:,0].argsort()[::-1]]  
        CON =  (Dim[:,0] >= 20 ) & (Dim[:,1] >= 3)       
        subNumber = Dim[CON]   
        subNumber = len(subNumber[:,0])
    return subNumber

 ###############################################################################
 #%% Event segmentation (Multivariate kernel density estimation )
    ###############################################################################
#    EVENTS ={}  
def KED_segmentation(No_c,NumCLust,CC3D_regions,T2):  
    # [X,Y,Z,R]=[[],[],[],[]] 
    # for No_c in range(len(NumCLust)):    

    # Characteristics of the Object storm 
    Clus_label = NumCLust[No_c,2] 
    Cluster = CC3D_regions[Clus_label ]              
    
    Object = Cluster.intensity_image                   
    y, x, z  = Cluster.coords.T
    Rain_cluster = Object[Object>= 1]         
    
    # 2. Connected labelling component T2
    #######################################                      
    subNumber = check_sub(Object,T2)  
     
    # 2. Kernel segmentation
    #######################################
    if subNumber >1: #find subtorm           

        Matrix = np.array([y, x, z,np.round(Rain_cluster,2)])    
        kde = stats.gaussian_kde(Matrix)
        density = kde(Matrix)      
    
        # select values based on the kernel threshold
        density[np.isnan(density)]=0
        
        # histogram               
        Hist = np.array([np.round(Rain_cluster),density]).T
        Hist = Hist[Hist[:,0].argsort()[::-1]]                  
        # prob_ = Hist[Hist[:,0] == T2,:]  
        closeT2 = np.abs(Hist[:,0] - T2)
        prob_ = Hist[closeT2==closeT2.min(),:]
        prob = np.max(prob_[:,1]) 
        
        mask = density > prob
        X = x[mask] 
        Y = y[mask] 
        Z = z[mask] 
        R = Rain_cluster[mask] 

        kde,density = [None,None]                
        print('kernel_' + str(No_c))
        
    else:  
        X = x 
        Y = y 
        Z = z 
        R = Rain_cluster   

    return X,Y,Z,R
        ###############################################################################
        #%% Extracting characteristics
        ###############################################################################  
def ST_extraction(No_c,P,Dates,LAT,LON,Psize,GPM_thresholded,CC3D, CC3D_regions, NumCLust):  
    
    # for No_c in range(len(NumCLust)):          
        
    Clus_label = NumCLust[No_c,2] 
    Cluster = CC3D_regions[Clus_label ]            
                
    Rain_cluster  = GPM_thresholded[CC3D == NumCLust[No_c,1]]    
    Region_box = Cluster.bbox
    Object = Cluster.intensity_image 
    
    y, x, z  = Cluster.coords.T
    Y = np.array([LAT[i] for i in y])
    X = np.array([LON[i] for i in x])   
     
    Event = [Y,X,z,Rain_cluster]
    Event_M= [y,x,z,Rain_cluster]     
    
    Centroid =np.round( Cluster.centroid).astype(int)
    Centroid_coord = [ LAT[Centroid[0]],LON[Centroid[1]], Centroid[2] ]
    Ut = np.unique(z) 								#Time
#    Ux = np.unique(x)  							#Space
#    Uy = np.unique(y) 							#Space      
    [a,b] = histc(z,Ut)
    MCspace = np.max(a)*Psize
    MCtime = len(a)
    MCmax = np.max(Rain_cluster)
    MCvol = np.nansum(Rain_cluster*0.001*Psize*1000000*1e-9)
    
    Evn_date = Dates[Ut.astype(int)]
    weighted_centroid = np.round(Cluster.weighted_centroid).astype(int)
    weighted_centroid = [weighted_centroid[0],weighted_centroid[1],weighted_centroid[2]]
    weighted_centroid_coord = [LAT[weighted_centroid[0]],LON[weighted_centroid[1]], weighted_centroid[2] ]
    #z= z* 154285714.2 # speed of ligth / pixel resolution
    
    Object = Cluster.intensity_image  
    
    event = {'Evn_date':Evn_date[0],'MCtime':MCtime, 'MCspace': MCspace, 'MCmax': MCmax, 'MCvol': MCvol,
             'centroid':Centroid,'centroid_coord':Centroid_coord,'weighted_centroid':weighted_centroid,
             'weighted_centroid_coord':weighted_centroid_coord,
             'Event': Event,'Event_M': Event_M,'Region_box': Region_box,'Object':Object}
    
    Date_init = str(Evn_date[0]).replace(' ','_')
    Date_init = Date_init.replace(':','')
    
# #            EVENTS[No_c] = event
    output = open(os.path.join(P['Storm_dir'],'pkl_files','MCS_' + Date_init + '.pkl'), 'wb')
    pkl.dump(event, output)
    output.close()
    
    # # Save compress pickle
    # Pick_name = os.path.join(P['Storm_dir'],'MCS_' + Date_init + '.gz')
    # dump(event, Pick_name, compression="lzma", set_default_extension=False)                 
    
    print('stcora. {0} from {1}'.format(No_c,len(NumCLust)))    
#    return EVENTS
