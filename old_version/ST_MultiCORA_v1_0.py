# -*- coding: utf-8 -*-
"""

######              ST-CORA             ############ 
######              version 2.1         ############
######              multi objects         ###########

@author: laverde-Barajas M.A.
@date: Aug 26 17:00:35 2019

input:
ST_CORA(GPM_event,T,T2,Prob,LAT,LON,Minsize,Psize)

Output:
EVENT = {MCtime, MCspace, Centroid, weighted_centroid, GPM_storm, Kernel, Event}
    
"""
import numpy as np
import skimage as skm 
from skimage import morphology
from statsmodels.nonparametric.kernel_density import KDEMultivariate
from cc3d import connected_components

def ST_CORA(GPM_event,T,T2,Prob,LAT,LON,Minsize,Psize):    
    
    GPM_event[np.isnan(GPM_event)]=0  
    connec=6 # only 26, 18, and 6 are allowed  6 is the strongest connection
    [lat,lon,t]=GPM_event.shape    
    ###############################################################################
    ############ 1. Clean maps (2D analysis)
    ###############################################################################
    
    BN_Clean = np.empty([lat,lon,t])
    GPM_eventZ = np.empty([lat,lon,t])
    for Step in range(t):
        
    #    Step=0
        GPM_step = GPM_event[:,:,Step]
        SEE = GPM_step > T
        
        SEE_clean = skm.morphology.closing(SEE, skm.morphology.square(3)) # close objects   to remove small dark spots
        CC2D  = skm.measure.label(SEE_clean, background=0)
        CC2D_clean = skm.morphology.remove_small_objects(CC2D, Minsize)
        GPM_eventZ[:,:,Step] = GPM_step
        BN_Clean[:,:,Step] = CC2D_clean
    
    temp = BN_Clean  > 0 
    GPM_eventClean = np.zeros([lat,lon,t])
    GPM_eventClean[temp] = GPM_eventZ[temp]
    GPM_eventClean[np.isnan(GPM_eventClean)] = 0 #remove isnan isinf
    GPM_eventClean[np.isinf(GPM_eventClean) ] = 0 
    
    ###############################################################################
    #%% 2. Event segmentation (Multivariate kernel density estimation )
    ###############################################################################
    BN = np.array(GPM_eventClean  > T, dtype=int)
    CC3D = connected_components(BN, connectivity=connec) 
    
    #Statistics
    CC3D_regions = skm.measure.regionprops(CC3D )
    Areas = [props.area for props in CC3D_regions]  
    Labels = [props.label for props in CC3D_regions]  
    
    NumCLust = np.array([Areas, Labels,range(len(Labels))])
    NumCLust =NumCLust.T
    NumCLust = NumCLust[NumCLust[:,0].argsort()[::-1]] 
    
    [X,Y,Z,R]=[[],[],[],[]]

    for No_c in range(len(NumCLust)):
        if NumCLust[No_c,0] > 50:  # clusters need to be bigger than 100
            Clus_label = NumCLust[No_c,2] 
            Cluster= CC3D_regions[Clus_label ]
            Cluste_m= CC3D == NumCLust[No_c,1]
            Rain_cluster  = GPM_eventClean[Cluste_m]
            y, x, z = Cluster.coords.T
    #        z= z* 154285714.2 # speed of ligth / pixel resolution
            Matrix = np.array([y, x, z,np.round(Rain_cluster,2)])       
            
            kde = KDEMultivariate(Matrix, 
                              var_type='cccc')
            density = 1 - kde.cdf(Matrix)        
            
            # select values based on the kernel threshold
            density[np.isnan(density)]=0
            mask = density > Prob
            X.append( x[mask] )
            Y.append( y[mask] )
            Z.append( z[mask] )
            R.append( Rain_cluster[mask] )

            kde = None
            density=None            
            
    X= np.concatenate(X)
    Y= np.concatenate(Y)
    Z= np.concatenate(Z)
    R= np.concatenate(R)
    
    GPM_thresholded = np.zeros([lat,lon,t])
    for i in range(len(Z)):
        GPM_thresholded[Y[i],X[i],Z[i]] = R[i]    
     
    ###############################################################################
    #%% 3. Extracting characteristics
    ###############################################################################   
    def histc(X, bins):
        map_to_bins = np.digitize(X,bins)
        r = np.zeros(bins.shape)
        for i in map_to_bins:
            r[i-1] += 1
        return [r, map_to_bins]    
   
    BN = np.array(GPM_thresholded  > T2, dtype=int)
    CC3D = connected_components(BN,connectivity= connec) 
    
    #Statistics
    CC3D_regions = skm.measure.regionprops(CC3D,GPM_thresholded )
    Areas = [props.area for props in CC3D_regions]  
    Labels = [props.label for props in CC3D_regions] 
    Region_box = np.array([list(props.bbox) for props in CC3D_regions] )
    
    Duration = Region_box[:,5] - Region_box[:,2]
    Extension = (Region_box[:,4]-Region_box[:,1]) * (Region_box[:,3]-Region_box[:,0]) 
    
    # select 
    NumCLust = np.array([Areas, Labels,range(len(Labels)),Duration,Extension]).T
    NumCLust = NumCLust[NumCLust[:,3].argsort()[::-1]] 
    
    EVENTS ={}
    
    for No_c in range(len(NumCLust)):          
        if NumCLust[No_c,3] > 2: # bigger than 2 hours
            Clus_label = NumCLust[No_c,2] 
            Cluster = CC3D_regions[Clus_label ]            
                        
            Rain_cluster  = GPM_thresholded[CC3D == NumCLust[No_c,1]]
            
            Region_box = Cluster.bbox
            
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
            
            weighted_centroid = np.round(Cluster.weighted_centroid).astype(int)
            weighted_centroid = [weighted_centroid[0],weighted_centroid[1],weighted_centroid[2]]
            weighted_centroid_coord = [LAT[weighted_centroid[0]],LON[weighted_centroid[1]], weighted_centroid[2] ]
            #z= z* 154285714.2 # speed of ligth / pixel resolution
            
            GPM_storm = np.empty([lat,lon,t])
            GPM_storm.fill(np.nan)
            for i in range(len(Rain_cluster)):
                GPM_storm[y[i],x[i],z[i]] = Rain_cluster[i]    
            
            event = {'MCtime':MCtime, 'MCspace': MCspace, 'MCmax': MCmax, 'MCvol': MCvol,
                     'centroid':Centroid,'centroid_coord':Centroid_coord,'weighted_centroid':weighted_centroid,
                     'weighted_centroid_coord':weighted_centroid_coord,'GPM_storm': GPM_storm,
                     'Event': Event,'Event_M': Event_M, 'Region_box': Region_box}
            
            EVENTS[No_c] = event
    
    return EVENTS