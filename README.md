# ST-CORA
Spatiotemporal Object-based Rainfall Analysis 


## Dependencies

numpy
skimage 
skimage
statsmodels
cc3d 


 ##   
input:
ST_CORA(GPM_event,T,T2,Prob,LAT,LON,Minsize,Psize)

T and T2 = convective thresholds
Prob = kernel thresholds
LAT = Array Latitude
LON = Array Longitude
Minsize = Filtering size object filtering
Psize =  Satellite resolution

Output:
EVENT = {MCtime, MCspace, Centroid, weighted_centroid, GPM_storm, Kernel, Event}
