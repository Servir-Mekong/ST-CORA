# ST-CORA
Spatiotemporal Object-based Rainfall Analysis 


## Dependencies

numpy
skimage 
skimage
statsmodels
cc3d 


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
