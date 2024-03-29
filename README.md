# ST-CORA with KDE storm segmentation
Spatiotemporal Object-based Rainfall Analysis 
* 3D and 4D object-based storm analysis 
* 4D analysis using a Multivariate kernel density estimation (KDE)
* multicore analysis

## Dependencies

- numpy==1.16.2  
- pandas==0.25.0
- scikit-image==0.14.2
- statsmodels==0.10.1 
- pickleshare==0.7.5
- datetime == 4.0.1
- skimage == 0.17.2
- scipy
-  cc3d  (https://github.com/seung-lab/connected-components-3d)
-  matplotlib
-  joblib
-  multiprocessing

@author: M. Laverde-Barajas
@email: m.laverde@un-ihe.org
        mlaverdeb@gmail.com

@company SERVIR- MEKONG program ;
        IHE Delft Institute for Water Education;
        Delft University of Technology        
## Paramaters
- T =  # wet values
- T2 =   # delineation
- Minsize =  # noise value
- kernel =   1= kernel segmentation 4D 0= 3D
- Psize =    # Min Object size
- pixel_value =  # resolution

- StartTime = datetime(YYY, M, D, h, 0, 0) #  (GMT)
- EndTime = datetime(YYY, M, D, h, 0, 0) 

- MATRIX = 3D rainfall matrix
- boundary = [Xmin_lbm,Xmax_lbm,Ymin_lbm,Ymax_lbm] 


REFERENCES

ST-CORA
Laverde-Barajas, M., Corzo, G., Bhattacharya, B., Uijlenhoet, R., & Solomatine, D. P. (2019). Spatiotemporal analysis of extreme rainfall events using an object-based approach. In Spatiotemporal Analysis of Extreme Hydrological Events (pp. 95-112). Elsevier.

ST-CORA with Multivariate kernel density estimation (KDE)
Laverde-Barajas, M., Corzo, G. A., Poortinga, A., Chishtie, F., Meechaiya, C., Jayasinghe, S., ... & Solomatine, D. P. (2020a). St-corabico: A spatiotemporal object-based bias correction method for storm prediction detected by satellite. Remote Sensing, 12(21), 3538.
