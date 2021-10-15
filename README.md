# ghc4vpc
Geostatistical Hierarchical Clustering (GHC) applied to Vertical Sand Proportion Curves (VPC)

## Data preparation
Two options:
* Use Flumy (https://flumy.minesparis.psl.eu) to generate 3D facies models, export virtual wells and VPC files
* Adapt the code to read files having another format

Currently, two VPC examples are provided:
* Loranca Basin (Spain) VPC (8 wells) : Stats1m
* Synthetic VPC (20 wells / 3 units) : 3u_20w

## Running the script
Simply execute the R code in ghc4vpc.R file.

## References

Romary T., Ors F., Rivoirard J., Deraisme J. 
Unsupervised classification of multivariate geostatistical data: Two algorithms.
Comput. Geosci. 85, 96–103 (2015).
https://doi.org/10.1016/j.cageo.2015.05.019

Bubnova, A., Ors, F., Rivoirard, J. et al.
Automatic Determination of Sedimentary Units from Well Data.
Math Geosci 52, 213–231 (2020).
https://doi.org/10.1007/s11004-019-09793-w
