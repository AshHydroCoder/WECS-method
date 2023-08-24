# WECS-method
WECS method or Hydest method is an empirical method of determining the peak discharge of a return period for given value of Drainage Area (Sq. Km.) and Peak Rainfall of same return period. This code is written to determine the coefficient and the exponents associated with the drainage area and rainfall. 
The code utilizes Least Square Technique to determine the parameters. 
It started with collecting the peak rainfall data and instanteneous discharge for several years and fitting it under the Gumbel distribution to determine the peak rainfall and instatenteous discharge for 50 and 100 year return period. The drainage basin was also delineated using NASA Earth data 30 m DEM.
This gave us the input parameters i.e., area, rainfall and discharge. 
The matrix formulation is shown in image. The final formula is in the form:
Q50 = α* R50^β* A^Γ


