# Guide to Python code for SOTE 2.0

The code contained in this repository corresponds to the SOTE 2.0 model, as described in [] (2021) (LINK).  

This code can be used to study the risk of soil degradation (as measured by changes in saturated hydraulic conductivity, $K_s$) for soils with different weight functions. The code also compares the risk computed using weight functions to that computed using the McNeal function (McNeal, 1968).

For more information regarding the default SOTE 2.0 code, as contained in this repository, users are referred to section 3.1.1 of Kramer et al. (2022). 

1. In order to run the code, execute run_SOTE.py

2. The following input parameters can be easily adjusted:

   - Soil type (must choose from following weight functions: "class_1_irr", "class_1_rev", "class_2_irr", "class_2_rev", "class_3_irr", and "class_3_rev")
   - Simulation length (*years*, float)
   - Number of simulations (*runs*, integer)
   - Rainfall season length (*days*, integer)
   - Probably of rain during rainy season (*rain_prob*, float)
   - Rainfall event height (shape and scale, according to Weibull distribution, float)
   - Irrigation water salinity (*Cirr*, mmol_c/L, float)
   - Irrigation water sodicity fraction (*Eirr*, nondimensional, float)

3. The following parameters must be adjusted through the function_McNeal.py and function_SOTE.py files:

   - Initial soil solution salinity (*C_init*, mmol_c/L, float)

   - Initial soil sodicity fraction (*E_init*, nondimensiona, floatl)

   - Initial relative soil water content (*s_init*, nondimensiona, floatl)

   - Rain water salinity (*Crain*, mmol_c/L, float)

   - Rain water sodicity fraction (*Erain*, nondimensional, float)

   - Minimum evapotranspiration rate (*ET_w*, mm, float)

   - Time step (*dt*, day, float)
   - Ratio of Irrigation to ET: (*ET_ratio*, integer, float)

4. Soil properties can be controlled through the soil_properties.py file.

## References
