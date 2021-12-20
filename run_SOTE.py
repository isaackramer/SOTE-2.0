"""SOTE 2.0"""

# import necessary functions
import sys
sys.path.insert(0, './lib')
import function_SOTE
import function_McNeal
import concurrent.futures


# soils to be tested. choose from:
    # "class_1_rev", "class_1_irr",
    # "class_2_irr", "class_2_rev",
    # "class_3_rev", "class_3_irr"
soils = ["class_1_irr", "class_1_rev",
         "class_2_irr", "class_2_rev",
         "class_3_irr", "class_3_rev"]

# number of soils
num_soils = len(soils)

# simulation length
years = [1.5] * num_soils

# number of iterations
runs = [100] * num_soils

# rainfall properties
rain_prob = [0.30] * num_soils
shape = [0.47] * num_soils
scale = [3.25] * num_soils
days = [130] * num_soils

# water quality parameters
Eirr = [0.5] * num_soils
Cirr = [15] * num_soils


if __name__ ==  '__main__':
    # run SOTE 2.0 model using parallel processes
    with concurrent.futures.ProcessPoolExecutor() as executor:
        executor.map(function_SOTE.SOTEplus,
                      soils,
                      years,
                      runs,
                      rain_prob,
                      shape,
                      scale,
                      days,
                      Eirr,
                      Cirr)

    # compare SOTE 2.0 to McNeal function
    function_McNeal.McNeal(soils[0],
                            years[0],
                            runs[0],
                            rain_prob[0],
                            shape[0],
                            scale[0],
                            days[0],
                            Eirr[0],
                            Cirr[0])
