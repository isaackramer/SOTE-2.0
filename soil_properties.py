import numpy as np

def soil_props(soil_type, Zr, C_max):
    """Parameters c, Ks, n, s_h, s_w, s_bal, s_fc, bulk_d:
    Laio et al., 2001, Plants in water-controlled ecosystems:"""

    # default soil parameters
    params = {'c':6.8, 'Ks':1000.0, 'n':0.43, 's_h':0.14,
              's_w':0.18, 's_bal':0.46, 's_fc':0.56, 'bulk_d':1.5,
              'CEC': 150}

    soil_dict = {**params}

    # choose weight function
    weights_C = np.load('weights/'+soil_type+'_C.npy')
    weights_E = np.load('weights/'+soil_type+'_E.npy')

    # set maximum C and E value
    E_max = 1.0

    # other soil properties
    gapon = 0.01475
    mass = soil_dict['bulk_d']*Zr
    soil_dict.update(Kg=gapon,
                     Zr = Zr,
                     Msoil = mass,
                     C_max = C_max,
                     E_max = E_max)

    return soil_dict, weights_C, weights_E
