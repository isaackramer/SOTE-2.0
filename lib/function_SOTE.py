"""Function for running SOTE using the Preisach Framework"""
import numpy as np
import pandas as pd
import class_SOTEplus


# Preisach framework function
def SOTEplus(soils,
             years,
             runs,
             rain_prob,
             shape,
             scale,
             days,
             Eirr,
             Cirr):

    # time conditions
    years = years # simulation length
    dt = 0.1 # time step (in days)
    time_array = np.arange(0, years*365+dt, dt)

    # set random seet
    np.random.seed(30)

    # instantiate class and function
    eq = class_SOTEplus.SOTE(Cirr=Cirr,
                             Eirr=Eirr,
                             Crain=0.0,
                             Erain=0.0,
                             C_init=20.0,
                             E_init=0.02,
                             s_init=0.3,
                             Zr=350,
                             soil_type=soils,
                             C_max = 100,
                             rain_prob=rain_prob,
                             days=days,
                             shape=shape,
                             scale=scale,
                             ET_w=.1,
                             ET_ratio=1.1,
                             dt=dt,
                             runs=runs,
                             hys_dens=100)

    # bounds of rainy season
    rain_start = 365 - eq.par['days']
    rain_mid = 365 - eq.par['days']/2
    sin_shift = rain_mid - 273.5
    rain_end = 365

    # inital values based on class initalization
    s0 = eq.par['s_init']          # soil moisture [nondimensional]
    w0 = s0 * eq.soil_parms['nZr'] # water content [mm]
    C0 = eq.par['C_init']          # salt concentration [mmol_c/L]
    q0 = eq.par['C_init']*w0       # salt mass content [mmol_c]
    E0 = eq.par['E_init']          # exchangeable sodium fraction [nondimensional]
    Krel0 = 1                      # relative Ks [nondimensional]

    # we will run an ensemble of N=runs simulations at once
    # vectorize inital values
    s = np.full((eq.par['runs'],), s0)
    w = np.full((eq.par['runs'],), w0)
    q = np.full((eq.par['runs'],), q0)
    E = np.full((eq.par['runs'],), E0)

    # array to record average results (t, s, C, E, Krel)
    results_avg = np.array([[time_array[0], s0, C0, E0, Krel0]])

    # array to record Ks at each time step
    Krel = np.full((eq.par['runs'],), Krel0)
    Krel_results = np.array([Krel])

    # main loop
    for tim in time_array[1:]:
        # keep old values of C and E for Ks calculation
        C_old = q/w
        E_old = E

        # ET_max dependent on time of year
        eq.ETmax = 2.0 * np.sin((tim - sin_shift) * ((2 * np.pi) / 365)) + 5

        # stochastic rainfall
        DOY = (np.around(tim,1) % 365) # day of year since 1 January
        if DOY == 0:
            days = eq.par['days']
            rain_start = 365 - days
            rain_end = rain_start + days
        if (DOY >= rain_start and DOY < rain_end):
            # if True --> rainy season
            eq.rain_height()
            eq.Irr = 0.0
        else:
            # if False --> dry season
            eq.event_height = np.zeros(eq.par['runs'],)
            eq.rain_rate = np.zeros(eq.par['runs'],)
            eq.Irr = eq.par['ET_ratio'] * eq.ETmax

        # calculates net change in water content + salinity and sodicity of input and output water
        eq.water_net(s, q/w, E)

        # integration
        eq.derivs()
        new_values = eq.rk4_step([q, E, s])
        q, E, s = new_values[0, :]

        # water content must be less than 1 and greater than 0
        s = np.where(s > 1.0, 1.0, s)
        s = np.where(s < 0.0, eq.soil_parms['s_h'], s)
        w = s * eq.soil_parms['nZr']

        # calculate new Ks, based on change in C and E
        eq.Ksat = eq.set_K(q/w, -E, C_old, -E_old)

        # update array to track results (t, s, C, E, Krel)
        step = np.array([[tim,
                          np.mean(s),
                          np.mean(q/w),
                          np.mean(E),
                          np.mean(eq.Ksat/eq.soil_parms['Ks'])]])
        results_avg = np.concatenate((results_avg, step))
        Krel_results = np.concatenate((Krel_results,
                                       np.array([eq.Ksat])/eq.soil_parms['Ks']))


    # save output
    results_avg = pd.DataFrame(results_avg,
                                columns = ['time', 'water', 'salinity',
                                          'sodicity','hyd_cond'])
    results_avg.to_csv("output_csvs/deg_"+soils+'.csv')
    ends = np.array([s, q/w, E, eq.Ksat/eq.soil_parms['Ks']])
    np.savetxt('output_csvs/deg_'+soils+'_ends.csv',
                ends,
                delimiter=",")
    np.savetxt('output_csvs/deg_'+soils+'_Ks.csv',
                Krel_results,
                delimiter=",")
