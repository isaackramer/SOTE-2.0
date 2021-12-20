"""Class used to run the system. Includes integration of ds/dt term.
Includes soil degradation related to declines in hydraulic conductivity."""

import numpy as np
import sympy as sym
import soil_properties

class SOTE(object):
    def __init__(self, 
                 Cirr=0.0,
                 Eirr=0.0,
                 Crain=0.0,
                 Erain=0.0, 
                 C_init=0.0, 
                 E_init=0.0, 
                 s_init=0.0, 
                 Zr=0.0, 
                 soil_type="soil", 
                 C_max=0,
                 rain_prob=0.0, 
                 days=0.0, 
                 shape=0.0, 
                 scale=0.0,
                 ET_w=0.1,                 
                 ET_ratio=0.0, 
                 dt=0.0, 
                 runs=0.0, 
                 hys_dens=0.0):
        
        # dictionary with input parameters
        self.par = {}
        self.par.update(Cirr=Cirr, 
                        Eirr=Eirr, 
                        Crain=Crain, 
                        Erain=Erain,
                        C_init=C_init, 
                        E_init=E_init, 
                        s_init=s_init,
                        rain_prob=rain_prob, 
                        days=days, 
                        shape=shape, 
                        scale=scale,
                        ET_w=ET_w,
                        ET_ratio=ET_ratio, 
                        dt=dt,
                        runs=runs, 
                        hys_dens=hys_dens)
        
        # create dictionary with soil properties and load weight functions
        (self.soil_parms, 
         self.weights_C, 
         self.weights_E) = soil_properties.soil_props(soil_type=soil_type,
                                                      Zr=Zr,
                                                      C_max=C_max)

        # initial values for saturated hydrualic conductivity and water dynamics
        self.Ksat =  np.full((self.par['runs'],),
                             self.soil_parms['Ks'])
        self.soil_parms['nZr'] = self.soil_parms['n']*self.soil_parms['Zr']
        
        # initate functions
        self.functions_def()

    def functions_def(self):
        # define parts needed to calculate dE/dt
        self.K1 = lambda q, E, w: (self.dgsdC(q/w, E) * self.dsdt * 
                                   np.power(q,2)) * (self.soil_parms['nZr'])
        self.K2 = lambda q, E, w: self.dgsdC(q/w, E) * q * self.dqdt * w
        self.K3 = lambda w: (self.Cinput * self.Einput * 
                             self.s_In_rate * np.power(w,2))
        self.K4 = lambda q, E, w: self.Lw * q * w * self.gs(q/w, E)
        self.K5 = lambda q, E, w: self.dqdt * np.power(w,2) * self.gs(q/w, E)
        self.K6 = lambda w: (self.soil_parms['CEC'] * 
                             self.soil_parms['Msoil'] * np.power(w,2))
        self.K7 = lambda q, E, w: self.dgsdE(q/w, E) * q * np.power(w,2)

        # define Gapon equation and partial derivatives
        self.set_g()
        
        # height of preisach triangles equal to number of runs
        z_values = np.arange(0, self.par['runs'])
        
        # salinity triangle
        cmax = self.soil_parms['C_max']
        hps = self.par['hys_dens']
        C_beta_grid_values = np.linspace(0, cmax, hps)
        C_alpha_grid_values = C_beta_grid_values[::-1]
        (self.beta_grid_C,
         self.alpha_grid_C,
         self.runs_grid) = np.meshgrid(C_beta_grid_values,
                                       C_alpha_grid_values,
                                       z_values, indexing='ij')
        self.triangle_C = np.where(self.alpha_grid_C>=self.beta_grid_C,
                                   True,
                                   False)
        
        # sodicity triangle
        E_beta_grid_values = np.linspace(-1, 0, hps)
        E_alpha_grid_values = E_beta_grid_values[::-1]
        (self.beta_grid_E,
         self.alpha_grid_E,
         self.runs_grid) = np.meshgrid(E_beta_grid_values,
                                       E_alpha_grid_values,
                                       z_values, indexing='ij')
        self.triangle_E = np.where(self.alpha_grid_E>=self.beta_grid_E,
                                   True,
                                   False)
        
    def derivs(self):
        # define dq/dt, dE/dt, and ds/dt
        self.dqdt = self.input_salts_rate - self.output_salts_rate
        self.dEdt = lambda q, E, w: ((self.K1 (q, E, w) - 
                                      self.K2(q, E, w) + 
                                      self.K3(w) - 
                                      self.K4(q, E, w) - 
                                      self.K5(q, E, w))/(self.K6(w) + 
                                                         self.K7(q, E, w)))  
        self.dsdt = ((self.s_In_rate - self.s_Out_rate)/
                     (self.soil_parms['n']*self.soil_parms['Zr']))

        # put all derivatives into a single function
        self.rhs = lambda data: np.array([self.dqdt,
                                          self.dEdt(data[0], 
                                                    data[1], 
                                                    (data[2]*
                                                     self.soil_parms['n']*
                                                     self.soil_parms['Zr'])),
                                          self.dsdt])

    def set_g(self):
        # define gapon exchange and partials
        C, E = sym.symbols('E C')
        gs_eq = lambda C, E: (2.0 / 
                              (1.0 + (1.0 + 8.0 * 
                                      (self.soil_parms['Kg']) ** 2 * 
                                      C * (1.0 / E - 1.0) ** 2) ** 0.5))
        self.gs = lambda C, E: gs_eq(C, E)
        diffgsdE = sym.diff(self.gs(C, E), E)
        diffgsdC = sym.diff(self.gs(C, E), C)
        self.dgsdE = sym.lambdify((C, E), diffgsdE, "numpy")
        self.dgsdC = sym.lambdify((C, E), diffgsdC, "numpy")

    def set_K(self, C, E, C_old, E_old):
        ## calculate relative hydralulic conductivity (Kramer et. al., 2020)
        # increase in C
        self.triangle_C = np.where(np.logical_and(C > C_old, 
                                                  C > self.alpha_grid_C),
                                   True, 
                                   self.triangle_C) 
        # decrease in C
        self.triangle_C = np.where(np.logical_and(C < C_old,
                                                  C < self.beta_grid_C),
                                   False, 
                                   self.triangle_C) 
        # increase in E
        self.triangle_E = np.where(np.logical_and(E > E_old,
                                                  E > self.alpha_grid_E),
                                   True,
                                   self.triangle_E)
        # decrease in E
        self.triangle_E = np.where(np.logical_and(E < E_old,
                                                  E < self.beta_grid_E),
                                   False,
                                   self.triangle_E)
            
        # C_weights depend on value of E
        E_index = -np.rint(np.minimum(E*(self.par['hys_dens']-1)/
                                      self.soil_parms['E_max'],
                                      self.par['hys_dens']-1)).astype(int)
        C_weights = self.weights_C[:,:,E_index]
        
        # E_weights depend on value of C
        C_index = np.rint(np.minimum(C*(self.par['hys_dens']-1)/
                                     self.soil_parms['C_max'],
                                     self.par['hys_dens']-1)).astype(int)
        E_weights = self.weights_E[:,:,C_index]
        
        # total 'mass' of the weight fucntions
        mass = (np.nansum(C_weights, (0,1)) + 
                np.nansum(E_weights, (0,1)))
        
        # normalize weight functions (fully saturated state corresponds to 1)
        C_weights = np.where(mass > 0, 
                             C_weights/mass, 
                             C_weights)
        E_weights = np.where(mass > 0, 
                             E_weights/mass, 
                             E_weights)
        
        # integrate over weight functions
        weighted_prism = (self.triangle_C*C_weights + 
                          self.triangle_E*E_weights)
        factor = np.nansum(weighted_prism, (0,1)) 
        
        # determine Ks value based on integrated value; 
        # if mass == 0, retain previous value
        Ks = np.where(mass>0, 
                      factor*self.soil_parms['Ks'], 
                      self.Ksat)
        return Ks
        
        
    def input_q(self):
        # calculate the salinity of the input water
        irr_salts_rate = self.irrigation_rate*self.par['Cirr']
        rain_salts_rate = self.rain_rate*self.par['Crain']
        self.input_salts_rate = irr_salts_rate + rain_salts_rate
        with np.errstate(divide='ignore', invalid='ignore'):
            self.Cinput = np.where(self.s_In_rate == 0, 
                                   0, 
                                   self.input_salts_rate/(self.s_In_rate))

    def input_sodicity(self):
        # calculate the sodicity of the input water
        irrigation_sodicity = self.irrigation_rate*self.par['Eirr']
        rain_sodicity = self.rain_rate*self.par['Erain']
        self.input_sodium = irrigation_sodicity + rain_sodicity
        with np.errstate(divide='ignore', invalid='ignore'):
            self.Einput = np.where(self.s_In_rate == 0, 
                                   0, 
                                   self.input_sodium/(self.s_In_rate))

    def water_input(self):
        # calculate total input water and salinity and sodicity of input water
        self.irrigation_rate = self.Irr
        self.s_In_rate = self.irrigation_rate + self.rain_rate
        self.input_q()
        self.input_sodicity()

    def rain_height(self):
        # calculate rain height based on exponential distribution
        rand = np.random.rand(self.par['runs'],)
        event_height = np.where(rand <= self.par['rain_prob']*self.par['dt'],
                                (np.random.weibull(self.par['shape'], 
                                                  size = self.par['runs'])*
                                 self.par['scale']), 
                                0)
        self.rain_rate = event_height/self.par['dt']
        self.rain_rate = np.where(self.rain_rate>self.Ksat,
                                  self.Ksat,
                                  self.rain_rate)

    def water_loss(self, s, C):
        # calculate water loss 
        self.ET_act = np.piecewise(s,
                                   [s <= self.soil_parms['s_h'],
                                    (s <= self.soil_parms['s_w']) & (s > self.soil_parms['s_h']),
                                    (s <= self.soil_parms['s_bal']) & (s > self.soil_parms['s_w'])],
                                    [lambda s: 0,
                                     lambda s: (self.par['ET_w']*
                                                ((s - self.soil_parms['s_h'])/
                                                 (self.soil_parms['s_w'] -
                                                  self.soil_parms['s_h']))),
                                     lambda s: (self.par['ET_w'] + 
                                                (self.ETmax - self.par['ET_w']) *
                                                ((s - self.soil_parms['s_w'])/
                                                 (self.soil_parms['s_bal'] - 
                                                  self.soil_parms['s_w']))),
                                     lambda s: self.ETmax])

        # calculate leaching
        self.Lw = self.Ksat*np.power(s,(self.soil_parms['c']))

        # calculate total loss rate
        self.s_Out_rate = self.ET_act + self.Lw

        # calculate the loss of salts
        self.output_salts_rate = (self.Lw)*C

    def water_net(self, s, C, E):
        # net change in water content
        self.water_input()
        self.water_loss(s, C)

    def rk4_step(self, data):
        # Integrate using Runge-Kutta
        dt = self.par['dt']
        k1 = dt * self.rhs(data)
        k2 = dt * self.rhs(data + 0.5 * k1)
        k3 = dt * self.rhs(data + 0.5 * k2)
        k4 = dt * self.rhs(data + k3)
        return np.array([data + (k1 + 2.0 * (k2 + k3) + k4) / 6.0])