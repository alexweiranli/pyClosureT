import numpy as np
import scipy

## Parameter values for three diffusion geometries 
# A values from Dodson, 1973
A_geom = {'sphere':55,'cylinder':27,'plane':8.7}
# M values read off Fig. 2 of Ganguly & Tirone, 1999 
M_geom = {'sphere':0.4,'cylinder':0.8,'plane':1.1}

class cal_t:

    def __init__(self,E,D0,s,a,geometry):
        self.E = E      # activation energy, in kJ/mol
        self.D0 = D0    # pre-exponential factor, in m2/s
        self.s = s      # cooling rate, in C/year 
        self.a = a      # crystal radius, in micron

        geometry = geometry.lower()
        if geometry in list(A_geom):
            self.geometry = geometry
        else:
            print('Check input "geometry" name')

        self.A = A_geom[geometry]
        self.M = M_geom[geometry]
    
    ## calculation of mean closure temperature (Dodson, 1973)
    def tc_mean(self): 
        R = 8.314
        def Tc_eq(x):
            y = (self.E*1e3)/(R*x) - np.log(self.A*R*(x**2)*self.D0/((self.E*1e3)*(self.s/(3600*24*365))*((self.a*1e-6)**2)))
            return y
        Tc = scipy.optimize.fsolve(Tc_eq,1200) -273.15 
        return Tc[0]   # in degree Celsius 

    ## calculation of peak temperature T0 (Faak et al. 2014, Eq. (6))
    def tpeak(self): 
        R = 8.314
        def t0_eq(x):
            D_ = self.D0*np.exp(-(self.E*1000)/(R*x))
            y = np.sqrt(self.M*(self.a*1e-6)**2*(self.E*1e3)*self.s/(3600*24*365)/(D_*R)) - x
            return y
        T0 = scipy.optimize.fsolve(t0_eq,300) - 273.15 
        return T0[0]      # in degree Celsius 
    