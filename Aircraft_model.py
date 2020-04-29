import Engine as e
import SystemAutomaticControl_release as sys
import Aerodinamic as aero
import atmosphere as atm
import numpy as np
from math import sqrt
class Aircraft:
    def __init__(self,H,V,angle,dt):
        self.engine = e.Engine()
        self.atmos = atm.Atmosphere()
        self.atmos.set_H(H)
        ro = self.atmos.get_density()
        g = self.atmos.get_accel_of_gravity()
        self.aerodynamic = aero.Aerodynamics(H,V,g,ro,angle,dt)
        self.control = sys.Control()
        
        
    def run(self,v_zad,gama_zad,ny_zad,dt):
       angl,X,w,V,n,e =  self.aerodynamic.Get_data()
       self.atmos.set_H(X[1])
       ro = self.atmos.get_density()
       g = self.atmos.get_accel_of_gravity()
       self.control.set_data(ny_zad,n[1],gama_zad,angl[0],w[0])
       eliv,eleron = self.control.get_data()
       V_m=sqrt(V[0]*V[0]+V[1]*V[1]+V[2]*V[2])
       self.engine.Set_data(V_m,v_zad,ro)
       P,M = self.engine.Get_data()
       self.aerodynamic.Set_data(ro,eliv,0,eleron,P,g)
       return n
H=2000
angle=np.array([0,0,0])
V = np.array([50,0,0])
plane = Aircraft(H,V,angle,0.02)
T=10
t=0
while(T>t):
    X=plane.run(75,0.0,2.0,0.02)
    t+=0.02
print(X)