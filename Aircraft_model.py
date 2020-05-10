import Engine as e
import SystemAutomaticControl_release as sys
import Aerodinamic as aero
import atmosphere as atm
import numpy as np
from math import sqrt
import matplotlib.pyplot as plt
class Aircraft:
    def __init__(self,H,V,angle,dt):
        self.engine = e.Engine()
        self.atmos = atm.Atmosphere()
        self.atmos.set_H(H)
        ro = self.atmos.get_density()
        g = self.atmos.get_accel_of_gravity()
        self.aerodynamic = aero.Aerodynamics()
        self.control = sys.Control()
        
        
    def run(self,v_zad,gama_zad,ny_zad,dt):
       n,alpha,betta,V_a,w_a,V,w,angle,X =  self.aerodynamic.Get_data()
       self.atmos.set_H(X[1])
       ro = self.atmos.get_density()
       g = self.atmos.get_accel_of_gravity()
       self.control.set_data(ny_zad,n[1],gama_zad,angle[0],w[0])
       eliv,eleron = self.control.get_data()
       V_m=np.sqrt(np.dot(V,V))
       P = 1
       eliv = 0.0
       self.engine.Set_data(V_m,v_zad,ro)
       P,M = self.engine.Get_data()
       P = np.array([5000,0,0])
       eliv = -0.17
       eleron = 0.05
       self.aerodynamic.Set_data(P,dt,eliv,0,eleron,ro,g)
       return angle[0]
H=2000
angle=np.array([0,0,0])
V = np.array([50,0,0])
plane = Aircraft(H,V,angle,0.02)
T=10
t=0
X=[]
TT=[]
while(T>t):
    x=plane.run(55,0.0,0.0,0.02)
    t+=0.02

    X.append(x)
    TT.append(t)
#X=np.array(X)
plt.plot(TT,X)
plt.grid()
plt.show()