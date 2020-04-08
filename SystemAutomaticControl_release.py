import numpy as np
from model_function_upr import *
from matplotlib import pyplot as plt

#eps_ny = 0.01       # |(ny_now - ny_spec)| > eps_ny
dt = 0.02
k_elev = 0.001
k_eleron = 1.5
  

class Control:

    # Функция получения положения руля высоты и текущей перегрузки 
    #   ny_spec - заданное значение перегрузки 
    def get_elev_and_ny_new (self, ny_spec) :
             
        self.ny_now = 0 # начальное значение перегрузки
        self.elev_new = 0 # начальное положение угла руля высоты
        
        self.t = 0

        #while ((abs(self.ny_now - self.ny_spec) > eps_ny) and (t <= 10)) :
        while (self.t <= 10) :
        
            self.ny_now = get_ny(self.elev_new) / ACCELERATION_OF_GRAVITY # получение значения перегрузки в данный момент времени
            self.elev_new = self.elev_new + k_elev * (ny_spec - self.ny_now) # получение отклонения рулей высоты
            #elev_new =  0.1 * (ny_spec - ny_now) # получение отклонения рулей высоты

            if (self.elev_new * 180.0/np.pi >= 26) :
                self.elev_new = 26/180.0*np.pi

            if (self.elev_new * 180.0/np.pi <= -28) :
                self.elev_new = -28/180.0*np.pi
        

            self.t = self.t + dt
            
        return self.elev_new

    # Функция получения положения элерона и крена самолета 
    #   gamma_spec - заданное значение крена самолета (желаемое)
    def get_GammaAngle_and_eleron_now (self, gamma_spec) :
        
        self.eleron_now = 0.0 # начальное положение элерона
        self.gamma_now = 0 # начальный угол крена самолета
        
        
        self.t = 0

        #while ((abs(gamma_spec - gamma_now) > 0.01) and (t <= 20)) :
        while  (self.t <= 20):
            
            self.gamma_now, self.w0 = get_gamma(self.eleron_now)

            self.eleron_now = k_eleron*(gamma_spec - self.gamma_now) - 3*self.w0

            if (self.eleron_now * 180.0/np.pi >= 20) :
                self.eleron_now = 20/180.0*np.pi

            if (self.eleron_now * 180.0/np.pi <= -15) :
                self.eleron_now = -15/180.0*np.pi

            self.t = self.t + dt
        return self.eleron_now


"""
control = Control()

ny = 1
a1 = control.get_elev_and_ny_new(ny)
print(a1)

gamma1 = 0.52
a2 = control.get_GammaAngle_and_eleron_now(gamma1)
print(a2)
"""
