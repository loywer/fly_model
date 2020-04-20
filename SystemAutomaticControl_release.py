import numpy as np
from model_function_upr import *
from matplotlib import pyplot as plt

#eps_ny = 0.01       # |(ny_now - ny_spec)| > eps_ny
dt = 0.02
k_elev = 0.01
k_eleron = 3.1

t = 0
T_elev = 20
T_eleron = 2.7 

class Control:

    def set_data (self, ny_spec, ny_now, gamma_spec, gamma_now, w0) :

        self.ny_spec = ny_spec
        self.ny_now = ny_now
        self.gamma_spec = gamma_spec
        self.gamma_now = gamma_now
        self.w0 = w0

    def get_data (self) :
        return self.elev_new, self.eleron_now


    def aperiodic_link (self, T) :
        dt = 0.02
        return (1 - np.exp(-dt/T))


    # Функция получения положения руля высоты и текущей перегрузки 
    #   ny_spec - заданное значение перегрузки 
    def get_elev_and_ny_new (self, ny_spec, ny_now) :
     
        #self.ny_now = get_ny(self.elev_new) / ACCELERATION_OF_GRAVITY # получение значения перегрузки в данный момент времени
        #self.elev_new = self.elev_new + k_elev * (ny_spec - ny_now) # получение отклонения рулей высоты
        #self.elev_new = self.elev_new + ((ny_spec - ny_now) * k_elev * (1 - np.exp(-dt/T_elev))) # получение отклонения рулей высоты
        self.transition_function = self.aperiodic_link(T_elev)
        self.elev_new = self.elev_new + ((ny_spec - ny_now) * k_elev * self.transition_function)
        #elev_new =  0.1 * (ny_spec - ny_now) # получение отклонения рулей высоты

        if (self.elev_new * 180.0/np.pi >= 26) :
             self.elev_new = 26/180.0*np.pi

        if (self.elev_new * 180.0/np.pi <= -28) :
            self.elev_new = -28/180.0*np.pi
        
        return self.elev_new

    # Функция получения положения элерона и крена самолета 
    #   gamma_spec - заданное значение крена самолета (желаемое)
    def get_GammaAngle_and_eleron_now (self, gamma_spec, gamma_now, w0) :
           
        #self.gamma_now, self.w0 = get_gamma(self.eleron_now)

        #self.eleron_now = k_eleron * (gamma_spec - gamma_now) - 3*w0
        #self.eleron_now = k_eleron * (gamma_spec - gamma_now) * (1 - np.exp(-dt/T_eleron)) - 3*w0

        self.transition_function = self.aperiodic_link(T_eleron)
        self.eleron_now = k_eleron * (gamma_spec - gamma_now) * self.transition_function(T_eleron) - 3*w0

        # t == dt

        if (self.eleron_now * 180.0/np.pi >= 20) :
            self.eleron_now = 20/180.0*np.pi

        if (self.eleron_now * 180.0/np.pi <= -15) :
            self.eleron_now = -15/180.0*np.pi

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
