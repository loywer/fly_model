import numpy as np
from model_function_upr import *
from matplotlib import pyplot as plt

array_elev = []     # массив значений отклонения руля (угол отклонения руля [град])
array_ny_now = []   # массив значений перегрузок в данный момент времени
array_ny_spec = []  # массив значений заданной перегрузки (постоянная)

array_t = []        # массив значений времени для рулей высоты

array_time_eleron = []  # массив значений времени для элеронов
array_eleron_now = []   # массив значений положения элерона
array_gamma_now = []    # массив значений крена самолета
array_gamma_spec = []   # заданный крен самолета

eps_ny = 0.01       # |(ny_now - ny_spec)| > eps_ny
dt = 0.02


def aperiodic_link (T) :
    dt = 0.02
    return (1 - np.exp(-dt/T))

# Функция получения положения руля высоты и текущей перегрузки 
#   ny_spec - заданное значение перегрузки   
def get_elev_and_ny_new (ny_spec) :
    
    global array_elev
    global array_ny_now
    global array_ny_spec
        
    array_elev.clear()
    array_ny_now.clear()
    array_ny_spec.clear()
    
    array_t.clear()
     
    ny_now = 0 # начальное значение перегрузки
    elev_new = 0 # начальное положение угла руля высоты
     
    k = 3
    t = 0
    T = 3.5

    #while ((abs(ny_now - ny_spec) > eps_ny) and (t <= 100)) :
    while (t <= 20) :
    #while (ny_now != ny_spec) :
    
        array_t.append(t)
        array_elev.append(elev_new*180.0/np.pi)
        array_ny_now.append(ny_now)
        array_ny_spec.append(ny_spec)

        ny_now = get_ny(elev_new) / ACCELERATION_OF_GRAVITY # получение значения перегрузки в данный момент времени
        #elev_new = get_ny_2(Cy_now) / ACCELERATION_OF_GRAVITY # получение значения перегрузки в данный момент времени
        
        #elev_new =  0.1 * (ny_spec - ny_now) # получение отклонения рулей высоты
        #elev_new = elev_new + k * (ny_spec - ny_now) # получение отклонения рулей высоты
        #elev_new = elev_new + ((ny_spec - ny_now) * k * (1 - np.exp(-t/T))) # получение отклонения рулей высоты
        transition_function = aperiodic_link(T)
        elev_new = elev_new + ((ny_spec - ny_now) * k * transition_function)

        if (elev_new * 180.0/np.pi >= 26) :
            elev_new = 26/180.0*np.pi

        if (elev_new * 180.0/np.pi <= -28) :
            elev_new = -28/180.0*np.pi
    
        t = t + dt

    return elev_new, array_elev

# Функция получения положения элерона и крена самолета 
#   gamma_spec - заданное значение крена самолета (желаемое)
def get_GammaAngle_and_eleron_now (gamma_spec) :
    
    global array_time_eleron
    
    global array_eleron_now
    global array_gamma_now
    global array_gamma_spec

    array_eleron_now.clear()
    array_gamma_now.clear()
    array_gamma_spec.clear()
    array_time_eleron.clear()

    eleron_now = 0.0 # начальное положение элерона
    gamma_now = 0 # начальный угол крена самолета
    
    k = 250
    t = 0
    dt = 0.02
    T = 2.5

    #while ((abs(gamma_spec - gamma_now) > 0.01) and (t <= 20)) :
    while  (t <= 10):
        
        array_time_eleron.append(t)
        array_eleron_now.append(eleron_now*180.0/np.pi)
        array_gamma_now.append(gamma_now*180.0/np.pi)
        array_gamma_spec.append(gamma_spec*180.0/np.pi)

        #print ("gamma = ", gamma)

        gamma_now, w0 = get_gamma(eleron_now)
        #print ("gamma_now = ", gamma_now)

        #eleron_now = k*(gamma_spec - gamma_now) - 3*w0
        #eleron_now = k*(gamma_spec - gamma_now) * (1 - np.exp(-t/T)) - 3*w0
        transition_function = aperiodic_link(T)
        eleron_now =  k*(gamma_spec - gamma_now) * transition_function - 3*w0 # получение отклонения рулей высоты
        #print ("eleron_now = ", eleron_now)

        if (eleron_now * 180.0/np.pi >= 20) :
            eleron_now = 20/180.0*np.pi

        if (eleron_now * 180.0/np.pi <= -15) :
            eleron_now = -15/180.0*np.pi

        t = t + dt
    return eleron_now, gamma_now

# Функция записи в файл полученных результатов ny и elev
#   filename - имя файла для записи
#   ny_spec - заданное значение перегрузки 
def WriteToFile_ny_elev (filename, ny_spec) :
   
    get_elev_and_ny_new (ny_spec)

    np.savetxt(filename, 
        np.column_stack((array_t, array_elev, array_ny_now)), 
        fmt = '%5.2f | %8.5f | %8.5f', 
        header = ('Для перегрузки = ' + str(ny_spec) + '\nt   | elev_new |  ny_now\n-------------------------'), 
        delimiter = ' | ')
    
    return

# Функция записи в файл полученных результатов eleron и gamma
#   filename - имя файла для записи
#   gamma_spec - желаеиый угол крена 
def WriteToFile_eleron_gamma (filename, gamma_spec) :
   
    get_GammaAngle_and_eleron_now (gamma_spec)
   
    np.savetxt(filename, 
        np.column_stack((array_time_eleron, array_eleron_now, array_gamma_now)), 
        fmt = '%5.2f | %10.5f | %8.5f', 
        header = ('Желаемый угол крена = ' + str(gamma_spec*180.0/np.pi) + '\nt   | eleron_now | gamma_now\n-----------------------------'), 
        delimiter = ' | ')
    
    return

fig1 = plt.figure()
plt.title("Изменение отклонения руля высоты")
plt.xlabel ("Время наблюдения, [c]")
plt.ylabel ("Угол текущего положения рулей, [град]")

#for i in range(1, 6, 2): # цикл от 1 до 5 с шагом 2, 
    #i - значения перегрузки
get_elev_and_ny_new (0.5)
plt.plot (array_t, array_elev, label = 'Заданная перегрузка = ' + str(0.5))
    #print("\ni = ", i)

plt.legend()
plt.grid()


fig2 = plt.figure()
get_elev_and_ny_new (0.5)
plt.title("Изменение значения перегрузки")
plt.xlabel ("Время наблюдения, [c]")
plt.ylabel ("Перегрузка, []")
plt.plot (array_t, array_ny_now, label = 'Текущая перегрузка')
plt.plot (array_t, array_ny_spec, label = 'Заданная перегрузка')
plt.legend()
plt.grid()


fig3 = plt.figure()
get_GammaAngle_and_eleron_now (0.52)
plt.title("ЭЛЕРОНЫ")
plt.xlabel ("Время наблюдения, [c]")
plt.ylabel ("Угол отклонения элеронов, [град]")
plt.plot (array_time_eleron, array_eleron_now)
plt.grid()


fig4 = plt.figure()
#get_GammaAngle_and_eleron_now (0.52)
plt.title("КРЕН")
plt.xlabel ("Время наблюдения, [c]")
plt.ylabel ("Крен самолета, [град]")
plt.plot (array_time_eleron, array_gamma_now, label = 'Текущий крен самолета')
plt.plot (array_time_eleron, array_gamma_spec, label = 'Заданный крен самолета')
plt.legend()
plt.grid()

plt.show()

WriteToFile_ny_elev ("output_ny_elev1.txt", 1)
WriteToFile_ny_elev ("output_ny_elev3.txt", 3)
WriteToFile_ny_elev ("output_ny_elev5.txt", 5)
