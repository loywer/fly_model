import numpy as np
from scipy.interpolate import RectBivariateSpline
from math import sin,cos
from matplotlib import pyplot as plt
#------------------------------------------------Начальные данные-------------------------------------------------------
#Геометрия   
Sw = 16.2     #м2 #Площадь крыла
b = 10.91184  #м  #Размах
c = 1.49352   #м  #Хорда

#Масса и инерция
mass = 2300*0.453592
inertia = (np.diag([948, 1346, 1967]))*1.35581795

#Аэродинамика    
alpha_data = np.deg2rad(np.array([-7.5,-5,-2.5,0,2.5,5,7.5,10,15,17,18,19.5])) #рад
deltaElev_data = np.deg2rad(np.array([-26, -20, -10, -5, 0, 7.5, 15, 22.5, 28])) #рад
deltaAile_data = np.deg2rad(np.array([-15, -10, -5, -2.5, 0, 5, 10, 15, 20])) #рад

#Входные данные для коэффициентов
Cx_data = np.array([0.044,0.034,0.03,0.03,0.036,0.048,0.067,0.093,0.15,0.169,0.177,0.184])                                                  #CD
Cx_DeltaElev_data = np.array([[0.0135,0.0119,0.0102,0.00846,0.0067,0.0049,0.00309,0.00117,-0.0033,-0.00541,-0.00656,-0.00838],              #CDDeltaElev
                    [0.0121,0.0106,0.00902,0.00738,0.00574,0.00406,0.00238,0.00059,-0.00358,-0.00555,-0.00661,-0.00831],
                    [0.00651,0.00552,0.00447,0.00338,0.00229,0.00117,0.0000517,-0.00114,-0.00391,-0.00522,-0.00593,-0.00706],
                    [0.00249,0.002,0.00147,0.000931,0.000384,-0.000174,-0.000735,-0.00133,-0.00272,-0.00337,-0.00373,-0.00429],
                    [0,0,0,0,0,0,0,0,0,0,0,0],
                    [-0.00089,-0.00015,0.00064,0.00146,0.00228,0.00311,0.00395,0.00485,0.00693,0.00791,0.00844,0.00929],
                    [0.00121,0.00261,0.00411,0.00566,0.00721,0.00879,0.0104,0.0121,0.016,0.0179,0.0189,0.0205],
                    [0.00174,0.00323,0.00483,0.00648,0.00814,0.00983,0.0115,0.0133,0.0175,0.0195,0.0206,0.0223],
                    [0.00273,0.00438,0.00614,0.00796,0.0098,0.0117,0.0135,0.0155,0.0202,0.0224,0.0236,0.0255]])

Cy_data = np.array([-0.571,-0.321,-0.083,0.148,0.392,0.65,0.918,1.195,1.659,1.789,1.84,1.889])                                              #CL
Cy_Q_data = np.array([7.282,7.282,7.282,7.282,7.282,7.282,7.282,7.282,7.282,7.282,7.282,7.282])                                             #CLQ
Cy_DdeltaElev_data = np.array([-0.132,-0.123,-0.082,-0.041,0,0.061,0.116,0.124,0.137])                                                      #CLDeltaElev

Cz_Betta_data = np.array([-0.268,-0.268,-0.268,-0.268,-0.268,-0.268,-0.268,-0.268,-0.268,-0.268,-0.268,-0.268])                             #CYBeta                
Cz_p_data = np.array([-0.032, -0.0372, -0.0418, -0.0463, -0.051, -0.0563, -0.0617, -0.068, -0.0783, -0.0812, -0.0824, -0.083])              #CYp
Cz_r_data = np.array([0.2018, 0.2054, 0.2087, 0.2115, 0.2139, 0.2159, 0.2175, 0.2187, 0.2198, 0.2198, 0.2196, 0.2194])                      #CYr
Cz_DeltaRud_data = (-1)*np.array([0.561, 0.561, 0.561, 0.561, 0.561, 0.561, 0.561, 0.561, 0.561, 0.561, 0.561, 0.561])                      #CYDEltaRud

Mx_betta_data = np.array([-0.178, -0.186, -0.1943, -0.202, -0.2103, -0.219, -0.2283, -0.2376, -0.2516, -0.255, -0.256, -0.257])             #Clbeta
Mx_p_data = np.array([-0.4968, -0.4678, -0.4489, -0.4595, 0.487, -0.5085, -0.5231, -0.4916, -0.301, -0.203, -0.1498, -0.0671])              #Clp
Mx_r_data = np.array([-0.09675, -0.05245, -0.01087, 0.02986, 0.07342, 0.1193, 0.1667, 0.2152, 0.2909, 0.3086, 0.3146, 0.3197])              #Clr
Mx_DeltaRud_data = (-1)*np.array([0.091, 0.082, 0.072, 0.063, 0.053, 0.0432, 0.0333, 0.0233, 0.0033, -0.005, -0.009, -0.015])               #ClDeltaRud
Mx_DeltaAile_data = np.array([-0.078052, -0.059926, -0.036422, -0.018211, 0, 0.018211, 0.036422, 0.059926, 0.078052])                       #ClDeltaAile

My_betta_data=np.array([0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126])                    #CNbeta
My_p_data=np.array([0.03, 0.016, 0.00262, -0.0108, -0.0245, -0.0385, -0.0528, -0.0708, -0.113, -0.1284, -0.1356, -0.1422])                 #CNp
My_r_data=np.array([-0.028, -0.027, -0.027, -0.0275, -0.0293, -0.0325, -0.037, -0.043, -0.05484, -0.058, -0.0592, -0.06015])               #CNr
My_DeltaRud_data=(-1)*np.array([-0.211, -0.215, -0.218, -0.22, -0.224, -0.226, -0.228, -0.229, -0.23, -0.23, -0.23, -0.23])                #CNDeltaRud
My_DeltaAile_data=np.array([[-0.004321, -0.002238, -0.0002783, 0.001645, 0.003699, 0.005861, 0.008099, 0.01038, 0.01397, 0.01483, 0.01512, 0.01539],      #CNDeltaAile
                [-0.003318, -0.001718, -0.0002137, 0.001263, 0.00284, 0.0045, 0.006218, 0.00797, 0.01072, 0.01138, 0.01161, 0.01181],
                [-0.002016, -0.001044, -0.000123, 0.0007675, 0.00173, 0.002735, 0.0038, 0.004844, 0.00652, 0.00692, 0.00706, 0.0072],
                [-0.00101, -0.000522, -0.0000649, 0.000384, 0.000863, 0.00137, 0.0019, 0.00242, 0.00326, 0.00346, 0.00353, 0.0036],
                [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                [0.00101, 0.00052, 0.000065, -0.000384, -0.00086, -0.0014, -0.002, -0.002422, -0.00326, -0.00346, -0.00353, -0.0036],
                [0.00202, 0.001044, 0.00013, -0.0008, -0.00173, -0.002735, -0.0038, -0.004844, -0.00652, -0.00692, -0.00706, -0.0072],
                [0.00332, 0.00172, 0.000214, -0.001263, -0.00284, -0.0045, -0.00622, -0.008, -0.01072, -0.01138, -0.01161, -0.01181],
                [0.004321, 0.00224, 0.00028, -0.001645, -0.0037, -0.00586, -0.0081, -0.0104, -0.014, -0.01483, -0.01512, -0.0154]])

Mz_data=np.array([0.0597, 0.0498, 0.0314, 0.0075, -0.0248, -0.068, -0.1227, -0.1927, -0.3779, -0.4605, -0.5043, -0.5496])                  #CM
Mz_q_data=np.array([-6.232, -6.232, -6.232, -6.232, -6.232, -6.232, -6.232, -6.232, -6.232, -6.232, -6.232, -6.232])                       #CMq
Mz_DeltaElev_data=np.array([0.3302, 0.3065, 0.2014, 0.1007, -0.0002, -0.1511, -0.2863, -0.3109, -0.345])                                    #CMDeltaElev

class Aerodynamics:

    def __init__ (self):

        self.w = [0,0,0]
        self.angl = [0,0,0]    #Крен, рыскание,тангаж (Последовательность)
        self.V = [50,0,0]
        H = 500
        self.koordinat = [0,H,0]   # X,Y,Z
        self.dt = 0.02
        self.ro=1.25    #Плотность
        self.DeltaElev =0
        self.DeltaRud=0
        self.DeltaAile = 0
        self.P=0
        self.g=9.81
        self.G=[0,-self.g,0]

    def Set_data (self,P,dt,DeltaElev,DeltaRud,DeltaAile,ro,g):

        self.dt = dt
        self.ro=ro    #Плотность
        self.DeltaElev =DeltaElev 
        self.DeltaRud=DeltaRud
        self.DeltaAile = DeltaAile
        self.P=P
        self.g=g
        self.G=[0,-g,0]

    def speed_abs (self):

        Modul_V = np.sqrt((self.V[0]**2)+(self.V[1]**2)+(self.V[2]**2))

        return Modul_V

# Вычисление углов атаки и скольжения
    def finding_angles_of_attack_and_slip (self):
        
        Modul_V = self.speed_abs()

        if self.V[0] >= 0:
            alpha = (-1)*np.arcsin((self.V[1])/(np.sqrt((self.V[0]**2)+(self.V[1]**2))))
        else:
            alpha = (-1)*np.pi + np.arcsin((abs(self.V[1]))/(np.sqrt(self.V[0]**2+self.V[1]**2))*np.sign(self.V[1]))
        betta = np.arcsin(self.V[2]/Modul_V)

        return alpha,betta

    #Матрица перехода из НСК в ССК
    def matrix_NSK_CCK (self):

        rotation_matrix_around_OY_axis = np.array([[cos(self.angl[1]),0,(-1)*sin(self.angl[1])],              #Поворот вокруг оси OY на угол рыскания (phi)
            [0,1,0],
            [sin(self.angl[1]),0,cos(self.angl[1])]])
        rotation_matrix_around_OZ_axis = np.array([[cos(self.angl[2]),sin(self.angl[2]),0],                   #Поворот вокруг оси OZ на угол тангажа (v)
            [(-1)*sin(self.angl[2]),cos(self.angl[2]),0],
            [0,0,1]])
        rotation_matrix_around_OX_axis = np.array([[1,0,0],                                                             #Поворот вокруг оси OX на угол крена (gamma)
            [0,cos(self.angl[0]),sin(self.angl[0])],
            [0,(-1)*sin(self.angl[0]),cos(self.angl[0])]])
        transition_matrix_NSK_CCK = (rotation_matrix_around_OX_axis.dot(rotation_matrix_around_OZ_axis)).dot(rotation_matrix_around_OY_axis)

        return transition_matrix_NSK_CCK

#Матрица перехода из СкСК в ССК
    def matrix_CkCK_CCK (self):

        alpha,betta = self.finding_angles_of_attack_and_slip ()

        transition_matrix_CkCK_CCK = np.array([[(cos(alpha) * cos(betta)) , (sin(alpha)) , ((-1.0) * cos(alpha) * sin(betta))],
                    [((-1.0) * sin(alpha) * cos(betta)) , (cos(alpha)) , (sin(alpha) * sin(betta))],
                    [(sin(betta)) , 0 , (cos(betta))]])

        return transition_matrix_CkCK_CCK


    def Aero_Forces_coefficient (self,Modul_V,alpha,betta):

    #Вычисление Cx
        CD_delta_elev_interp = RectBivariateSpline(deltaElev_data,alpha_data,Cx_DeltaElev_data)
        cx_1 = CD_delta_elev_interp(self.DeltaElev,alpha)[0, 0]
        cx_2 = np.interp(alpha,alpha_data,Cx_data)
        cx = cx_1 + cx_2
    #Вычисление Cy
        cy_1 = np.interp(alpha,alpha_data,Cy_data)
        cy_2 = np.interp(alpha,alpha_data,Cy_Q_data)
        cy_3 = np.interp(self.DeltaElev,deltaElev_data,Cy_DdeltaElev_data)
        cy = cy_1+(c/(2*Modul_V))*cy_2*self.w[2] + cy_3                   
    #Вычисление Cz
        cz_1 = np.interp(alpha,alpha_data, Cz_Betta_data) 
        cz_2 = np.interp(alpha,alpha_data,Cz_p_data) 
        cz_3 = np.interp(alpha,alpha_data,Cz_r_data)
        cz_4 = np.interp(alpha,alpha_data,Cz_DeltaRud_data)
        cz = cz_1*betta+(b/(2*Modul_V))*(cz_2*self.w[0]+cz_3*self.w[1]) + cz_4*self.DeltaRud

        Aero_Forces_coefficient = np.array([cx,cy,cz])

        return Aero_Forces_coefficient

    def Aero_Forces(self):

        Modul_V = self.speed_abs()
        alpha,betta = self.finding_angles_of_attack_and_slip ()
        Aero_Forces_coefficient = self.Aero_Forces_coefficient(Modul_V,alpha,betta)

        Fx=(-1)*Aero_Forces_coefficient[0] * self.ro * (Modul_V**2) * (Sw / 2)
        Fy=Aero_Forces_coefficient[1] * self.ro * (Modul_V**2) * (Sw / 2)
        Fz=Aero_Forces_coefficient[2] * self.ro * (Modul_V**2) * (Sw / 2)

        Forces=np.array([Fx,Fy,Fz])

        return Forces

    def Aero_Moment_coefficient (self,Modul_V,alpha,betta):

    #Вычисление mx
        mx_1 = np.interp(alpha, alpha_data,Mx_betta_data)
        mx_2 = np.interp(alpha,alpha_data,Mx_p_data)
        mx_3 = np.interp(alpha,alpha_data,Mx_r_data)
        mx_4 = np.interp(alpha,alpha_data,Mx_DeltaRud_data)
        mx_5 = np.interp(self.DeltaAile,deltaAile_data,Mx_DeltaAile_data)
        mx = 0.1*mx_1 * betta + (b/(2*Modul_V))*(mx_2*self.w[0]+mx_3*self.w[1])+0.075*mx_4*self.DeltaRud+mx_5    
    #Вычисление my
        my_1 = np.interp(alpha,alpha_data,My_betta_data)
        my_2 = np.interp(alpha,alpha_data,My_p_data)
        my_3 = np.interp(alpha,alpha_data,My_r_data)
        my_4 = np.interp(alpha,alpha_data,My_DeltaRud_data)
        CNDeltaAile_interp = RectBivariateSpline(deltaAile_data,alpha_data,My_DeltaAile_data)
        my_5 = CNDeltaAile_interp(self.DeltaAile,alpha)[0, 0]
        my = my_1*betta+(b/(2*Modul_V))*(my_2*self.w[0]+my_3*self.w[1])+0.075*my_4*self.DeltaRud + my_5 
    #Вычисление mz
        mz_1 = np.interp(alpha,alpha_data,Mz_data)
        mz_2 = np.interp(alpha,alpha_data,Mz_q_data)
        mz_3 = np.interp(self.DeltaElev,deltaElev_data,Mz_DeltaElev_data)
        mz = mz_1+(c/(2*Modul_V))*2*mz_2*self.w[2]+mz_3  
         
        Aero_Moment_coefficient = np.array([mx,my,mz])

        return Aero_Moment_coefficient

    def Aero_Moment (self):

        Modul_V = self.speed_abs()
        alpha,betta = self.finding_angles_of_attack_and_slip ()
        Aero_Moment_coefficient = self.Aero_Moment_coefficient(Modul_V,alpha,betta)

        Mx = Aero_Moment_coefficient[0] * self.ro * b *(Modul_V**2) * (Sw / 2)
        My =Aero_Moment_coefficient[1] * self.ro * b *(Modul_V**2) * (Sw / 2)
        Mz = Aero_Moment_coefficient[2] * self.ro * c *(Modul_V**2) * (Sw / 2)

        Moment = np.array([Mx,My,Mz])

        return Moment

    #Перевод сил и ускорения в ССК
    def Translation_Forces_and_G_to_CCK (self):

        transition_matrix_NSK_CCK = self.matrix_NSK_CCK ()
        transition_matrix_CkCK_CCK = self.matrix_CkCK_CCK ()
        Forces = self.Aero_Forces ()

        Forces_CCK = np.dot(transition_matrix_CkCK_CCK,Forces)
        G_CCK = np.dot(transition_matrix_NSK_CCK,self.G)  
        return Forces_CCK,G_CCK

#Функция вычисляющая силы действующие на самолет
    def Forces_all  (self):

        Forces_CCK,G_CCK = self.Translation_Forces_and_G_to_CCK ()

        Forces_all = ((Forces_CCK + self.P) / mass) + G_CCK

        return Forces_all

    def overload (self):

        Forces_CCK = self.Translation_Forces_and_G_to_CCK ()[0]

        n =((Forces_CCK + self.P) / (mass*self.g))

        return n

# Вычисление линейных и угловых ускорений
    def finding_accelerations (self):

        Forces_all = self.Forces_all ()
        Moment = self.Aero_Moment ()

        V_a = (Forces_all - np.cross(self.w,self.V))                # Линейное ускорение V_a                           

        J_w = np.dot(self.w,inertia)
        w_J_W = np.cross(self.w,J_w)
        M_w_J_w = Moment - w_J_W
        inertia_obr = np.linalg.inv(inertia)

        w_a = np.dot(inertia_obr,M_w_J_w)                            # Угловое ускорение w_a

        return V_a,w_a

    def right_side_of_the_DU_angl (self):

        pitchin = self.w[1]*sin(self.angl[0]) + self.w[2]*cos(self.angl[0])
        rolling = self.w[0] - (self.w[1]*cos(self.angl[0]) - self.w[2]*sin(self.angl[0]))*np.tan(self.angl[2])
        yaming = (self.w[1]*cos(self.angl[0]) - self.w[2]*sin(self.angl[0]))/cos(self.angl[0])
        angl_prav = np.array([rolling,yaming,pitchin])

        return angl_prav

    def right_side_of_the_DU_koord (self):

        transition_matrix_NSK_CCK = self.matrix_NSK_CCK ()

        koordinat_prav = np.dot(transition_matrix_NSK_CCK.T,self.V)

        return koordinat_prav

#Интегратор
    def Integrator (self):  

        right_side_of_the_DU_V_and_w = self.finding_accelerations ()
        right_side_of_the_DU_angl = self.right_side_of_the_DU_angl ()
        right_side_of_the_DU_koord = self.right_side_of_the_DU_koord ()

        V_new = self.V + right_side_of_the_DU_V_and_w[0]*self.dt
        w_new = self.w + right_side_of_the_DU_V_and_w[1]*self.dt
        angl_new = self.angl + right_side_of_the_DU_angl*self.dt
        koordinat_new = self.koordinat + right_side_of_the_DU_koord*self.dt

        return V_new,w_new,angl_new,koordinat_new

#Функция вызова Get
    def Get_data (self):

        n = self.overload()
        alpha,betta = self.finding_angles_of_attack_and_slip()
        V_a,w_a = self.finding_accelerations()
        self.V,self.w,self.angl,self.koordinat = self.Integrator()

        return n,alpha,betta,V_a,w_a,self.V,self.w,self.angl,self.koordinat
