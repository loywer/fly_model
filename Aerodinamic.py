import numpy as np
from numpy import linalg as LA
class Aerodynamics () :
    # Входные данные
    def __init__(self, V, alpha, betta, DeltaElev, DeltaRud, DeltaAile, Alphadot,dt):
    #Масса и инерция
        self.mass = 1043.2
        self.inertia = np.diag([948, 1346, 1967])
        #_______________________Комменты убрать
        #self.Ix=self.inertia[0,0]
        #self.Iy=self.inertia[1,1]
        #self.Iz=self.inertia[2,2]
        self.ang_vel=np.array([0,0,0])
        self.q=self.ang_vel[0,0]
        self.p=self.ang_vel[0,1]
        self.r=self.ang_vel[0,2]
    #Аэродинамика    
        self.alpha_1 =np.rad2deg(np.array([-7.5,-5,-2.5,0,2.5,5,7.5,10,15,17,18,19.5])) #рад
        self.betta_1 = np.rad2deg(np.array([0,1,2,3,4,5,6,7,8,9,10,11]))  #рад
        self.deltaElev = np.rad2deg(np.array([-26, -20, -10, -5, 0, 7.5, 15, 22.5, 28]))  #рад
        self.deltaAile = np.rad2deg(np.array([-15, -10, -5, -2.5, 0, 5, 10, 15, 20])) #рад
        #Входные данные для Cx
        self.CD = np.array([0.044,0.034,0.03,0.03,0.036,0.048,0.067,0.093,0.15,0.169,0.177,0.184])
        self.CDDeltaElev = np.array([[0.0135,0.0119,0.0102,0.00846,0.0067,0.0049,0.00309,0.00117,-0.0033,-0.00541,-0.00656,-0.00838],
                            [0.0121,0.0106,0.00902,0.00738,0.00574,0.00406,0.00238,0.00059,-0.00358,-0.00555,-0.00661,-0.00831],
                            [0.00651,0.00552,0.00447,0.00338,0.00229,0.00117,0.0000517,-0.00114,-0.00391,-0.00522,-0.00593,-0.00706],
                            [0.00249,0.002,0.00147,0.000931,0.000384,-0.000174,-0.000735,-0.00133,-0.00272,-0.00337,-0.00373,-0.00429],
                            [0,0,0,0,0,0,0,0,0,0,0,0],
                            [-0.00089,-0.00015,0.00064,0.00146,0.00228,0.00311,0.00395,0.00485,0.00693,0.00791,0.00844,0.00929],
                            [0.00121,0.00261,0.00411,0.00566,0.00721,0.00879,0.0104,0.0121,0.016,0.0179,0.0189,0.0205],
                            [0.00174,0.00323,0.00483,0.00648,0.00814,0.00983,0.0115,0.0133,0.0175,0.0195,0.0206,0.0223],
                            [0.00273,0.00438,0.00614,0.00796,0.0098,0.0117,0.0135,0.0155,0.0202,0.0224,0.0236,0.0255]])
        #Входные данные для Cy
        self.CL = np.array([-0.571,-0.321,-0.083,0.148,0.392,0.65,0.918,1.195,1.659,1.789,1.84,1.889])
        self.CLAlphadot = np.array([2.434,2.362,2.253,2.209,2.178,2.149,2.069,1.855,1.185,0.8333,0.6394,0.4971])
        self.CLQ = np.array([7.282,7.282,7.282,7.282,7.282,7.282,7.282,7.282,7.282,7.282,7.282,7.282])
        self.CLDdeltaElev = np.array([-0.132,-0.123,-0.082,-0.041,0,0.061,0.116,0.124,0.137])   
        #Входные данные для Cz
        self.CYBeta = np.array([-0.268,-0.268,-0.268,-0.268,-0.268,-0.268,-0.268,-0.268,-0.268,-0.268,-0.268,-0.268])
        self.CYp = np.array([-0.032, -0.0372, -0.0418, -0.0463, -0.051, -0.0563, -0.0617, -0.068, -0.0783, -0.0812, -0.0824, -0.083])
        self.CYr = np.array([0.2018, 0.2054, 0.2087, 0.2115, 0.2139, 0.2159, 0.2175, 0.2187, 0.2198, 0.2198, 0.2196, 0.2194])
        self.CYDeltaRud = (-1)*np.array([0.561, 0.561, 0.561, 0.561, 0.561, 0.561, 0.561, 0.561, 0.561, 0.561, 0.561, 0.561])
        #Входные данные для mx
        self.Clbeta=np.array([-0.178, -0.186, -0.1943, -0.202, -0.2103, -0.219, -0.2283, -0.2376, -0.2516, -0.255, -0.256, -0.257])
        self.Clp=np.array([-0.4968, -0.4678, -0.4489, -0.4595, 0.487, -0.5085, -0.5231, -0.4916, -0.301, -0.203, -0.1498, -0.0671])
        self.Clr=np.array([-0.09675, -0.05245, -0.01087, 0.02986, 0.07342, 0.1193, 0.1667, 0.2152, 0.2909, 0.3086, 0.3146, 0.3197])
        self.ClDeltaRud=(-1)*np.array([0.091, 0.082, 0.072, 0.063, 0.053, 0.0432, 0.0333, 0.0233, 0.0033, -0.005, -0.009, -0.015])
        self.ClDeltaAile=np.array([-0.078052, -0.059926, -0.036422, -0.018211, 0, 0.018211, 0.036422, 0.059926, 0.078052])
        #Входные данные для my
        self.CM=np.array([0.0597, 0.0498, 0.0314, 0.0075, -0.0248, -0.068, -0.1227, -0.1927, -0.3779, -0.4605, -0.5043, -0.5496])
        self.CMq=np.array([-6.232, -6.232, -6.232, -6.232, -6.232, -6.232, -6.232, -6.232, -6.232, -6.232, -6.232, -6.232])
        self.CMAlphadot=np.array([-6.64, -6.441, -6.146, -6.025, -5.942, -5.861, -5.644, -5.059, -3.233, -2.273, -1.744, -1.356])
        self.CMDeltaElev=np.array([0.3302, 0.3065, 0.2014, 0.1007, -0.0002, -0.1511, -0.2863, -0.3109, -0.345])
        #Входные данные для mz
        self.CNbeta=np.array([0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126, 0.0126])
        self.CNp=np.array([0.03, 0.016, 0.00262, -0.0108, -0.0245, -0.0385, -0.0528, -0.0708, -0.113, -0.1284, -0.1356, -0.1422])
        self.CNr=np.array([-0.028, -0.027, -0.027, -0.0275, -0.0293, -0.0325, -0.037, -0.043, -0.05484, -0.058, -0.0592, -0.06015])
        self.CNDeltaRud=(-1)*np.array([-0.211, -0.215, -0.218, -0.22, -0.224, -0.226, -0.228, -0.229, -0.23, -0.23, -0.23, -0.23])
        self.CNDeltaAile=np.array([[-0.004321, -0.002238, -0.0002783, 0.001645, 0.003699, 0.005861, 0.008099, 0.01038, 0.01397, 0.01483, 0.01512, 0.01539],
                        [-0.003318, -0.001718, -0.0002137, 0.001263, 0.00284, 0.0045, 0.006218, 0.00797, 0.01072, 0.01138, 0.01161, 0.01181],
                        [-0.002016, -0.001044, -0.000123, 0.0007675, 0.00173, 0.002735, 0.0038, 0.004844, 0.00652, 0.00692, 0.00706, 0.0072],
                        [-0.00101, -0.000522, -0.0000649, 0.000384, 0.000863, 0.00137, 0.0019, 0.00242, 0.00326, 0.00346, 0.00353, 0.0036],
                        [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],
                        [0.00101, 0.00052, 0.000065, -0.000384, -0.00086, -0.0014, -0.002, -0.002422, -0.00326, -0.00346, -0.00353, -0.0036],
                        [0.00202, 0.001044, 0.00013, -0.0008, -0.00173, -0.002735, -0.0038, -0.004844, -0.00652, -0.00692, -0.00706, -0.0072],
                        [0.00332, 0.00172, 0.000214, -0.001263, -0.00284, -0.0045, -0.00622, -0.008, -0.01072, -0.01138, -0.01161, -0.01181],
                        [0.004321, 0.00224, 0.00028, -0.001645, -0.0037, -0.00586, -0.0081, -0.0104, -0.014, -0.01483, -0.01512, -0.0154]])
    #Геометрия   
        self.Sw = 16.2     #м2 #Площадь крыла
        self.b = 10.91184  #м  #Размах
        self.c = 1.49352   #м  #Хорда
    #Окружающая среда
        self.ro=1.225      #Плотность
    #Данные полученные из других программ
        self.V=V
        self.alpha = alpha
        self.betta = betta  
        self.DeltaElev = DeltaElev
        self.DeltaRud = DeltaRud
        self.DeltaAile = DeltaAile
        self.Alphadot = Alphadot
        self.dt=dt

# Функция вычисляющая коэффициенты cx,cy,cz
    def Aero_Forces_coeff (self):
    #Вычисление Cx
        S=[]
        for i in range(12):
            CDDD = np.array(np.interp( self.DeltaElev, self.deltaElev, self.CDDeltaElev[:,i]))
            S.append( CDDD)
        S = np.asarray(S)      
        CDD = np.interp(self.alpha, self.alpha_1, S)
        cx_1 = np.interp(self.alpha, self.alpha_1, self.CD)
        self.cx = CDD + cx_1
    #Вычисление Cy
        cy_1 = np.interp(self.alpha,self.alpha_1,self.CL)
        cy_2 = np.interp(self.alpha,self.alpha_1,self.CLAlphadot)
        cy_3 = np.interp(self.alpha,self.alpha_1,self.CLQ)
        cy_4 = np.interp(self.DeltaElev,self.deltaElev,self.CLDdeltaElev)
        self.cy = cy_1 +((self.c)/2*self.V)*(cy_2*self.Alphadot+cy_3 * self.q) + cy_4 * self.DeltaElev
    #Вычисление Cz
        cz_1 = np.interp(self.betta,self.betta_1,self.CYBeta) 
        cz_2 = np.interp(self.alpha,self.alpha_1,self.CYp) 
        cz_3 = np.interp(self.alpha,self.alpha_1,self.CYr)
        cz_4 = np.interp(self.alpha,self.alpha_1,self.CYDeltaRud)
        self.cz = (cz_1 * self.betta)+((self.b)/2*self.V)*(cz_2 * self.p+cz_3 * self.r)+ +(cz_4*self.DeltaRud)
        koeff = np.array([self.cx,self.cy,self.cz])
        return koeff

#Вычисление аэродинамических сил
    def Aero_Forces(self,koeff):
        Fx=koeff[0,0] * self.ro * (self.V**2) * (self.Sw)/2
        Fy=koeff[0,1] * self.ro * (self.V**2) * (self.Sw)/2
        Fz=koeff[0,2] * self.ro * (self.V**2) * (self.Sw)/2
        F=np.array([Fx,Fy,Fz])
        return F

#Вычисление коэффициентов аэродинамических моментов 
    def Aero_Moment_coeff (self):
    #Вычисление mx
        mx_1 = np.interp(self.alpha, self.alpha_1, self.Clbeta)
        mx_2 = np.interp(self.alpha, self.alpha_1, self.Clp)
        mx_3 = np.interp(self.alpha, self.alpha_1, self.Clr)
        mx_4 = np.interp(self.alpha, self.alpha_1, self.ClDeltaRud)
        mx_5 = np.interp(self.DeltaAile, self.deltaAile, self.ClDeltaAile)
        self.mx = (mx_1 * self.betta) + ((self.c) / 2 * self.V) * (mx_2 * self.p + mx_3 * self.r)+(mx_4 * self.DeltaRud)+(mx_5 * self.DeltaAile)
    #Вычисление my
        my_1 = np.interp(self.alpha, self.alpha_1, self.CM)
        my_2 = np.interp(self.alpha, self.alpha_1, self.CMq)
        my_3 = np.interp(self.alpha, self.alpha_1, self.CMAlphadot)
        my_4 = np.interp(self.DeltaElev, self.deltaElev, self.CMDeltaElev)
        self.my = (my_1 * self.alpha) +((self.c) / 2 * self.V) * (my_2 * self.q+my_3 * self.Alphadot) + my_4 * self.DeltaElev
    #Вычисление mz
        mz_1 = np.interp(self.betta, self.betta_1, self.CNbeta)
        mz_2 = np.interp(self.alpha, self.alpha_1, self.CNp)
        mz_3 = np.interp(self.alpha, self.alpha_1, self.CNr)
        mz_4 = np.interp(self.alpha, self.alpha_1, self.CNDeltaRud)
        k = []
        for i in range(12):
            CNNN = np.array(np.interp(self.DeltaAile, self.deltaAile, self.CNDeltaAile[:,i]))
            k.append(CNNN)
        k = np.asarray(k)     
        mz_5 = np.interp(self.alpha, self.alpha_1, k)
        self.mz = (mz_1 * self.betta) + ((self.c) / 2 * self.V) * (mz_2 * self.p + mz_3 * self.r) + (mz_4 * self.DeltaRud) + mz_5 * self.DeltaAile
        moment = np.array([self.mx,self.my,self.mz])
        return moment

#Вычисление аэродинамических моментов
    def Aero_Moment (self,moment):
        Mx = moment[0,0] * self.ro * self.b *(self.V**2) * self.Sw/2
        My = moment[0,1] * self.ro * self.c *(self.V**2) * self.Sw/2
        Mz = moment[0,2] * self.ro * self.b *(self.V**2) * self.Sw/2
        M = np.array([[Mx],[My],[Mz]])
        return M

    def Integrator (self,F,dt):
        V_new=self.V+F*self.dt
        return V_new
        
# Вычисление угловых скоростей
 #   def Angular_vel(self):
        #self.w=np.array([r],[q],[w])
        #self.k = self.inertia*self.w
        #dw_dt=LA.inv(self.inertia)*(self.M-np.dot(self.w,self.k)
