from direct.showbase.ShowBase import ShowBase
from direct.showbase import DirectObject
from direct.task import Task
from direct.interval.IntervalGlobal import *
from direct.gui.OnscreenText import OnscreenText
from direct.gui.DirectGui import *
from panda3d.core import *
from numpy import append
from sys import exit



class MyApp(ShowBase, DirectObject.DirectObject):
    def __init__(self):
        ShowBase.__init__(self)
        # create buttons for menu
        self.menuLbl = DirectLabel (text = "MENU", pos = Vec3(0, 0, 0.9), scale = 0.1, textMayChange = 1)
        self.phiLbl = DirectLabel(text = "Enter latitude", pos = Vec3(0, 0, 0.8), scale = 0.08, textMayChange = 1)
        self.phiEnt = DirectEntry(scale = 0.04, pos = Vec3(-0.2, 0, 0.72))
        self.lambdaLbl = DirectLabel(text = "Enter longitude", pos = Vec3(0, 0, 0.6), scale = 0.08, textMayChange = 1)
        self.lambdaEnt = DirectEntry(scale = 0.04, pos = Vec3(-0.2, 0, 0.5))
        self.heightLbl = DirectLabel(text = "Enter height", pos = Vec3(0, 0, 0.4), scale = 0.08, textMayChange = 1)
        self.heightEnt = DirectEntry(scale = 0.04, pos = Vec3(-0.19, 0, 0.3))
        self.speedLbl = DirectLabel(text = "Enter speed", pos = Vec3(0, 0, 0.2), scale = 0.08, textMayChange = 1)
        self.speedEnt = DirectEntry(scale = 0.04, pos = Vec3(-0.19, 0, 0.1))
        self.rollLbl = DirectLabel(text = "Enter roll angle", pos = Vec3(0, 0, 0), scale = 0.08, textMayChange = 1)
        self.rollEnt = DirectEntry(scale = 0.04, pos = Vec3(-0.19, 0, -0.1))
        self.pitchLbl = DirectLabel(text = "Enter pitch angle", pos = Vec3(0, 0, -0.2), scale = 0.08, textMayChange = 1)
        self.pitchEnt = DirectEntry(scale = 0.04, pos = Vec3(-0.19, 0, -0.3))
        self.yawingLbl = DirectLabel(text = "Enter yawing angle", pos = Vec3(0, 0, -0.4), scale = 0.08, textMayChange = 1)
        self.yawingEnt = DirectEntry(scale = 0.04, pos = Vec3(-0.19, 0, -0.5))
        self.startBtn = DirectButton(text = "Start", scale = 0.1, command = self.setScene, pos = Vec3(0, 0, -0.7))
        self.points = []
        # binding keys
        self.accept("mouse1", self.set_coords)
        self.accept("escape", exit)
        self.accept("time-a-repeat", self.inc_roll)
        self.accept("time-d-repeat", self.dec_roll)
        self.accept("time-+-repeat", self.inc_speed)
        self.accept("time---repeat", self.dec_speed)
        self.accept("time-w-repeat", self.inc_overload)
        self.accept("time-s-repeat", self.dec_overload)
        self.overload = 0.0 # should be the result of the some function
   
    def setScene(self):
        self.acceptDlg = YesNoDialog(text = "Are you sure?", command = self.createScene)


    def createScene(self, clickedYes):
        if clickedYes:
            # hide menu elements
            self.acceptDlg.hide()
            self.menuLbl.hide()
            self.phiEnt.hide()
            self.phiLbl.hide()
            self.lambdaLbl.hide()
            self.lambdaEnt.hide()
            self.heightLbl.hide()
            self.heightEnt.hide()
            self.speedLbl.hide()
            self.speedEnt.hide()
            self.rollLbl.hide()
            self.rollEnt.hide()
            self.pitchLbl.hide()
            self.pitchEnt.hide()
            self.yawingLbl.hide()
            self.yawingEnt.hide()
            self.startBtn.hide()
            


            self.disableMouse()
            self.plane = loader.loadModel("/c/Panda3D-1.10.6-x64/models/boeing707.egg")
            self.plane.setScale(0.005, 0.005, 0.005)
            self.plane.setPos(0,0,0)
            self.cam.setPos(25.3, 2.26, 2.46)
            self.cam.lookAt(self.plane)

            self.plane.reparentTo(self.render)
            self.taskMgr.add(self.plane_coordiantes, "plane_coordiantes")
            
            # posInterval = time to move, finalPosition, startPosition
            posInterval1 = self.plane.posInterval(5, Point3(0, -6, -2), startPos=Point3(0,6,2))
            posInterval2 = self.plane.posInterval(5, Point3(0, 6, 2), startPos=Point3(0,-6,-2))
            self.get_Var() # read input from textboxes 
            self.planePace = Sequence(posInterval1, posInterval2, name = "planePace")
            self.planePace.loop()
        else:
            exit()



    def plane_coordiantes(self, task):
        cam_coords = []
        cam_coords.append(self.plane.getPos())
        return Task.cont





    def setNameLabel(self):
        # read input values
        self.phi = self.phiEnt.get()
        self.lambd = self.lambdaEnt.get()
        self.height = self.heightEnt.get()
        self.speed = self.speedEnt.get()
        self.roll = self.rollEnt.get()
        self.pitch = self.pitchEnt.get()
        self.yawing = self.yawingEnt.get()
        # show new values
        # TODO

    def get_Var(self):
        show = Func(self.setNameLabel)
        show.start()

    def set_coords(self):
        points = []
        if base.mouseWatcherNode.hasMouse():
            x = base.mouseWatcherNode.getMouseX()
            y = base.mouseWatcherNode.getMouseY()
            points = append(x,y)
        print(points)




    def inc_roll(self, when):
        self.roll = float(self.roll) + 0.5
        print(self.roll)
        return Task.cont


    def inc_speed(self, when):
        self.speed = float(self.speed) + 0.5
        print(self.speed)
        return Task.cont   

    def inc_overload(self, when):
        self.overload = float(self.overload) + 0.1
        print(self.overload)
        return Task.cont

    def dec_roll(self, when):
        self.roll = float(self.roll) - 0.5
        print(self.roll)
        return Task.cont
    
    def dec_speed(self, when):
        self.speed = float(self.speed) - 0.5
        print(self.speed)
        return Task.cont
    def dec_overload(self, when):
        self.overload = float(self.overload) - 0.1
        print(self.overload)
        return Task.cont



app = MyApp()
app.run()
