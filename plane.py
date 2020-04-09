from direct.showbase.ShowBase import ShowBase
from direct.showbase import DirectObject
from direct.task import Task
from panda3d.core import LineSegs, NodePath, Point3D
from direct.interval.IntervalGlobal import Sequence
from panda3d.core import Point3
from numpy import append
from sys import exit



class MyApp(ShowBase):
    def __init__(self):
        ShowBase.__init__(self)
        self.plane = loader.loadModel("/c/Panda3D-1.10.6-x64/models/boeing707.egg")
        self.plane.setScale(0.05, 0.05, 0.05)
        self.plane.setPos(0,0,0)
        self.cam.setPos(25.3, 2.26, 2.46)
        self.cam.lookAt(self.plane)

        self.plane.reparentTo(self.render)
        #self.taskMgr.add(self.set_coords, "get_coords")
        self.taskMgr.add(self.plane_coordiantes, "plane_coordiantes")
        
        #posInterval = time to move, finalPosition, startPosition
        posInterval1 = self.plane.posInterval(5, Point3(0, -6, -2), startPos=Point3(0,6,2))
        posInterval2 = self.plane.posInterval(5, Point3(0, 6, 2), startPos=Point3(0,-6,-2))

        self.planePace = Sequence(posInterval1, posInterval2, name = "planePace")
        self.planePace.loop()


    def set_coords(self, task):
        global points
        points = []
        if self.mouseWatcherNode.hasMouse():
            x = self.mouseWatcherNode.getMouseX()
            y = self.mouseWatcherNode.getMouseY()
            points = append(x,y)
        print(points)
        return Task.cont


    def plane_coordiantes(self, task):
        cam_coords = []
        cam_coords.append(self.plane.getPos())
        print(cam_coords)
        return Task.cont

class Events(DirectObject.DirectObject):
    def __init__(self):
        global points
        self.points = []
        self.accept("mouse1", self.set_coords)
        self.accept("escape", exit)
        #self.accept("mouse3", self.set_coords)



    def set_coords(self):
        points = []
        if base.mouseWatcherNode.hasMouse():
            x = base.mouseWatcherNode.getMouseX()
            y = base.mouseWatcherNode.getMouseY()
            points = append(x,y)
        print(points)
    
    
    
    def get_coords(self):
        points = self.get_coords()
        return points




app = MyApp()
h = Events()
app.run()
