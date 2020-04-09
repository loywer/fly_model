from direct.showbase.ShowBase import ShowBase
from direct.showbase import DirectObject
from direct.task import Task
from panda3d.core import LineSegs, NodePath, Point3D
from numpy import append



class MyApp(ShowBase):
    def __init__(self):
        ShowBase.__init__(self)
        self.taskMgr.add(self.get_coords, "get_coords")

    def get_coords(self, task):
        global points
        points = []
        if self.mouseWatcherNode.hasMouse():
            x = self.mouseWatcherNode.getMouseX()
            y = self.mouseWatcherNode.getMouseY()
            points = append(x,y)
        print(points)
        return Task.cont


class Hello(DirectObject.DirectObject):
    def __init__(self):
        global points
        self.points = []
        self.accept("mouse1", self.get_coords)
        #self.accept("mouse3", self.set_coords)



    def get_coords(self):
        points = []
        if base.mouseWatcherNode.hasMouse():
            x = base.mouseWatcherNode.getMouseX()
            y = base.mouseWatcherNode.getMouseY()
            points = append(x,y)
        print(points)
        return points
    
    
    
    def set_coords(self):
        points = self.get_coords()
        return points




app = MyApp()
h = Hello()
app.run()
