from direct.showbase.ShowBase import ShowBase
from direct.showbase import DirectObject
from panda3d.core import LineSegs, NodePath, Point3D
from math import cos, sin, pi
from numpy import append

class Hello(DirectObject.DirectObject):
    def __init__(self):
        global points
        self.points = []
        self.accept("mouse1", self.get_coords)
        self.accept("mouse3", self.set_coords)



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




base = ShowBase()
test = Hello()
base.run()