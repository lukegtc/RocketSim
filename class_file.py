from math import *
import numpy as np
'''
    Inputs:
        radius: radius of the planet
        g0: gravity of planet
        layers: altitudes in meters where the atmosperic layers start
        avals: a values for the atmospheric equations for each layer
        base temp: temperature at ground level
        base_press: temp at base of planet

'''
class Planet:
    def __init__(self,radius,g0 = 9.80665,layers = [0,11,20,32,47,51,71,86], avals = [-6.5,0,1,2.8,0,-2.8,-2.],base_temp = 0,base_press = 0,r = 0,GM = 0):

        self.radius = radius
        self.g0 = g0
        self.layers = layers
        self.avals = avals
        self.base_temp = base_temp
        self.base_press = base_press
        self.r = r
        self.GM = GM







class Rocket:
    def __init__(self,cD = 0.3,diameter = 3.7,thrust = np.array([0,0,0]),cgpos = np.array([0,0,0]),mass = 0,stblzrforce = np.array([0,0,0]),surface_area = 0,height = 0,specific_heat = 1000):
        self.thrust = thrust
        self.cgpos = cgpos
        self.mass_empty = mass
        self.stblzrforce = stblzrforce
        self.surface_area = surface_area
        self.height = height
        self.cD = cD
        self.diameter = diameter
        self.area = pi*(self.diameter/2)**2
        self.specific_heat = specific_heat #just looked up some similar vals of material similar to the falcon 9 skin

    class GridFin:
        def __init__(self, angle, size):
            self.angle = angle
            self.size = size

    class Fuel:
        def __init__(self,mass):
            self.mass = mass
            self.volume = mass
            pass


    class Engine:
        def __init__(self,thrust = 0,mdot = 0, pos = 0, isp = 330,gimbal = np.array([0.,0.,0.])):
            self.thrust = thrust
            self.isp = isp
            self.mdot = mdot
            self.pos = pos
            self.gimbal = gimbal

        def engine2bodygimbal(self,cgpos):



            bodygimbal = np.array([0,0,0])
            return bodygimbal
class Matrices:

    def roll_matrix(self,gamma):
        return np.array([[1,0,0],
                         [0,cos(gamma),-sin(gamma)],
                         [0,sin(gamma),cos(gamma)]])
    def pitch_matrix(self,beta):
        return np.array([[cos(beta),0,sin(beta)],
                         [0,1,0],
                         [-sin(beta),0,cos(beta)]])

    def yaw_matrix(self,alpha):
        return np.array([[cos(alpha),-sin(alpha),0],
                        [sin(alpha),cos(alpha),0],
                        [0,0,1]])