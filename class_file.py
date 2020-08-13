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





#height is pulled from a NASA estimate

class Rocket:
    def __init__(self,cD = 0.3,diameter = 3.7,thrust = np.array([0,0,0]),mass = 1.,
                 stblzrforce = np.array([0,0,0]),surface_area = 0,height = 40.9,specific_heat = 1000,
                 payload_cg = np.array([0.,0.,0.]),payload_mass = 0,prop_cg = np.array([0.,0.,0.])):
        self.fuel = self.Fuel()
        self.thrust = thrust

        self.mass_empty = mass
        self.stblzrforce = stblzrforce
        self.surface_area = surface_area
        self.height = height
        self.cD = cD
        self.diameter = diameter
        self.prop_cg = prop_cg
        self.prop_mass = self.fuel.mass
        self.payload_mass = payload_mass
        self.area = pi*(self.diameter/2)**2
        self.cgpos_empty = np.array([0.,0.,self.height/2])
        self.specific_heat = specific_heat #just looked up some similar vals of material similar to the falcon 9 skin
        self.tot_cgpos_empty = (self.mass_empty*self.cgpos_empty + self.payload_mass*payload_cg)/(self.payload_mass+self.mass_empty)
        self.tot_cgpos = (self.tot_cgpos_empty*(self.mass_empty+self.payload_mass)+self.prop_cg*self.prop_mass)/(self.mass_empty+self.payload_mass + self.prop_mass)
        #around cg at all times (I have no clue if this is correct i sucked at this angular stuff help)
        self.moi_mat = np.array([[(self.prop_mass+self.mass_empty+self.payload_mass)*(3*(self.diameter/2)**2 + self.height**2)/12 + (self.mass_empty+self.payload_mass + self.prop_mass)*(self.tot_cgpos[1]**2 + self.tot_cgpos[2]**2),0.,0.],
                  [0.,(self.prop_mass+self.mass_empty+self.payload_mass)*(3*(self.diameter/2)**2 + self.height**2)/12 + (self.mass_empty+self.payload_mass + self.prop_mass)*(self.tot_cgpos[0]**2 + self.tot_cgpos[2]**2),0.],
                  [0.,0.,0.5*(self.prop_mass+self.mass_empty+self.payload_mass)*(self.diameter/2)**2 + (self.mass_empty+self.payload_mass + self.prop_mass)*(self.tot_cgpos[1]**2 + self.tot_cgpos[0]**2)]])

    def angular_accel_mat(self,torque):

        angular_accel = np.matmul(np.linalg.inv(self.moi_mat),torque)
        return angular_accel





    class GridFin:
        def __init__(self, angle, size):
            self.angle = angle
            self.size = size
    #loxtoprop is the ratio of liquid oxygen to the total fuel mass (which is lox + kerosene in the case of the F9)
    class Fuel:
        def __init__(self,mass = 0,loxtotot = 0.699,loxtemp = 66,fueltemp = 266.5,radius = 3.7/2):
            self.mass = mass

            self.loxtemp = loxtemp
            self.fueltemp = fueltemp
            self.loxdensity = 1255. #for LOX at 66K
            self.loxmass = loxtotot*mass
            self.propmass = (1-loxtotot)*mass
            self.propdensity = 818 #for RP1 at 266.5
            self.loxvol = self.loxmass/self.loxdensity
            self.propvol = self.propmass/self.propdensity
            self.radius = radius
         #to find the length of the container for the prop and LOX. This is not the total length but simply part of it
            self.loxlen = (self.loxvol-pi*self.radius**3 * 4/3)/(pi*self.radius**2)
            self.proplen = (self.propvol-pi*self.radius**3*2/3-(pi*self.radius**3-pi*self.radius**3*2/3))/(pi*self.radius**2)

        def prop_cg(self):
         #   v = pi*h*(3*r**2 - h**2)/3#formula for partial volume of a hemisphere
            if self.loxvol> self.loxlen*pi*self.radius**2 + pi*self.radius**3*2/3 and self.loxvol<=self.loxlen*pi*self.radius**2 + pi*self.radius**3*4/3:
                loxvol_upper = self.loxvol - self.loxlen*pi*self.radius**2 + pi*self.radius**3*2/3
                h = self.radius-3*loxvol_upper/(2*pi*self.radius**2)
                loxupper_mass = loxvol_upper*self.loxdensity
                zbar_upper = (3*(2*self.radius-h)**2/(4*(3*self.radius-h)))#cg pos if theres still lox in the upper hemisphere of the tank
                cg_lox_upper = np.array([0.,0.,((self.radius**4*self.loxdensity*pi/4)-(self.loxdensity*2*pi*self.radius**3/3-loxupper_mass)*zbar_upper)/loxupper_mass])

                lox_cg = ((cg_lox_upper+self.loxlen+self.radius)*loxupper_mass+(self.loxlen*pi*self.radius**2*self.loxdensity)*(self.loxlen/2 + self.radius) + (self.radius - 3*self.radius/8)*(pi*self.radius**3*2/3*self.loxdensity))/self.loxmass
            elif self.loxvol> pi*self.radius**3*2/3 and self.loxvol<= self.loxlen*pi*self.radius**2 + pi*self.radius**3*2/3:
                self.loxlen = (self.loxvol - pi * self.radius ** 3 * 2 / 3) / (pi * self.radius ** 2)

                lox_cg = np.array([0.,0.,((self.loxlen/2)*(self.loxvol - pi * self.radius ** 3 * 2 / 3)*self.loxdensity+(pi * self.radius ** 3 * 2 / 3)*self.loxdensity*(self.radius - 3*self.radius/8))/self.loxmass])

            elif self.loxvol > 0. and self.loxvol<=pi*self.radius**3*2/3:
                loxvol_lower = pi * self.radius ** 3 * 2 / 3
                h = self.radius - 3 * loxvol_lower / (2 * pi * self.radius ** 2)
                lox_cg = np.array([0.,0.,3*self.radius/8-(3*(2*self.radius-h)**2/(4*(3*self.radius-h)))])
            else:
                lox_cg = np.array([0.,0.,0.])

            if self.propvol> self.proplen*pi*self.radius**2 + pi*self.radius**3*2/3 and self.propvol <= self.proplen*pi*self.radius**2 + pi*self.radius**3*2/3 + self.radius**3*pi/3:
                pass







    class Engine:
        def __init__(self,thrust = 0,mdot = 0, pos = 0, isp = 330,gimbal = np.array([0.,0.,0.])):
            self.thrust = thrust
            self.isp = isp
            self.mdot = mdot
            self.pos = pos
            self.gimbal = gimbal

        def engine2bodygimbal(self,cgpos,pitch_mat,roll_mat,yaw_mat):
            new_thrust = np.matmul(yaw_mat,np.matmul(pitch_mat,np.matmul(roll_mat,self.thrust)))

            moments = np.matmul(np.array([[0,new_thrust[2],new_thrust[1]],
                                          [new_thrust[2],0,new_thrust[0]],
                                          [new_thrust[1],new_thrust[0],0]]),
                                abs(cgpos))

            return moments


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