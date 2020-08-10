from class_file import *
import numpy as np
import matplotlib.pyplot as plt
from Atmosphere import *
from math import *
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D

#rocket created
rocket1 = Rocket()
rocket1.mass_empty = 27200.
#engine Created
engine1 = rocket1.Engine(mdot=845)
engine2 = rocket1.Engine(mdot=845)
#rocket thrust force
rocket1.thrust = np.array([0,0,7605e3])
#exit velocity
ve = np.array([0,0,700])
dt = 0.1
t = 0

#atmosphere values
layers = [0.,11.,20.,32.,47.,51.,71.,86.]
avals = [-6.5,0,1,2.8,0,-2.8,-2.]
base_temp = 288.15
base_press = 101325.
r = 287.
#Planet created
g0=9.80665
GM = 3.986004418e14
earth = Planet(6371e3,g0,layers,avals,base_temp,base_press,r,GM)
engine1.pos = np.array([0.,0.,0. + earth.radius])
#fuel created
fuel = rocket1.Fuel(418700)


running =True
zpos = [6371e3,6371e3]
ypos = [0,0]
xpos = [0,0]
time1 = [0,0]
rads = []
time = []
ttlforce = []
dragforce = []
densitys = []
temps = []
skin_temps = []
presss = []
rocket_vel = np.array([0.,0.,0.])
off_pad = False
#values for circular orbit
orbital_alt = 500e3
orbit_vel = sqrt(earth.GM/(earth.r+orbital_alt))
radius = earth.radius
skin_temp = 0
angular_vel = np.array([0.,0.,0.])
total_gimbal = np.array([0.,0.,0.])
while running:
     #estimate from NASA estimated vals
    radius = sqrt(engine1.pos[0]**2 + engine1.pos[1]**2 + (engine1.pos[2])**2)
    # theta = np.arctan(sqrt(engine1.pos[0]**2 + engine1.pos[1]**2)/(engine1.pos[2]+earth.radius))
    # psi = np.arctan(engine1.pos[1]/engine1.pos[0])
   # cart2sphere = np.array([[sin(theta)*cos(psi),sin(theta)*sin(psi),cos(theta)],
   #                         [cos(theta)*cos(psi),cos(theta)*sin(psi),-sin(theta)],
    #                        [-sin(theta),cos(psi),0]])
    #calculates pressure, temperature and density using the atmossolver equation from the atmosphere file
    press,temp,density=atmossolver(earth.g0,earth.layers,earth.avals,radius-earth.radius,earth.r,earth.base_temp,earth.base_press)


    engine1.isp = ve/earth.g0
    #engine1 gimballing
    #Relation of gimballed engine to pitch, roll and yaw angle of rocket
    #This is the gimbal of the engine
    if engine1.gimbal[0] >=-np.pi:

        engine1.gimbal-= np.array([pi/2/1800,0,0])

    else:
        engine1.gimbal[0] = -pi/2
    if engine1.gimbal[2]>= -np.pi/2:
        engine1.gimbal-=np.array([0,0,pi/2/1800])

    else:
        engine1.gimbal[2] = -pi/2
    matrices = Matrices()
    roll_mat = matrices.roll_matrix(engine1.gimbal[0]) #around x axis
    pitch_mat = matrices.pitch_matrix(engine1.gimbal[1]) #around the y axis
    yaw_mat = matrices.yaw_matrix(engine1.gimbal[2]) #around the z axis

    #this is the gimbal of the rocket

    moment_from_gimbal_engine = engine1.engine2bodygimbal(rocket1.tot_cgpos,pitch_mat,roll_mat,yaw_mat)
    angular_accel = rocket1.angular_accel_mat(moment_from_gimbal_engine)
    angular_vel += angular_accel*dt
    total_gimbal +=angular_vel*dt

    tot_roll_mat = matrices.roll_matrix(total_gimbal[0])
    tot_pitch_mat = matrices.pitch_matrix(total_gimbal[1])
    tot_yaw_mat = matrices.yaw_matrix(total_gimbal[2])
    #fuel mass change
    if fuel.mass <=0:
        engine1.mdot = 0
        engine2.mdot = 0
        fuel.mass = 0

    thrust1 = engine1.mdot * ve
    thrust2 = engine2.mdot*ve
    fuel.mass=fuel.mass - dt*engine1.mdot-dt*engine2.mdot*8

    #Drag force change
    if radius<earth.radius:
        drag_force = np.array([0.,0.,0.])
        rocket_vel = np.array([0.,0.,0.])
    else:
        if np.linalg.norm(rocket_vel) != 0.:
            drag_force = -0.5*density*(np.linalg.norm(rocket_vel)**2)*rocket1.area*rocket1.cD*rocket_vel/(np.linalg.norm(rocket_vel))
        else:
            drag_force = np.array([0,0,0])
    # if np.linalg.norm(rocket_vel) == 0:
    #     norm_vel = 1
    # else:
    #     norm_vel = np.linalg.norm(rocket_vel)
    #summation of all forces acting on the rocket
    total_force = (rocket1.mass_empty+fuel.mass)*(engine1.pos/radius)*-9.80665+np.matmul(tot_yaw_mat,np.matmul(tot_pitch_mat,np.matmul(tot_roll_mat,thrust2*8+thrust1)))+drag_force


    rocket_accel = total_force/(rocket1.mass_empty+fuel.mass)

    rocket_vel +=rocket_accel*dt
    engine1.pos +=rocket_vel*dt
    temps.append(temp)

    t+=dt
    #change in temp due to atmosphere
    nose_radius = 0.7 #don't know if this is true but eh close enough for F9
    if radius<= 6371e3 + earth.layers[-1]*10**3:
        #added temp per dt
        fric_temp = 1.83*10**(-4)*25*np.linalg.norm(rocket_vel)**3*sqrt(density/nose_radius)/rocket1.specific_heat/(rocket1.mass_empty+fuel.mass) #the 25 sqm is the SA of the nose cone exposed to the air friction
        skin_temp+=fric_temp
    else:
        if skin_temp >=2:
            skin_temp-=  0.1 #honestly have no clue if this is how fast heat dissipates

    tot_skin_temp = temp + skin_temp


    if skin_temp >= 1923.15:  #This reentry temp is totally made up based on space shuttle values
        running = False
        print('Burnt up on reentry or ascent')

    if t > 7000 or radius<earth.radius:
        running = False
        break
    zpos.append(engine1.pos[2])
    ypos.append(engine1.pos[1])
    xpos.append((engine1.pos[0]))
    ttlforce.append((np.linalg.norm(total_force)))
    skin_temps.append(tot_skin_temp)
    dragforce.append(np.linalg.norm(drag_force))
    time1.append(t)
    time.append(t)
    densitys.append(density)
    temps.append(temp)
    presss.append(press)
    rads.append(radius-6371e3)
mpl.use('Qt5Agg')
fig = plt.figure()
# Make data
ax = fig.gca(projection='3d')
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 6371e3 * np.outer(np.cos(u), np.sin(v))
y = 6371e3 * np.outer(np.sin(u), np.sin(v))
z = 6371e3 * np.outer(np.ones(np.size(u)), np.cos(v))

# Plot the surface
ax.plot_surface(x, y, z, color='b')

ax.plot(xpos,ypos,zpos,label = 'position',color = 'r')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.show()
plt.title('ttlforce')
plt.plot(time,rads)
plt.show()
fig, ax = plt.subplots()
plt.title('Y-Z Plot')
circle1 = plt.Circle((0,0),6371e3,color='b')
ax.add_artist(circle1)
plt.plot(ypos,zpos,color = 'r')
plt.show()

fig1, ax1 = plt.subplots()
plt.title('X-Z Plot')
circle2 = plt.Circle((0,0),6371e3, color='b')

ax1.add_artist(circle2)
plt.plot(xpos,zpos,color = 'r')
plt.show()
plt.title('Altitude')
plt.plot(time,rads)
plt.show()
plt.title('Drag Force')
plt.plot(time,dragforce,color = 'r')
plt.show()
plt.title('Skin Temp')
plt.plot(time,skin_temps)
plt.show()
# plt.title('Density')
# plt.plot(time,densitys)
# plt.show()
# plt.title('temp')
# plt.plot(time,temps)
# plt.show()




#
# plt.show()



