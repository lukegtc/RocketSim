from class_file import *
import numpy as np
import matplotlib.pyplot as plt
from Atmosphere import *
from math import *
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
payload = 22800
#rocket created
rocket1 = Rocket(prop_mass=418700.,mass_empty=27200.)
rocket2 = Rocket(prop_mass=111000,mass_empty=5000,height=12.6,payload_mass=payload)
#prop created
propellant = Propellant(mass=rocket1.prop_mass,radius=rocket1.radius,height=rocket1.height)


propellant1 = Propellant(mass=rocket2.prop_mass,radius=rocket2.radius,height=rocket2.height)
#engine Created
engine1 = Engine(mdot=845)
engine2 = Engine(mdot=845)
engine3 = Engine(mdot=934)
#rocket thrust force
rocket1.thrust = np.array([0,0,7605e3])
rocket2.thrust = np.array([0.,0.,934000.])
actual_latitude = 0
latitude = (actual_latitude-90)*pi/180
longitude = 0*pi/180
#exit velocity

ve = np.array([0,0,1000.])
dt = 0.1
t = 0
engine1.thrust = ve*engine1.mdot
engine2.thrust = ve*engine2.mdot
engine3.thrust = ve*engine3.mdot
#atmosphere values
layers = [0.,11.,20.,32.,47.,51.,71.,86.]
avals = [-6.5,0,1,2.8,0,-2.8,-2.]
base_temp = 288.15
base_press = 101325.
r = 287.
#Planet created
g0=9.80665
GM = 3.986004418e14
  # 28 degrees latitude
earth = Planet(6371e3,g0,layers,avals,base_temp,base_press,r,GM)
engine1.pos = np.array([0.,0.,0. + earth.radius]) #Changed
rocket_pos = np.array([0.,0.,0. + earth.radius])
engine3.pos = np.array([0.,0.,0.+rocket1.height]) #Changed+earth.radius

running =True
zpos = []
ypos = []
xpos = []
time1 = []
rads = []
time = []
ttlforce = []
dragforce = []
points = []
densitys = []
temps = []
skin_temps = []
presss = []
vels = []
cg_vals = []
moments= []
prop_mass = []
vels2 = []
tot_mass = []

#2nd stage values-------
tot_mass2 = []
throwaway = []
s2moments = []
pos1 = []
xpos1 = []
ypos1 = []
zpos1 = []
rads1 = []
velocities = []
accels = []
total_force1 = []
rocket_pos_end = []
checks = []
rocket_vel_end = []
s2drag = []
accel_vel_pos = []
poss = []
rocket_vel = np.array([0.,0.,0.]) #Dont forget to change some stuff about relative velocity and position ya know
rocket_vel1 = np.array([0.,0.,0.])
rocket_accel = np.array([0.,0.,0.])
rocket_vels = []
rocket_vels.append(rocket_vel)
rocket_poss = []
rocket_poss.append(rocket_pos)
rocket_pos_end_final = np.array([0.,0.,0.])
off_pad = False
#values for circular orbit
orbital_alt = 500e3
orbit_vel = sqrt(earth.GM/(earth.r+orbital_alt))
radius = earth.radius
skin_temp = 0
angular_vel = np.array([0.,0.,0.])
total_gimbal = np.array([0.,0.,0.])
#rocket_pos_end_final  = np.array([0.,0.,0.])
angular_vel1 = np.array([0.,0.,0.])
total_gimbal1 = np.array([0.,0.,0.])
phase1 = False
radius1 = earth.radius
while running:

    #changing cg due to propellant draining

    #1ST STAGE-----------------------------------------
    tot_prop_cg = propellant.prop_cg(rocket1.prop_mass)
    propellant.tot_prop_cg = tot_prop_cg
    rocket1.tot_cgpos = rocket1.tot_cgpos1(tot_prop_cg)


    #MASS OF STAGE PAYLOAD-----------------------------

    rocket1.payload_mass = rocket2.prop_mass+rocket2.payload_mass
    rocket1.payload_cg = rocket2.tot_cgpos

    #Distance to the center of the planet--------------
    #estimate from NASA estimated vals
    radius = sqrt(rocket_poss[-1][0]**2 + rocket_poss[-1][1]**2 + (rocket_poss[-1][2])**2)
    #print(radius)
    press,temp,density=atmossolver(earth.g0,earth.layers,earth.avals,radius-earth.radius,earth.r,earth.base_temp,earth.base_press)



    # if rocket1.prop_mass<= 0.95*418700:
    #     engine1.gimbal = np.array([-pi/2/180,0,0])
    #
    # else:
    #     engine1.gimbal[0] = 0
#    if rocket2.prop_mass<= 0.75*111000:
#        engine3.gimbal = np.array([-5*pi/2/180,0.,0.])
#    else:
#       engine3.gimbal = np.array([0.,0.,0.])
    # if engine1.gimbal[2]>= -np.pi/2:
    #     engine1.gimbal-=np.array([0,0,pi/2/1800])
    #
    # else:
    #     engine1.gimbal[2] = -pi/2

    #Stage 1 engine pitch roll yaw matrices

    matrices = Matrices()
    roll_mat = matrices.roll_matrix(engine1.gimbal[0]) #around x axis
    pitch_mat = matrices.pitch_matrix(engine1.gimbal[1]) #around the y axis
    yaw_mat = matrices.yaw_matrix(engine1.gimbal[2]) #around the z axis



    #moment of inertia------------------------------
    moi_mat = matrices.moi_matrix(rocket1.prop_mass,rocket1.mass_empty,rocket1.payload_mass,rocket1.diameter/2,rocket1.height,rocket1.tot_cgpos)

    moment_from_gimbal_engine = engine1.engine2bodygimbal(rocket1.tot_cgpos, pitch_mat, roll_mat, yaw_mat)

    angular_accel = rocket1.angular_accel_mat(moi_mat,moment_from_gimbal_engine) #CHECK THIS

    angular_vel += angular_accel*dt
    total_gimbal +=angular_vel*dt

    tot_roll_mat = matrices.roll_matrix(total_gimbal[0])
    tot_pitch_mat = matrices.pitch_matrix(total_gimbal[1])
    tot_yaw_mat = matrices.yaw_matrix(total_gimbal[2])


    #1st Stage Thrust w gimbal in engine1------------------------
    thrust1 = engine1.mdot * ve
    thrust1 = np.matmul(yaw_mat,np.matmul(pitch_mat,np.matmul(roll_mat,thrust1)))
    thrust2 = engine2.mdot*ve



    rocket1.prop_mass = rocket1.prop_mass - dt * engine1.mdot - dt * engine2.mdot * 8


    if rocket1.prop_mass <=0.:

        # total_force1.append(rocket_accel)
        # rocket_vel_end.append(rocket_vel)
        # rocket_pos_end.append(engine1.pos)
        # if len(total_force1) == 1:
        #     total_accel_final = total_force1[0]
        # if len(rocket_vel_end) == 1:
        #     rocket_vel_end_final = rocket_vel_end[0]
        # if len(rocket_pos_end) == 1:
        #     rocket_pos_end_final = rocket_pos_end[0]
        #     print(rocket_pos_end_final)
        check = 1
        engine1.mdot = 0
        engine2.mdot = 0
        rocket1.prop_mass = 0
        #2nd stage and payload detach from stage 1
        rocket1.payload_mass = 0
    else:
        check = 0
    checks.append(check)
    if radius<earth.radius:
        drag_force = np.array([0.,0.,0.])
        rocket_vels.append(np.array([0.,0.,0.]))
    else:
        if np.linalg.norm(rocket_vels[-1]) != 0.:
            drag_force = -0.5*density*(np.linalg.norm(rocket_vels[-1])**2)*rocket1.area*rocket1.cD*rocket_vels[-1]/(np.linalg.norm(rocket_vels[-1]))
            #fix something with the very last bit of the drag equation because its giving the wrong drag vectors i think invert the lat long stuff

        else:
            drag_force = np.array([0,0,0])
   # print(drag_force)

    #summation of all forces acting on the rocket

    rotated_thrust = np.matmul(tot_yaw_mat,np.matmul(tot_pitch_mat,np.matmul(tot_roll_mat,thrust2*8+thrust1)))

    total_force =(rocket1.mass_empty+rocket1.prop_mass+rocket1.payload_mass)*(rocket_poss[-1]/radius)*-9.80665+rotated_thrust+drag_force #Fix coordinate stuff
   # print((rocket1.mass_empty+rocket1.prop_mass+rocket1.payload_mass)*(rocket_poss[-1]/radius)*-9.80665)


    rocket_accel = total_force/(rocket1.mass_empty+rocket1.prop_mass+rocket1.payload_mass)
    #print(rocket_accel)
    rocket_vels.append(rocket_accel*dt+rocket_vels[-1])

    rocket_poss.append((np.array([465.1*cos(latitude),465.1*sin(latitude),0.])+rocket_vels[-1])*dt + rocket_poss[-1])

    engine1.pos = rocket_poss[-1]
    temps.append(temp)

    accel_vel_pos.append([rocket_accel,rocket_vels[-1],rocket_poss[-1],check])

    #print(accel_vel_pos)


    #change in temp due to atmosphere


    nose_radius = 0.7 #don't know if this is true but eh close enough for F9
    skin_temp = friction_temp(rocket_vels[-1],density,nose_radius,rocket1.specific_heat,rocket1.payload_mass+rocket1.mass_empty+rocket1.prop_mass,radius,earth.radius,earth.layers,skin_temp)


    tot_skin_temp = temp + skin_temp


    if skin_temp >= 1923.15:  #This reentry temp is totally made up based on space shuttle values
        running = False
        print('Burnt up on reentry or ascent')

    long_yaw_mat = matrices.yaw_matrix(longitude)
    lat_roll_mat = matrices.roll_matrix(latitude)
    lat_pitch_mat = matrices.pitch_matrix(-pi/2)
    plotter_pos = np.matmul(long_yaw_mat,np.matmul(lat_pitch_mat,np.matmul(lat_roll_mat, (engine1.pos)))) #
    points.append(plotter_pos)


    if 0.<rocket1.prop_mass:
        s2drag.append(np.linalg.norm(drag_force))
        zpos1.append(plotter_pos[2])
        ypos1.append(plotter_pos[1])
        xpos1.append(plotter_pos[0])

#plotting section----------------------------------------------------------------------------------------------
    zpos.append(plotter_pos[2])
    ypos.append(plotter_pos[1])
    xpos.append(plotter_pos[0])
    prop_mass.append(rocket1.prop_mass)
    ttlforce.append((np.linalg.norm(total_force)))
    skin_temps.append(tot_skin_temp)
    dragforce.append(np.linalg.norm(drag_force))
    time1.append(t)
    time.append(t)
    moments.append(moment_from_gimbal_engine[0])
#    s2moments.append(moment_from_gimbal_engine1[0])
    densitys.append(density)
    tot_mass.append(rocket1.mass_empty+rocket1.payload_mass+rocket1.prop_mass)
  #  tot_mass2.append(rocket2.prop_mass+rocket2.mass_empty+rocket2.payload_mass)
    temps.append(temp)
    presss.append(press)
    rads.append(radius-6371e3)
 #   rads1.append(radius1)
    vels.append(np.linalg.norm(rocket_vels[-1]))
   # vels2.append(rocket_vel/np.linalg.norm(rocket_vel))
    cg_vals.append(rocket1.tot_cgpos[2])
    if t > 2500 or np.linalg.norm(rocket_poss[-1])<earth.radius:
        running = False
        break
    t += dt

mpl.use('Qt5Agg')
#fig = plt.figure()
# Make data
# ax = fig.gca(projection='3d')
# u = np.linspace(0, 2 * np.pi, 100)
# v = np.linspace(0, np.pi, 100)
# x = 6371e3 * np.outer(np.cos(u), np.sin(v))
# y = 6371e3 * np.outer(np.sin(u), np.sin(v))
# z = 6371e3 * np.outer(np.ones(np.size(u)), np.cos(v))
# ax.plot_surface(x, y, z, color='b')
# ax.plot(xpos,ypos,zpos,label = 'position',color = 'r')
# ax.set_xlabel('X')
# ax.set_ylabel('Y')
# ax.set_zlabel('Z')
# plt.show()

plt.title('Stage 1 & 2 Mass')
plt.plot(time,tot_mass)
#plt.plot(time1,tot_mass2)
plt.show()


plt.title('Total Height')
plt.plot(time,rads)
#plt.plot(time,rads1)
plt.show()


fig, ax = plt.subplots()
plt.title('Y-Z Plot')
circle1 = plt.Circle((0,0),6371e3,color='b')
ax.add_artist(circle1)
plt.plot(ypos,zpos,color = 'r')
#plt.plot(ypos1,zpos1,color = 'g')
plt.xlabel('Y')
plt.ylabel('Z')
plt.show()


fig1, ax1 = plt.subplots()
plt.title('X-Z Plot')
circle2 = plt.Circle((0,0),6371e3, color='b')
ax1.add_artist(circle2)
plt.plot(xpos,zpos,color = 'r')
#plt.plot(xpos1,zpos1,color = 'g')
plt.xlabel('X')
plt.ylabel('Z')
plt.show()


fig2, ax2 = plt.subplots()
plt.title('X-Y Plot')
circle3 = plt.Circle((0,0),6371e3, color='b')
ax2.add_artist(circle3)
plt.plot(xpos,ypos,color = 'r')
#plt.plot(xpos1,ypos1,color = 'g')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()


plt.title('moments')
plt.plot(time,moments)
#plt.plot(time,s2moments)
plt.show()


plt.title('Drag Force')
plt.plot(time,dragforce,color = 'r')
#plt.plot(time,s2drag,color = 'g')
plt.show()


plt.title('Skin Temp')
plt.plot(time,skin_temps)
plt.show()


fig3,ax3 = plt.subplots()
plt.title('Second Stage')
circle4 = plt.Circle((0,0),6371e3, color='b')
ax3.add_artist(circle4)
plt.plot(xpos1,ypos1,color = 'r')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()


