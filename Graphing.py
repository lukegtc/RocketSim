import matplotlib.pyplot as plt
from class_file import *
from rocket_pos_updater import *

xpos,ypos,zpos,xpos1,ypos1,zpos1,dragforce,dynamic_press,time,moments,rads,rads,tot_mass,accel_vel_pos,stage1,stage2,earth,latitude,longitude = rocket_func(
    Rocket(prop_mass=418700.,mass_empty=27200.,thrust=np.array([0,0,7605e3])),
    Rocket(prop_mass=111000,mass_empty=5000,height=12.6,payload_mass=22800,thrust=np.array([0.,0.,934000.])),
    Engine(mdot=845),Engine(mdot=845),(-90)*pi/180,0*pi/180,0,
    np.array([0,0,1000.]),[0.,11.,20.,32.,47.,51.,71.,86.],[-6.5,0,1,2.8,0,-2.8,-2.],288.15,101325,287,9.80665)

#Stage 2 Motion after disconnection from Stage 1
points,time_s2,rocket_vel_s2_set,xpos_s2,ypos_s2,zpos_s2 = staging_wo_1st_stage(stage2,np.array([0.,0.,1000.]),accel_vel_pos,latitude,longitude,earth)

mpl.use('Qt5Agg')

plt.title('Stage 1 & 2 Mass')
plt.plot(time,tot_mass)

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
plt.plot(time,dynamic_press,color = 'g')
#plt.plot(time,s2drag,color = 'g')
plt.show()

fig3,ax3 = plt.subplots()
plt.title('Second Stage')
circle4 = plt.Circle((0,0),6371e3, color='b')
ax3.add_artist(circle4)
plt.plot(xpos1,ypos1,color = 'r')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()

















#Stage Two Graphing Equations
fig_s2, ax_s2 = plt.subplots()
plt.title('X-Y Plot')
circle1 = plt.Circle((0,0),6371e3,color='b')
ax_s2.add_artist(circle1)
#plt.plot(xpos1,ypos1,color  = 'g')
plt.plot(xpos_s2,ypos_s2,color = 'r')
plt.plot(xpos,ypos,color = 'g')
#plt.plot(ypos1,zpos1,color = 'g')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()