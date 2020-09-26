import matplotlib.pyplot as plt
from class_file import *
from rocket_pos_updater import *
from thrust_gimbal_combo import timeset,random_val_test
mdot_vals,mdot_gimbal_vals,gimbal_vals = random_val_test(timeset)
points_s1,no_fuel_points_s1,dragforce,dynamic_press,time,moments,rads,tot_mass,accel_vel_pos,stage1,stage2,earth,latitude,longitude= rocket_func(
    Rocket(prop_mass=418700.,mass_empty=27200.,thrust=np.array([0,0,7605e3])),
    Rocket(prop_mass=111000,mass_empty=5000,height=12.6,payload_mass=22800,thrust=np.array([0.,0.,934000.])),
    Engine(mdot=845),Engine(mdot=845),mdot_vals,mdot_gimbal_vals,(-90)*pi/180,0*pi/180,gimbal_vals,
    np.array([0,0,1000.]),[0.,11.,20.,32.,47.,51.,71.,86.],[-6.5,0,1,2.8,0,-2.8,-2.],288.15,101325,287,9.80665,timeset)

#Stage 2 Motion after disconnection from Stage 1
points_s2,time_s2,rocket_vel_s2_set,eot,eot_t= staging_wo_1st_stage(stage2,np.array([0.,0.,1000.]),accel_vel_pos,latitude,longitude,earth,Engine(mdot=934),-30)


#This is the motion of the initial point where the craft launched from
target_track,init_time = moving_target(time,np.array([0.,0.,earth.radius]),latitude,longitude,earth.radius)



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
plt.plot(points_s1[:,1],points_s1[:,2],color = 'r')
plt.plot(target_track[:,1],target_track[:,2],color = 'g')
plt.xlabel('Y')
plt.ylabel('Z')
plt.show()

fig1, ax1 = plt.subplots()
plt.title('X-Z Plot')
circle2 = plt.Circle((0,0),6371e3, color='b')
ax1.add_artist(circle2)
plt.plot(points_s1[:,0],points_s1[:,2],color = 'r')
plt.plot(target_track[:,0],target_track[:,2],color = 'g')
plt.xlabel('X')
plt.ylabel('Z')
plt.show()

fig2, ax2 = plt.subplots()
plt.title('X-Y Plot')
circle3 = plt.Circle((0,0),6371e3, color='b')
ax2.add_artist(circle3)
plt.plot(points_s1[:,0],points_s1[:,1],color = 'r')
plt.plot(target_track[:,0],target_track[:,1],color = 'g')
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



#Stage Two Graphing Equations
fig_s2, ax_s2 = plt.subplots()
plt.title('X-Y Plot')
circle1 = plt.Circle((0,0),6371e3,color='b')
ax_s2.add_artist(circle1)
#plt.plot(xpos1,ypos1,color  = 'g')
plt.plot(points_s2[:,0],points_s2[:,1],color = 'r')
plt.plot(points_s1[:,0],points_s1[:,1],color = 'g')
plt.plot(eot[:,0],eot[:,1],color = 'b')
#plt.plot(ypos1,zpos1,color = 'g')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()