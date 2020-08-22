from class_file import *
import numpy as np
import matplotlib.pyplot as plt
from Atmosphere import *
from math import *
from rocket_pos_updater import earth,longitude,latitude,accel_vel_pos,rocket2,xpos1,ypos1,zpos1
import matplotlib as mpl
payload = 22800
#rocket2 = Rocket(prop_mass=111000,mass_empty=5000,height=12.6,payload_mass=payload)
matrices1 = Matrices()
propellant_s2 = Propellant(mass=rocket2.prop_mass,radius=rocket2.radius,height=rocket2.height)
ve = np.array([0,0,1000.])
engine3 = Engine(mdot=934)
engine3.gimbal = np.array([0.,0.,0.])
engine3.thrust = ve*engine3.mdot
rocket2.thrust = engine3.thrust
angular_vel_s2 = np.array([0.,0.,0.])
total_gimbal_s2 = np.array([0.,0.,0.])
skin_temp = 0
#From 1st stage
#initial_accel = np.array([0.,0.,0.])
#initial_velocity = np.array([0.,0.,0.])
#initial_pos = np.array([0.,0.,0. + earth.radius])
#engine3.pos = initial_pos
rocket_vel_s2_set = []
rocket_pos_s2_set = []
angular_vel_s2 = []
total_gimbal_s2  =[]
for item in accel_vel_pos:

    if item[-1] == 1:

        check = 1
        s2_initial = item
        initial_accel = s2_initial[0]
        initial_velocity =s2_initial[1]
        initial_pos = s2_initial[2]
        rocket_vel_s2_set.append(initial_velocity)
        rocket_pos_s2_set.append(initial_pos)
        initial_angular_accel = s2_initial[5]
        initial_angular_vel = s2_initial[3]
        angular_vel_s2.append(initial_angular_vel)
        #initial_gimbal= s2_initial[4]
        break
    else:
        check = 0
        initial_pos = np.array([0.,0.,0.])

long_yaw_mat1 = matrices1.yaw_matrix(longitude)
lat_roll_mat1 = matrices1.roll_matrix(latitude)
lat_pitch_mat1 = matrices1.pitch_matrix(-pi / 2)
initial_pos_rotated = np.matmul(long_yaw_mat1,np.matmul(lat_pitch_mat1,np.matmul(lat_roll_mat1,initial_pos)))
xpos_s2 = []
ypos_s2 = []
zpos_s2 = []
xpos_s2.append(initial_pos_rotated[0])
ypos_s2.append(initial_pos_rotated[1])
zpos_s2.append(initial_pos_rotated[2])
engine3.pos = initial_pos
running = True
points  =[]
time = []
dt = 0.1
t = 0
while running:
    if check == 0:
        break
    tot_prop_cg_s2 = propellant_s2.prop_cg(rocket2.prop_mass)
    propellant_s2.tot_prop_cg = tot_prop_cg_s2
    rocket2.tot_cgpos = rocket2.tot_cgpos1(tot_prop_cg_s2)

    radius_s2 = sqrt((rocket_pos_s2_set[-1][0])**2 + (rocket_pos_s2_set[-1][1])**2 + (rocket_pos_s2_set[-1][2])**2)
    #print(radius_s2)

    if rocket2.prop_mass<= 0.99*111000:
        engine3.gimbal = np.array([-1*pi/2/180,0.,0.])
    else:
       engine3.gimbal = np.array([0.,0.,0.])


    press_s2, temp_s2, density_s2 = atmossolver(earth.g0, earth.layers, earth.avals,radius_s2 - earth.radius, earth.r,earth.base_temp,earth.base_press)

    #Pitch, roll and yaw matrices of the engine due to the gimbal

    roll_mat1 = matrices1.roll_matrix(engine3.gimbal[0]) #around x axis
    pitch_mat1 = matrices1.pitch_matrix(engine3.gimbal[1]) #around the y axis
    yaw_mat1 = matrices1.yaw_matrix(engine3.gimbal[2])
    #MOI of 2nd Stage
    moi_mat_s2 = matrices1.moi_matrix(rocket2.prop_mass,rocket2.mass_empty,rocket2.payload_mass,rocket2.diameter/2,rocket2.height,rocket2.tot_cgpos)

    moment_from_gimbal_engine_s2 = engine3.engine2bodygimbal(rocket2.tot_cgpos, pitch_mat1, roll_mat1, yaw_mat1)

    angular_accel1 = rocket2.angular_accel_mat(moi_mat_s2,moment_from_gimbal_engine_s2) #CHECK THIS

    angular_vel_s2.append(angular_accel1*dt+angular_vel_s2[-1])
    total_gimbal_s2.append(angular_vel_s2[-1]*dt)
  # print(angular_vel_s2)
    tot_roll_mat1 = matrices1.roll_matrix(total_gimbal_s2[-1][0])
    tot_pitch_mat1 = matrices1.pitch_matrix(total_gimbal_s2[-1][1])
    tot_yaw_mat1 = matrices1.yaw_matrix(total_gimbal_s2[-1][2])


    #2nd Stage Thrust

    thrust3 = engine3.mdot*ve

    thrust3 = np.matmul(yaw_mat1,np.matmul(pitch_mat1,np.matmul(roll_mat1,thrust3)))



    rocket2.prop_mass = rocket2.prop_mass - dt * engine3.mdot
    if rocket2.prop_mass <=0.:
        engine3.mdot = 0.
        rocket2.prop_mass = 0.


    #Drag force change
    if radius_s2<earth.radius:
        drag_force_s2 = np.array([0.,0.,0.])
        rocket_vel_s2_set.append(np.array([0.,0.,0.]))
        break
    else:

        if np.linalg.norm(rocket_vel_s2_set[-1]) != 0.:
            drag_force_s2 = -0.5*density_s2*(np.linalg.norm(rocket_vel_s2_set[-1])**2)*rocket2.area*rocket2.cD*rocket_vel_s2_set[-1]/(np.linalg.norm(rocket_vel_s2_set[-1]))

        else:
            drag_force_s2 = np.array([0,0,0])

    rotated_thrust_s2 = np.matmul(tot_yaw_mat1, np.matmul(tot_pitch_mat1, np.matmul(tot_roll_mat1, thrust3)))

    total_force_s2 = (rocket2.mass_empty+rocket2.prop_mass+rocket2.payload_mass)*(rocket_pos_s2_set[-1]/radius_s2)*-9.80665+rotated_thrust_s2+drag_force_s2

    rock_accel_s2 = total_force_s2/(rocket2.mass_empty+rocket2.prop_mass+rocket2.payload_mass)
   # print(rock_accel_s2)
    rocket_vel_s2_set.append(rock_accel_s2*dt+rocket_vel_s2_set[-1])

    rocket_pos_s2_set.append((np.array([465.1*cos(latitude),465.1*sin(latitude),0.])+rocket_vel_s2_set[-1])*dt + rocket_pos_s2_set[-1])
    engine3.pos = rocket_pos_s2_set[-1]
   # print(rocket_pos_s2_set[-1])
    t+=dt

    nose_radius_s2 = 0.7 #don't know if this is true but eh close enough for F9
  #  skin_temp = friction_temp(rocket_vel_s2,density_s2,nose_radius_s2,rocket2.specific_heat,rocket2.payload_mass+rocket2.mass_empty+rocket2.prop_mass,radius_s2,earth.radius,earth.layers,skin_temp) +
  #  print(skin_temp)

    plotter_pos_s2 = np.matmul(long_yaw_mat1,np.matmul(lat_pitch_mat1,np.matmul(lat_roll_mat1, rocket_pos_s2_set[-1])))
    tot_skin_temp = temp_s2 + skin_temp

    points.append(plotter_pos_s2)
    time.append(t)
    xpos_s2.append(plotter_pos_s2[0])
    ypos_s2.append(plotter_pos_s2[1])
    zpos_s2.append(plotter_pos_s2[2])
  #  print(rocket_pos_s2_set[-1])
    if skin_temp >= 1923.15:  #This reentry temp is totally made up based on space shuttle values
        running = False
        print('Burnt up on reentry or ascent')

    if t > 2500 or radius_s2<earth.radius:
        running = False
        break
    t += dt


fig_s2, ax_s2 = plt.subplots()
plt.title('X-Y Plot')
circle1 = plt.Circle((0,0),6371e3,color='b')
ax_s2.add_artist(circle1)
#plt.plot(xpos1,ypos1,color  = 'g')
plt.plot(xpos_s2,ypos_s2,color = 'r')
#plt.plot(ypos1,zpos1,color = 'g')
plt.xlabel('X')
plt.ylabel('Y')
plt.show()