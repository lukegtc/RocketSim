from class_file import *
import numpy as np
import matplotlib.pyplot as plt
from Atmosphere import *
from math import *
from rocket_pos_updater import earth,matrices,longitude,latitude,total_accel_final,rocket_vel_end_final,rocket_pos_end_final
import matplotlib as mpl
payload = 22800
rocket2 = Rocket(prop_mass=111000,mass_empty=5000,height=12.6,payload_mass=payload)

propellant_s2 = Propellant(mass=rocket2.prop_mass,radius=rocket2.radius,height=rocket2.height)
ve = np.array([0,0,1000.])
engine3 = Engine(mdot=934)
engine3.thrust = ve*engine3.mdot
rocket2.thrust = engine3.thrust
angular_vel_s2 = np.array([0.,0.,0.])
total_gimbal_s2 = np.array([0.,0.,0.])
skin_temp = 0
#From 1st stage
initial_accel = total_accel_final
initial_velocity = rocket_vel_end_final
initial_pos = rocket_pos_end_final
engine3.pos = initial_pos
print(engine3.pos)
running = True
points  =[]
time = []
dt = 0.1
t = 0
while running:
    tot_prop_cg_s2 = propellant_s2.prop_cg(rocket2.prop_mass)
    propellant_s2.tot_prop_cg = tot_prop_cg_s2
    rocket2.tot_cgpos = rocket2.tot_cgpos1(tot_prop_cg_s2)

    radius_s2 = sqrt((engine3.pos[0])**2 + (engine3.pos[1])**2 + (engine3.pos[2])**2)
    print(engine3.pos)

    press_s2, temp_s2, density_s2 = atmossolver(earth.g0, earth.layers, earth.avals,
                                          radius_s2 - earth.radius, earth.r,
                                          earth.base_temp,
                                          earth.base_press)
    #Pitch, roll and yaw matrices of the engine due to the gimbal
    roll_mat1 = matrices.roll_matrix(engine3.gimbal[0]) #around x axis
    pitch_mat1 = matrices.pitch_matrix(engine3.gimbal[1]) #around the y axis
    yaw_mat1 = matrices.yaw_matrix(engine3.gimbal[2])
    #MOI of 2nd Stage
    moi_mat_s2 = matrices.moi_matrix(rocket2.prop_mass,rocket2.mass_empty,rocket2.payload_mass,rocket2.diameter/2,rocket2.height,rocket2.tot_cgpos)

    moment_from_gimbal_engine_s2 = engine3.engine2bodygimbal(rocket2.tot_cgpos, pitch_mat1, roll_mat1, yaw_mat1)

    angular_accel1 = rocket2.angular_accel_mat(moi_mat_s2,moment_from_gimbal_engine_s2) #CHECK THIS

    angular_vel_s2 += angular_accel1*dt
    total_gimbal_s2 +=angular_vel_s2*dt

    tot_roll_mat1 = matrices.roll_matrix(total_gimbal_s2[0])
    tot_pitch_mat1 = matrices.pitch_matrix(total_gimbal_s2[1])
    tot_yaw_mat1 = matrices.yaw_matrix(total_gimbal_s2[2])


    #2nd Stage Thrust

    thrust3 = engine3.mdot*ve

    thrust3 = np.matmul(yaw_mat1,np.matmul(pitch_mat1,np.matmul(roll_mat1,thrust3)))
    rotated_thrust_s2 = np.matmul(tot_yaw_mat1, np.matmul(tot_pitch_mat1, np.matmul(tot_roll_mat1, thrust3)))

    long_yaw_mat1 = matrices.yaw_matrix(longitude)
    lat_roll_mat1 = matrices.roll_matrix(latitude)
    lat_pitch_mat1 = matrices.pitch_matrix(-pi/2)

    #Drag force change
    if radius_s2<earth.radius:
        drag_force_s2 = np.array([0.,0.,0.])
        rocket_vel_s2 = np.array([0.,0.,0.])
    else:
        if np.linalg.norm(rocket_vel_s2) != 0.:
            drag_force_s2 = -0.5*density_s2*(np.linalg.norm(rocket_vel_s2)**2)*rocket2.area*rocket2.cD*rocket_vel_s2/(np.linalg.norm(rocket_vel_s2))

        else:
            drag_force_s2 = np.array([0,0,0])

    total_force_s2 = (rocket2.mass_empty+rocket2.prop_mass+rocket2.payload_mass)*(engine3.pos/radius_s2)*-9.80665+rotated_thrust_s2+drag_force_s2

    rock_accel_s2 = total_force_s2/(rocket2.mass_empty+rocket2.prop_mass+rocket2.payload_mass) + initial_accel

    rocket_vel_s2+=rock_accel_s2*dt

    engine3.pos +=(np.array([465.1*cos(latitude),465.1*sin(latitude),0.])+rocket_vel_s2+initial_velocity)*dt

    t+=dt

    nose_radius_s2 = 0.7 #don't know if this is true but eh close enough for F9
    skin_temp = friction_temp(rocket_vel_s2,density_s2,nose_radius_s2,rocket2.specific_heat,rocket2.payload_mass+rocket2.mass_empty+rocket2.prop_mass,radius_s2,earth.radius,earth.layers,skin_temp)

    plotter_pos_s2 = np.matmul(long_yaw_mat1,np.matmul(lat_pitch_mat1,np.matmul(lat_roll_mat1, (engine3.pos))))
    tot_skin_temp = temp_s2 + skin_temp
    print(plotter_pos_s2)
    points.append(plotter_pos_s2)
    time.append(t)
    if skin_temp >= 1923.15:  #This reentry temp is totally made up based on space shuttle values
        running = False
        print('Burnt up on reentry or ascent')

    if t > 2500 or radius_s2<earth.radius:
        running = False
        break
    t += dt

plt.plot(points[:][0],points[:][1])
plt.show()
