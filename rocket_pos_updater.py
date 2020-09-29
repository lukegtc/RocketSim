from class_file import *
import numpy as np
import matplotlib.pyplot as plt
from Atmosphere import *
from math import *
from tqdm import tqdm
import matplotlib as mpl
def rocket_func(stage1,stage2,gimbal_engine1,engine1,latitude,longitude,exit_vel,layers,avals,base_temp,base_press,r,g0,gimbal_throttle_pcts,throttle_pcts,timeset,gimbals):
    GM = 3.986004418e14
    earth = Planet(6371e3,g0,layers,avals,base_temp,base_press,r,GM)
    gimbal_engine1.pos = np.array([0.,0.,0. + earth.radius])

    dt = 0.1
    t = 0.
    i = 0
    propellant1 = Propellant(mass=stage1.prop_mass,radius=stage1.radius,height=stage1.height)
    running =True


    rads = []
    time = []
    ttlforce = []
    dragforce = []
    points = []
    s2_points_on_s1 = []
    densitys = []
    temps = []
    presss = []
    vels = []
    cg_vals = []
    moments= []
    prop_mass = []
    tot_mass = []

    #2nd stage values-------

    tang_vel = []
    checks = []
    s2drag = []
    accel_vel_pos = []
    angular_vel_mat = []
    total_gimbal_mat = []

    rocket_vel = np.array([0.,0.,0.]) #Dont forget to change some stuff about relative velocity and position ya know

    rocket_vels = []
    rocket_vels.append(rocket_vel)
    rocket_poss = []
    dynamic_press = []
    rocket_poss.append(gimbal_engine1.pos)
    prop_masses =[]
    rocket_calc_pos=[]
    rocket_calc_pos.append(gimbal_engine1.pos)
    angular_vel_mat.append(np.array([0.,0.,0.]))
    total_gimbal_mat.append(np.array([0.,0.,0.]))
    drag_force = np.array([0.,0.,0.])
    total_force = np.array([0.,0.,0.])
    while running:

        gimbal_engine1.thrust = exit_vel * gimbal_engine1.mdot*gimbal_throttle_pcts[i]
        engine1.thrust = exit_vel * engine1.mdot * throttle_pcts[i]
        #changing cg due to propellant draining

        #1ST STAGE-----------------------------------------
        tot_prop_cg = propellant1.prop_cg(stage1.prop_mass)
        propellant1.tot_prop_cg = tot_prop_cg
        stage1.tot_cgpos = stage1.tot_cgpos1(tot_prop_cg)


        #MASS OF STAGE PAYLOAD-----------------------------

        stage1.payload_mass = stage2.prop_mass+stage2.payload_mass
        stage1.payload_cg = stage2.tot_cgpos

        #Distance to the center of the planet--------------
        #estimate from NASA estimated vals
        radius = sqrt(rocket_poss[-1][0]**2 + rocket_poss[-1][1]**2 + (rocket_poss[-1][2])**2)
        tangential_vel = rocket_vels[-1] - np.dot(rocket_vels[-1],rocket_poss[-1])/np.sqrt(sum(rocket_poss[-1]**2)) * rocket_poss[-1]



        press,temp,density=atmossolver(earth.g0,earth.layers,earth.avals,radius-earth.radius,earth.r,earth.base_temp,earth.base_press)

        #Gimballing rules

        #if stage1.prop_mass<= 0.75*418700:
        gimbal_engine1.gimbal = gimbals[i]

      #  else:
        #    gimbal_engine1.gimbal[0] = 0



        #Stage 1 engine pitch roll yaw matrices

        matrices = Matrices()
        roll_mat = matrices.roll_matrix(gimbal_engine1.gimbal[0]) #around x axis
        pitch_mat = matrices.pitch_matrix(gimbal_engine1.gimbal[1]) #around the y axis
        yaw_mat = matrices.yaw_matrix(gimbal_engine1.gimbal[2]) #around the z axis



        #moment of inertia------------------------------
        moi_mat = matrices.moi_matrix(stage1.prop_mass,stage1.mass_empty,stage1.payload_mass,stage1.diameter/2,stage1.height,stage1.tot_cgpos)

        moment_from_gimbal_engine = gimbal_engine1.engine2bodygimbal(stage1.tot_cgpos, pitch_mat, roll_mat, yaw_mat)

        drag_pos = stage1.tot_cgpos-np.array([0,0,stage1.height])
        new_drag_force = np.matmul(yaw_mat,np.matmul(pitch_mat,np.matmul(roll_mat,drag_force)))

        drag_moment = np.matmul(np.array([[0,new_drag_force[2],new_drag_force[1]],
                                      [new_drag_force[2],0,new_drag_force[0]],
                                      [new_drag_force[1],new_drag_force[0],0]]),
                            abs(drag_pos))

        total_moment = moment_from_gimbal_engine#-drag_moment
        angular_accel = stage1.angular_accel_mat(moi_mat,total_moment) #CHECK THIS

        angular_vel_mat.append(angular_accel*dt+angular_vel_mat[-1])
        total_gimbal_mat.append(angular_vel_mat[-1]*dt + total_gimbal_mat[-1])

        tot_roll_mat = matrices.roll_matrix(total_gimbal_mat[-1][0])
        tot_pitch_mat = matrices.pitch_matrix(total_gimbal_mat[-1][1])
        tot_yaw_mat = matrices.yaw_matrix(total_gimbal_mat[-1][2])

        #Model moment the drag creates



        #1st Stage Thrust w gimbal in engine1------------------------
        thrust1 = gimbal_engine1.thrust
        thrust1 = np.matmul(yaw_mat,np.matmul(pitch_mat,np.matmul(roll_mat,thrust1)))
        thrust2 = engine1.thrust

        stage1.prop_mass = stage1.prop_mass - dt * gimbal_engine1.mdot*gimbal_throttle_pcts[i] - dt * engine1.mdot*throttle_pcts[i] * 8

        if stage1.prop_mass <=0.:

            check = 1
            gimbal_engine1.mdot = 0
            engine1.mdot = 0
            stage1.prop_mass = 0
            #2nd stage and payload detach from stage 1
            stage1.payload_mass = 0
        else:
            check = 0
        checks.append(check)
        if radius<earth.radius:
            drag_force = np.array([0.,0.,0.])
            rocket_vels.append(np.array([0.,0.,0.]))
            dynamic_press.append(0.)
        else:
            if np.linalg.norm(rocket_vels[-1]) != 0.:
                drag_force = -0.5*density*(np.linalg.norm(rocket_vels[-1])**2)*stage1.area*stage1.cD*rocket_vels[-1]/(np.linalg.norm(rocket_vels[-1]))
                #fix something with the very last bit of the drag equation because its giving the wrong drag vectors i think invert the lat long stuff
                dynamic_press.append(0.5*density*(np.linalg.norm(rocket_vels[-1])**2))
            else:
                drag_force = np.array([0,0,0])
                dynamic_press.append(0.)



        #summation of all forces acting on the rocket

        rotated_thrust = np.matmul(tot_yaw_mat,np.matmul(tot_pitch_mat,np.matmul(tot_roll_mat,thrust2*8+thrust1)))

        total_force =(stage1.mass_empty+stage1.prop_mass+stage1.payload_mass)*(rocket_calc_pos[-1]/radius)*-9.80665+rotated_thrust+drag_force #Fix coordinate stuff


        rocket_accel = total_force/(stage1.mass_empty+stage1.prop_mass+stage1.payload_mass)
        rocket_vels.append(rocket_accel*dt+rocket_vels[-1])

        rocket_poss.append((np.array([465.1*cos(latitude),465.1*sin(latitude),0.])+rocket_vels[-1])*dt + rocket_poss[-1])
        #initial_launch_point.append(np.array([465.1*cos(latitude),465.1*sin(latitude),0.])*dt+initial_launch_point[-1])
        #THIS ROCKET_CALC_POS IS USED FOR THE CALCULATIONS ONLY, I STILL NEED TO FIND A PROPER WAY TO
        rocket_calc_pos.append(rocket_vels[-1]*dt+rocket_calc_pos[-1])
        #used to be rocket_poss[-1]
        gimbal_engine1.pos = rocket_poss[-1]
        temps.append(temp)

        accel_vel_pos.append([rocket_accel,rocket_vels[-1],rocket_poss[-1],angular_vel_mat[-1],total_gimbal_mat[-1],angular_accel,check])
        tang_vel.append(np.linalg.norm(rocket_vels[-1]-(np.dot(rocket_vels[-1],rocket_poss[-1])/np.sqrt(rocket_poss[-1]**2))*rocket_poss[-1]))
        long_yaw_mat = matrices.yaw_matrix(longitude)
        lat_roll_mat = matrices.roll_matrix(latitude)
        lat_pitch_mat = matrices.pitch_matrix(-pi/2)
        plotter_pos = np.matmul(long_yaw_mat,np.matmul(lat_pitch_mat,np.matmul(lat_roll_mat, (rocket_poss[-1]))))
       # init_point = np.matmul(long_yaw_mat,np.matmul(lat_pitch_mat,np.matmul(lat_roll_mat, (initial_launch_point[-1]))))
        points.append(plotter_pos)

        if 0.<stage1.prop_mass:

            s2drag.append(np.linalg.norm(drag_force))
            s2_points_on_s1.append(plotter_pos)

        #plotting section----------------------------------------------------------------------------------------------

        prop_mass.append(stage1.prop_mass)
        ttlforce.append((np.linalg.norm(total_force)))

        dragforce.append(np.linalg.norm(drag_force))

        time.append(t)
        moments.append(total_moment[0])
        densitys.append(density)
        tot_mass.append(stage1.mass_empty+stage1.payload_mass+stage1.prop_mass)
        prop_masses.append(stage1.prop_mass)
        temps.append(temp)
        presss.append(press)
        rads.append(radius-6371e3)
        vels.append(np.linalg.norm(rocket_vels[-1]))
        cg_vals.append(stage1.tot_cgpos[2])

        if timeset[i] > 2500 or np.linalg.norm(rocket_poss[-1])<earth.radius:
            running = False
            break
        i+=1
        t+=dt
    return np.array(points),np.array(s2_points_on_s1),dragforce,dynamic_press,time,moments,rads,tot_mass,accel_vel_pos,stage1,stage2,earth,latitude,longitude,np.array(tang_vel),np.array(prop_masses)

def staging_wo_1st_stage(stage,exit_vel,initial_vals,latitude,longitude,planet,engine_s2,tot_gimbal):

    matrices1 = Matrices()
    propellant_s2 = Propellant(mass=stage.prop_mass,radius=stage.radius,height=stage.height)
    ve = exit_vel

    engine_s2.gimbal = np.array([0.,0.,0.])
    engine_s2.thrust = ve*engine_s2.mdot
    stage.thrust = engine_s2.thrust

    rocket_vel_s2_set = []
    rocket_pos_s2_set = []
    rocket_pos_s2_set1 = []
    angular_vel_s2 = []
    total_gimbal_s2  =[]
    for item in initial_vals:

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
    rocket_pos_s2_set1.append(initial_pos)
    engine_s2.pos = initial_pos
    running = True
    points  =[]
    points.append(initial_pos_rotated)
    time = []
    eot = []
    eot_t = []
    dt = 0.1
    t = 0
    while running:

        if check == 0:
            break
        tot_prop_cg_s2 = propellant_s2.prop_cg(stage.prop_mass)
        propellant_s2.tot_prop_cg = tot_prop_cg_s2
        stage.tot_cgpos = stage.tot_cgpos1(tot_prop_cg_s2)

        radius_s2 = sqrt((rocket_pos_s2_set[-1][0])**2 + (rocket_pos_s2_set[-1][1])**2 + (rocket_pos_s2_set[-1][2])**2)

        if stage.prop_mass<= 0.75*111000:
            engine_s2.gimbal = np.array([tot_gimbal*pi/2/180,0.,0.])
        else:
           engine_s2.gimbal = np.array([0.,0.,0.])

        press_s2, temp_s2, density_s2 = atmossolver(planet.g0, planet.layers, planet.avals,radius_s2 - planet.radius, planet.r,planet.base_temp,planet.base_press)

        #Pitch, roll and yaw matrices of the engine due to the gimbal

        roll_mat1 = matrices1.roll_matrix(engine_s2.gimbal[0]) #around x axis
        pitch_mat1 = matrices1.pitch_matrix(engine_s2.gimbal[1]) #around the y axis
        yaw_mat1 = matrices1.yaw_matrix(engine_s2.gimbal[2])

        #MOI of 2nd Stage
        moi_mat_s2 = matrices1.moi_matrix(stage.prop_mass,stage.mass_empty,stage.payload_mass,stage.diameter/2,stage.height,stage.tot_cgpos)

        moment_from_gimbal_engine_s2 = engine_s2.engine2bodygimbal(stage.tot_cgpos, pitch_mat1, roll_mat1, yaw_mat1)

        angular_accel1 = stage.angular_accel_mat(moi_mat_s2,moment_from_gimbal_engine_s2) #CHECK THIS

        angular_vel_s2.append(angular_accel1*dt+angular_vel_s2[-1])
        total_gimbal_s2.append(angular_vel_s2[-1]*dt)

        tot_roll_mat1 = matrices1.roll_matrix(total_gimbal_s2[-1][0])
        tot_pitch_mat1 = matrices1.pitch_matrix(total_gimbal_s2[-1][1])
        tot_yaw_mat1 = matrices1.yaw_matrix(total_gimbal_s2[-1][2])

        #2nd Stage Thrust

        thrust3 = engine_s2.mdot*ve

        thrust3 = np.matmul(yaw_mat1,np.matmul(pitch_mat1,np.matmul(roll_mat1,thrust3)))

        stage.prop_mass = stage.prop_mass - dt * engine_s2.mdot
        if stage.prop_mass <=0.:
            engine_s2.mdot = 0.
            stage.prop_mass = 0.
            eot.append(plotter_pos_s2)

            eot_t.append(t)

        #Drag force change
        if radius_s2<planet.radius:
            drag_force_s2 = np.array([0.,0.,0.])
            rocket_vel_s2_set.append(np.array([0.,0.,0.]))
            break
        else:
            if np.linalg.norm(rocket_vel_s2_set[-1]) != 0.:
                drag_force_s2 = -0.5*density_s2*(np.linalg.norm(rocket_vel_s2_set[-1])**2)*stage.area*stage.cD*rocket_vel_s2_set[-1]/(np.linalg.norm(rocket_vel_s2_set[-1]))

            else:
                drag_force_s2 = np.array([0,0,0])

        rotated_thrust_s2 = np.matmul(tot_yaw_mat1, np.matmul(tot_pitch_mat1, np.matmul(tot_roll_mat1, thrust3)))

        total_force_s2 = (stage.mass_empty+stage.prop_mass+stage.payload_mass)*(rocket_pos_s2_set1[-1]/radius_s2)*-9.80665+rotated_thrust_s2+drag_force_s2



        #print(total_force_s2)
        rock_accel_s2 = total_force_s2/(stage.mass_empty+stage.prop_mass+stage.payload_mass)

        rocket_vel_s2_set.append(rock_accel_s2*dt+rocket_vel_s2_set[-1])

        rocket_pos_s2_set.append((np.array([465.1*cos(latitude),465.1*sin(latitude),0.])+rocket_vel_s2_set[-1])*dt + rocket_pos_s2_set[-1])

        rocket_pos_s2_set1.append(rocket_vel_s2_set[-1]*dt + rocket_pos_s2_set1[-1])

        engine_s2.pos = rocket_pos_s2_set1[-1]

        t+=dt
        plotter_pos_s2 = np.matmul(long_yaw_mat1,np.matmul(lat_pitch_mat1,np.matmul(lat_roll_mat1, rocket_pos_s2_set[-1])))
        points.append(plotter_pos_s2)
        time.append(t)

        if t > 2500 or radius_s2<planet.radius:
            running = False
            break
        t += dt

    return np.array(points),time,rocket_vel_s2_set,np.array(eot),eot_t

def moving_target(time_steps,init_pos,latitude,longitude,radius):
    dt = 0.1

    track = []
    matrices = Matrices()
    long_yaw_mat = matrices.yaw_matrix(longitude)
    lat_roll_mat = matrices.roll_matrix(latitude)
    lat_pitch_mat = matrices.pitch_matrix(-pi / 2)
    init_pos = (np.array([0.,-100000.,0.])+init_pos)/np.linalg.norm(init_pos)*radius
    tot_points = [init_pos]
    for i in range(len(time_steps)):

        new = (np.array([465.1 * cos(latitude), 465.1 * sin(latitude), 0.])) * dt + tot_points[-1]

        tot_points.append(new)
        plotter_pos = np.matmul(long_yaw_mat, np.matmul(lat_pitch_mat, np.matmul(lat_roll_mat, (tot_points[-1]))))
        plotter_pos = plotter_pos/np.linalg.norm(plotter_pos)*radius
        track.append(plotter_pos)


    return np.array(track),time_steps