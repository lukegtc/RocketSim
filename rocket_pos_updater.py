from class_file import *
import numpy as np
import matplotlib.pyplot as plt
from Atmosphere import *
from math import *
import matplotlib as mpl
def rocket_func(stage1,stage2,gimbal_engine1,engine1,latitude,longitude,total_gimbal,exit_vel,layers,avals,base_temp,base_press,r,g0):
    GM = 3.986004418e14
    earth = Planet(6371e3,g0,layers,avals,base_temp,base_press,r,GM)
    gimbal_engine1.pos = np.array([0.,0.,0. + earth.radius])

    dt = 0.1
    t = 0.
    gimbal_engine1.thrust = exit_vel * engine1.mdot
    engine1.thrust = exit_vel * engine1.mdot

    propellant1 = Propellant(mass=stage1.prop_mass,radius=stage1.radius,height=stage1.height)
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
    presss = []
    vels = []
    cg_vals = []
    moments= []
    prop_mass = []

    tot_mass = []

    #2nd stage values-------


    xpos1 = []
    ypos1 = []
    zpos1 = []

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

    angular_vel_mat.append(np.array([0.,0.,0.]))
    total_gimbal_mat.append(np.array([0.,0.,0.]))

    while running:

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
        #print(radius)
        press,temp,density=atmossolver(earth.g0,earth.layers,earth.avals,radius-earth.radius,earth.r,earth.base_temp,earth.base_press)

        #Gimballing rules

        if stage1.prop_mass<= 0.95*418700:
            gimbal_engine1.gimbal = np.array([-total_gimbal*pi/2/180,0,0])

        else:
            gimbal_engine1.gimbal[0] = 0



        #Stage 1 engine pitch roll yaw matrices

        matrices = Matrices()
        roll_mat = matrices.roll_matrix(gimbal_engine1.gimbal[0]) #around x axis
        pitch_mat = matrices.pitch_matrix(gimbal_engine1.gimbal[1]) #around the y axis
        yaw_mat = matrices.yaw_matrix(gimbal_engine1.gimbal[2]) #around the z axis



        #moment of inertia------------------------------
        moi_mat = matrices.moi_matrix(stage1.prop_mass,stage1.mass_empty,stage1.payload_mass,stage1.diameter/2,stage1.height,stage1.tot_cgpos)

        moment_from_gimbal_engine = gimbal_engine1.engine2bodygimbal(stage1.tot_cgpos, pitch_mat, roll_mat, yaw_mat)

        angular_accel = stage1.angular_accel_mat(moi_mat,moment_from_gimbal_engine) #CHECK THIS

        angular_vel_mat.append(angular_accel*dt+angular_vel_mat[-1])
        total_gimbal_mat.append(angular_vel_mat[-1]*dt + total_gimbal_mat[-1])

        tot_roll_mat = matrices.roll_matrix(total_gimbal_mat[-1][0])
        tot_pitch_mat = matrices.pitch_matrix(total_gimbal_mat[-1][1])
        tot_yaw_mat = matrices.yaw_matrix(total_gimbal_mat[-1][2])


        #1st Stage Thrust w gimbal in engine1------------------------
        thrust1 = gimbal_engine1.mdot * exit_vel
        thrust1 = np.matmul(yaw_mat,np.matmul(pitch_mat,np.matmul(roll_mat,thrust1)))
        thrust2 = engine1.mdot*exit_vel



        stage1.prop_mass = stage1.prop_mass - dt * gimbal_engine1.mdot - dt * engine1.mdot * 8


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

        total_force =(stage1.mass_empty+stage1.prop_mass+stage1.payload_mass)*(rocket_poss[-1]/radius)*-9.80665+rotated_thrust+drag_force #Fix coordinate stuff


        rocket_accel = total_force/(stage1.mass_empty+stage1.prop_mass+stage1.payload_mass)
        rocket_vels.append(rocket_accel*dt+rocket_vels[-1])

        rocket_poss.append((np.array([465.1*cos(latitude),465.1*sin(latitude),0.])+rocket_vels[-1])*dt + rocket_poss[-1])

        gimbal_engine1.pos = rocket_poss[-1]
        temps.append(temp)

        accel_vel_pos.append([rocket_accel,rocket_vels[-1],rocket_poss[-1],angular_vel_mat[-1],total_gimbal_mat[-1],angular_accel,check])

        long_yaw_mat = matrices.yaw_matrix(longitude)
        lat_roll_mat = matrices.roll_matrix(latitude)
        lat_pitch_mat = matrices.pitch_matrix(-pi/2)
        plotter_pos = np.matmul(long_yaw_mat,np.matmul(lat_pitch_mat,np.matmul(lat_roll_mat, (gimbal_engine1.pos)))) #
        points.append(plotter_pos)


        if 0.<stage1.prop_mass:

            s2drag.append(np.linalg.norm(drag_force))
            zpos1.append(plotter_pos[2])
            ypos1.append(plotter_pos[1])
            xpos1.append(plotter_pos[0])

        #plotting section----------------------------------------------------------------------------------------------
        zpos.append(plotter_pos[2])
        ypos.append(plotter_pos[1])
        xpos.append(plotter_pos[0])
        prop_mass.append(stage1.prop_mass)
        ttlforce.append((np.linalg.norm(total_force)))

        dragforce.append(np.linalg.norm(drag_force))
        time1.append(t)
        time.append(t)
        moments.append(moment_from_gimbal_engine[0])
        densitys.append(density)
        tot_mass.append(stage1.mass_empty+stage1.payload_mass+stage1.prop_mass)
        temps.append(temp)
        presss.append(press)
        rads.append(radius-6371e3)
        vels.append(np.linalg.norm(rocket_vels[-1]))
        cg_vals.append(stage1.tot_cgpos[2])

        if t > 2500 or np.linalg.norm(rocket_poss[-1])<earth.radius:
            running = False
            break
        t += dt
    return xpos,ypos,zpos,xpos1,ypos1,zpos1,dragforce,dynamic_press,time,moments,rads,rads,tot_mass,accel_vel_pos,stage1,stage2,earth,latitude,longitude



def staging_wo_1st_stage(stage,exit_vel,initial_vals,latitude,longitude,planet):

    matrices1 = Matrices()
    propellant_s2 = Propellant(mass=stage.prop_mass,radius=stage.radius,height=stage.height)
    ve = exit_vel
    engine3 = Engine(mdot=934)
    engine3.gimbal = np.array([0.,0.,0.])
    engine3.thrust = ve*engine3.mdot
    stage.thrust = engine3.thrust

    rocket_vel_s2_set = []
    rocket_pos_s2_set = []
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
    engine3.pos = initial_pos
    running = True
    points  =[]
    points.append(initial_pos_rotated)
    time = []
    dt = 0.1
    t = 0
    while running:
        if check == 0:
            break
        tot_prop_cg_s2 = propellant_s2.prop_cg(stage.prop_mass)
        propellant_s2.tot_prop_cg = tot_prop_cg_s2
        stage.tot_cgpos = stage.tot_cgpos1(tot_prop_cg_s2)

        radius_s2 = sqrt((rocket_pos_s2_set[-1][0])**2 + (rocket_pos_s2_set[-1][1])**2 + (rocket_pos_s2_set[-1][2])**2)


        if stage.prop_mass<= 0.99*111000:
            engine3.gimbal = np.array([0*pi/2/180,0.,0.])
        else:
           engine3.gimbal = np.array([0.,0.,0.])


        press_s2, temp_s2, density_s2 = atmossolver(planet.g0, planet.layers, planet.avals,radius_s2 - planet.radius, planet.r,planet.base_temp,planet.base_press)

        #Pitch, roll and yaw matrices of the engine due to the gimbal

        roll_mat1 = matrices1.roll_matrix(engine3.gimbal[0]) #around x axis
        pitch_mat1 = matrices1.pitch_matrix(engine3.gimbal[1]) #around the y axis
        yaw_mat1 = matrices1.yaw_matrix(engine3.gimbal[2])
        #MOI of 2nd Stage
        moi_mat_s2 = matrices1.moi_matrix(stage.prop_mass,stage.mass_empty,stage.payload_mass,stage.diameter/2,stage.height,stage.tot_cgpos)

        moment_from_gimbal_engine_s2 = engine3.engine2bodygimbal(stage.tot_cgpos, pitch_mat1, roll_mat1, yaw_mat1)

        angular_accel1 = stage.angular_accel_mat(moi_mat_s2,moment_from_gimbal_engine_s2) #CHECK THIS

        angular_vel_s2.append(angular_accel1*dt+angular_vel_s2[-1])
        total_gimbal_s2.append(angular_vel_s2[-1]*dt)

        tot_roll_mat1 = matrices1.roll_matrix(total_gimbal_s2[-1][0])
        tot_pitch_mat1 = matrices1.pitch_matrix(total_gimbal_s2[-1][1])
        tot_yaw_mat1 = matrices1.yaw_matrix(total_gimbal_s2[-1][2])


        #2nd Stage Thrust

        thrust3 = engine3.mdot*ve

        thrust3 = np.matmul(yaw_mat1,np.matmul(pitch_mat1,np.matmul(roll_mat1,thrust3)))



        stage.prop_mass = stage.prop_mass - dt * engine3.mdot
        if stage.prop_mass <=0.:
            engine3.mdot = 0.
            stage.prop_mass = 0.


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

        total_force_s2 = (stage.mass_empty+stage.prop_mass+stage.payload_mass)*(rocket_pos_s2_set[-1]/radius_s2)*-9.80665+rotated_thrust_s2+drag_force_s2

        rock_accel_s2 = total_force_s2/(stage.mass_empty+stage.prop_mass+stage.payload_mass)

        rocket_vel_s2_set.append(rock_accel_s2*dt+rocket_vel_s2_set[-1])

        rocket_pos_s2_set.append((np.array([465.1*cos(latitude),465.1*sin(latitude),0.])+rocket_vel_s2_set[-1])*dt + rocket_pos_s2_set[-1])
        engine3.pos = rocket_pos_s2_set[-1]

        t+=dt
        plotter_pos_s2 = np.matmul(long_yaw_mat1,np.matmul(lat_pitch_mat1,np.matmul(lat_roll_mat1, rocket_pos_s2_set[-1])))
        points.append(plotter_pos_s2)
        time.append(t)
        xpos_s2.append(plotter_pos_s2[0])
        ypos_s2.append(plotter_pos_s2[1])
        zpos_s2.append(plotter_pos_s2[2])


        if t > 2500 or radius_s2<planet.radius:
            running = False
            break
        t += dt

    return points,time,rocket_vel_s2_set,xpos_s2,ypos_s2,zpos_s2


