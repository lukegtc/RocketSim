import numpy as np
from math import *
def rocket_func(stage1,stage2,propellant1,gimbal_engine1,engine1,engine2,s1thrust:np.array,s2thrust:np.array,latitude,longitude,exit_vel,
                s1engine_pos,s2engine_pos,rocket_pos,payload,planet):
    dt = 0.1
    t = 0.
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
#    skin_temps = []
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
    rocket_poss.append(rocket_pos)

    skin_temp = 0
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
        press,temp,density=atmossolver(planet.g0,planet.layers,planet.avals,radius-planet.radius,planet.r,planet.base_temp,planet.base_press)



        if stage1.prop_mass<= 0.95*418700:
            gimbal_engine1.gimbal = np.array([-0*pi/2/180,0,0])

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

        angular_vel_mat.append( angular_accel*dt+angular_vel_mat[-1])
        total_gimbal_mat.append(angular_vel_mat[-1]*dt + total_gimbal_mat[-1])

        tot_roll_mat = matrices.roll_matrix(total_gimbal_mat[-1][0])
        tot_pitch_mat = matrices.pitch_matrix(total_gimbal_mat[-1][1])
        tot_yaw_mat = matrices.yaw_matrix(total_gimbal_mat[-1][2])


        #1st Stage Thrust w gimbal in engine1------------------------
        thrust1 = gimbal_engine1.mdot * exit_vel
        thrust1 = np.matmul(yaw_mat,np.matmul(pitch_mat,np.matmul(roll_mat,thrust1)))
        thrust2 = engine2.mdot*exit_vel



        stage1.prop_mass = stage1.prop_mass - dt * gimbal_engine1.mdot - dt * engine2.mdot * 8


        if stage1.prop_mass <=0.:


            check = 1
            gimbal_engine1.mdot = 0
            engine2.mdot = 0
            stage1.prop_mass = 0
            #2nd stage and payload detach from stage 1
            stage1.payload_mass = 0
        else:
            check = 0
        checks.append(check)
        if radius<planet.radius:
            drag_force = np.array([0.,0.,0.])
            rocket_vels.append(np.array([0.,0.,0.]))
        else:
            if np.linalg.norm(rocket_vels[-1]) != 0.:
                drag_force = -0.5*density*(np.linalg.norm(rocket_vels[-1])**2)*stage1.area*stage1.cD*rocket_vels[-1]/(np.linalg.norm(rocket_vels[-1]))
                #fix something with the very last bit of the drag equation because its giving the wrong drag vectors i think invert the lat long stuff

            else:
                drag_force = np.array([0,0,0])

        #summation of all forces acting on the rocket

        rotated_thrust = np.matmul(tot_yaw_mat,np.matmul(tot_pitch_mat,np.matmul(tot_roll_mat,thrust2*8+thrust1)))

        total_force =(stage1.mass_empty+stage1.prop_mass+stage1.payload_mass)*(rocket_poss[-1]/radius)*-9.80665+rotated_thrust+drag_force #Fix coordinate stuff


        rocket_accel = total_force/(stage1.mass_empty+stage1.prop_mass+stage1.payload_mass)
        rocket_vels.append(rocket_accel*dt+rocket_vels[-1])

        rocket_poss.append((np.array([465.1*cos(latitude),465.1*sin(latitude),0.])+rocket_vels[-1])*dt + rocket_poss[-1])

        gimbal_engine1.pos = rocket_poss[-1]
        temps.append(temp)

        accel_vel_pos.append([rocket_accel,rocket_vels[-1],rocket_poss[-1],angular_vel_mat[-1],total_gimbal_mat[-1],angular_accel,check])




        #change in temp due to atmosphere
       # nose_radius = 0.7 #don't know if this is true but eh close enough for F9
        #skin_temp = friction_temp(rocket_vels[-1],density,nose_radius,rocket1.specific_heat,rocket1.payload_mass+rocket1.mass_empty+rocket1.prop_mass,radius,earth.radius,earth.layers,skin_temp)
        #tot_skin_temp = temp + skin_temp


        if skin_temp >= 1923.15:  #This reentry temp is totally made up based on space shuttle values
            running = False
            print('Burnt up on reentry or ascent')

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

        if t > 2500 or np.linalg.norm(rocket_poss[-1])<planet.radius:
            running = False
            break
        t += dt
    return xpos,ypos,zpos,xpos1,ypos1,zpos1,dragforce,time,moments,rads,rads,tot_mass