import numpy as np
import random
from random import randint
from class_file import *
from math import *
import pickle
import cv2
from tqdm import tqdm
import matplotlib as mpl
from rocket_pos_updater import rocket_func
import matplotlib.pyplot as plt
timeset = np.arange(0.,2000.1,0.1)
def random_val_test(timeset):
    mdot_vals = []
    mdot_gimbal_vals = []
    gimbal_vals = []

    for i in range(len(timeset)):
        if i<=len(timeset)*0.001:
            mdot_vals.append(1)
            mdot_gimbal_vals.append(1)
            gimbal_vals.append(np.array([0,0,0]))
        else:
            mdot_vals.append(1)
            mdot_gimbal_vals.append(1)
            gimbal_vals.append(np.array([0,0,0]))
            # mdot_vals.append(random.uniform(0.39,1))
            # mdot_gimbal_vals.append(random.uniform(0.39,1))
            # gimbal_vals.append(np.array([random.uniform(-10*pi/180,10*pi/180),random.uniform(-10*pi/180,10*pi/180),random.uniform(-10*pi/180,10*pi/180)]))

    return mdot_vals,mdot_gimbal_vals,gimbal_vals

def initial_set(timeset):
    first_set = []
    for i in range(len(timeset)):
        if i < 200:
            first_set.append([1,1,1,1,1,0])
        else:
            first_set.append([randint(0,2),randint(0,2),randint(0,2),randint(0,2),randint(0,2),0])

    return np.array(first_set)


def adjust(num, val,step):
    if num == 0:
        val-=step
    if num == 1:
        val = val
    if num == 2:
        val+=step
    return val


def control_val_creator(dist_from_surf,tang_rocket_vel,drag,fuel_mass,set):
    running = True
    while running:

        if len(dist_from_surf)<len(set):

            dist_from_surf = np.append(dist_from_surf,0)
            tang_rocket_vel = np.append(tang_rocket_vel,0)
            drag = np.append(drag,0)
            fuel_mass = np.append(fuel_mass,0)
        elif len(dist_from_surf) == len(set):
            running = False
        # dist_from_surf = np.array(dist_from_surf)
        # tang_rocket_vel = np.array(tang_rocket_vel)
        # drag = np.array(drag)
        # fuel_mass = np.array(fuel_mass)
    for i in tqdm(range(200,len(set))):
        if fuel_mass[i] <=0.:
            set[i,:] = 1
        else:
            set[i,:] = 0
            if dist_from_surf[i]-dist_from_surf[i-1] <= 0:
                set[i,0] +=2
                set[i,1] +=2
            else:
                if dist_from_surf[i] >= 80500.:
                    set[i,5] = 1
                else:
                    set[i,0] -=1
                    set[i,1] -=1

            if tang_rocket_vel[i]-tang_rocket_vel[i-1] <=0:
                set[i,2] +=2
                set[i,1]+=1
            else:
                set[i,2] -=1

            if drag[i]-drag[i-1] >0:
                set[i,0] -=1
            elif drag[i]-drag[i-1] <=0:
                set[i,1] +=1
        if max(set[i,:]) == 0:
            set[i] = 0
        else:
            set[i] = set[i,:]/max(set[i,:])
        for j in range(len(set[i])):
            if set[i,j]<0.33:
                set[i,j] = 0
            elif set[i,j]>0.33 and set[i,j]<0.66:
                set[i, j] = 1
            else:
                set[i,j] = 2
    return set
def alg_tester():
    set = initial_set(timeset)



    gimbal_throttle_pcts = np.array([1]*len(set))
    throttle_pcts = np.array([1]*len(set))
    gimbals = np.array([np.array([0,0,0])]*len(set))
    mpl.use('Qt5Agg')
    for i in range(2):
        print(set)
        for j in tqdm(range(len(set))):
            throttle_pcts[j] = adjust(set[j,0],throttle_pcts[j],0.01)

            if throttle_pcts[j] >1:
                throttle_pcts[j] = 1
            if throttle_pcts[j]<0.39:
                throttle_pcts[j] = 0.39
            gimbal_throttle_pcts[j] = adjust(set[j,1],gimbal_throttle_pcts[j],0.01)

            if gimbal_throttle_pcts[j]>1:
                gimbal_throttle_pcts[j]= 1
            if gimbal_throttle_pcts[j]<0.39:
                gimbal_throttle_pcts[j] = 0.39

            gimbals[j,0] = adjust(set[j,2],gimbals[j,0],-0.1*pi/180)

            if gimbals[j,0] >5*pi/180:
                gimbals[j,0] = 5*pi/180
            if gimbals[j,0] < -5*pi/180:
                gimbals[j,0] = -5*pi/180
            if set[j,5] == 1:
                throttle_pcts[j] = 0.
                gimbal_throttle_pcts[j] = 0.

        points_s1, no_fuel_points_s1, dragforce, dynamic_press, time, moments, rads, tot_mass, accel_vel_pos, stage1, stage2, earth, latitude, longitude,tang_rocket_vel,fuel_mass = rocket_func(
            Rocket(prop_mass=418700., mass_empty=27200., thrust=np.array([0, 0, 7605e3])),
            Rocket(prop_mass=111000, mass_empty=5000, height=12.6, payload_mass=22800, thrust=np.array([0., 0., 934000.])),
            Engine(mdot=266), Engine(mdot=266), (-90) * pi / 180, 0 * pi / 180,
            np.array([0, 0, 2550]), [0., 11., 20., 32., 47., 51., 71., 86.], [-6.5, 0, 1, 2.8, 0, -2.8, -2.], 288.15,
            101325, 287, 9.80665, gimbal_throttle_pcts, throttle_pcts, timeset, gimbals)
     #   print(np.linalg.norm(tang_rocket_vel))
        set = control_val_creator(rads,tang_rocket_vel[:,0],dragforce,fuel_mass,set)
        print(i)
    plt.plot(timeset,set)
    plt.show()

    plt.plot(points_s1[:,0],points_s1[:,1])
    plt.show()

    plt.plot(timeset,throttle_pcts)
    plt.show()

    plt.plot(time,rads)
    plt.show()