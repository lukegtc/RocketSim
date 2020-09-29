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

def initial_set(timeset,pop):

    first_pop = []
    for i in range(pop):
        first_set = []
        for i in range(len(timeset)):
            if i<200:
                first_set.append([1,1,1,1,1,0])
            else:

                first_set.append([randint(0,2),randint(0,2),randint(0,2),randint(0,2),randint(0,2),0])
        first_pop.append(first_set)
    return np.array(first_pop)


def adjust(num, val,step):
    if num == 0:
        val-=step
    if num == 1:
        val = val
    if num == 2:
        val+=step
    return val

def fitness_num(dist_from_surf,tang_rocket_vel,drag,fuel_mass):
    #set is a group of mdot_vals,mdot_gimbal_vals,gimbal_vals
    dist_fit = max(dist_from_surf)
    tang_vel_fit = max(tang_rocket_vel)
    for i in range(len(drag)):

        drag_fit = max(drag)
        if max(dist_from_surf)>81000:
            if dist_from_surf[i] >= 78000 and dist_from_surf[i]<=81000:
                fuel_mass_fit = min(fuel_mass)
        else:
            fuel_mass_fit = min(fuel_mass)
    return -dist_fit-tang_vel_fit+drag_fit-fuel_mass_fit

def new_pop(old_pop1,old_pop2,pop_size):
    new_popu = []
    set_check = []
    for j in range(len(old_pop1)):
        count = 0
        for k in range(len(old_pop1[j])):
            if old_pop1[j,k] == old_pop2[j,k]:
                count+=1
        if count == len(old_pop2[j]):
            set_check.append(1)
        else:
            set_check.append(0)
    for i in range(pop_size):
        new_unit = []
        for k in range(len(set_check)):
            if set_check[k] == 1:
                new_unit.append(old_pop1[k])
            elif set_check[k] == 0:
                new_unit.append([randint(0,2),randint(0,2),randint(0,2),randint(0,2),randint(0,2),0])
        new_popu.append(new_unit)


    return np.array(new_popu)
set = initial_set(timeset,100)


print('done')

mpl.use('Qt5Agg')
running = True

final = 0
while running:
    fit_nums = []
    for i in range(len(set)):
        print(final,i)
        gimbal_throttle_pcts = np.array([1] * len(timeset))
        throttle_pcts = np.array([1] * len(timeset))
        gimbals = np.array([np.array([0, 0, 0])] * len(timeset))
        for j in tqdm(range(len(timeset))):

            throttle_pcts[j] = adjust(set[i,j,0],throttle_pcts[j],0.01)

            if throttle_pcts[j] >1:
                throttle_pcts[j] = 1
            if throttle_pcts[j]<0.39:
                throttle_pcts[j] = 0.39
            gimbal_throttle_pcts[j] = adjust(set[i,j,1],gimbal_throttle_pcts[j],0.01)

            if gimbal_throttle_pcts[j]>1:
                gimbal_throttle_pcts[j]= 1
            if gimbal_throttle_pcts[j]<0.39:
                gimbal_throttle_pcts[j] = 0.39

            gimbals[j,0] = adjust(set[i,j,2],gimbals[j,0],-0.1*pi/180)

            if gimbals[j,0] >5*pi/180:
                gimbals[j,0] = 5*pi/180
            if gimbals[j,0] < -5*pi/180:
                gimbals[j,0] = -5*pi/180
            if set[i,j,5] == 1:
                throttle_pcts[j] = 0.
                gimbal_throttle_pcts[j] = 0.

        points_s1, no_fuel_points_s1, dragforce, dynamic_press, time, moments, rads, tot_mass, accel_vel_pos, stage1, stage2, earth, latitude, longitude,tang_rocket_vel,fuel_mass = rocket_func(
            Rocket(prop_mass=418700., mass_empty=27200., thrust=np.array([0, 0, 7605e3])),
            Rocket(prop_mass=111000, mass_empty=5000, height=12.6, payload_mass=22800, thrust=np.array([0., 0., 934000.])),
            Engine(mdot=266), Engine(mdot=266), (-90) * pi / 180, 0 * pi / 180,
            np.array([0, 0, 2550]), [0., 11., 20., 32., 47., 51., 71., 86.], [-6.5, 0, 1, 2.8, 0, -2.8, -2.], 288.15,
            101325, 287, 9.80665, gimbal_throttle_pcts, throttle_pcts, timeset, gimbals)
        fit_num = fitness_num(rads,tang_rocket_vel,dragforce,fuel_mass)
        fit_nums.append(fit_num)
        print(points_s1)
    best1 = fit_nums.index(min(fit_nums))
    fit_nums.remove(fit_nums[best1])

    best2 = fit_nums.index(min(fit_nums))

    set = new_pop(set[best1],set[best2],100)

    final+=1
    if final == 100:
        running = False


# plt.plot(timeset,set)
# plt.show()

plt.plot(points_s1[:,0],points_s1[:,1])
plt.show()

plt.plot(timeset,throttle_pcts)
plt.show()

plt.plot(time,rads)
plt.show()