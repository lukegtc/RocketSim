import numpy as np
import random
from math import *
import pickle
import cv2
timeset = np.arange(0.,2000.1,0.1)
def random_val_test(times):
    mdot_vals = []
    mdot_gimbal_vals = []
    gimbal_vals = []

    for i in range(len(timeset)):
        if i<=len(timeset)*0.001:
            mdot_vals.append(1)
            mdot_gimbal_vals.append(1)
            gimbal_vals.append(np.array([0,0,0]))
        else:
            mdot_vals.append(random.uniform(0.39,1))
            mdot_gimbal_vals.append(random.uniform(0.39,1))
            gimbal_vals.append(np.array([random.uniform(-10*pi/180,10*pi/180),random.uniform(-10*pi/180,10*pi/180),random.uniform(-10*pi/180,10*pi/180)]))

    return mdot_vals,mdot_gimbal_vals,gimbal_vals



num_episodes = 25000
throttle_pen = 3
gimbal_pen = 1
crash_penalty = 1000
#The higher the altitude the better

#The less fuel used the better

epsilon = 0.9
eps_decay = 0.9998
show_every = 2500
start_q_table = None
learning_rate = 0.1
discount = 0.95
