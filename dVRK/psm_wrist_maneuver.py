from pyscurve import ScurvePlanner, plot_trajectory
import rospy
import dvrk
import math
import numpy as np
import matplotlib.pyplot as plt
import time

p = dvrk.psm("PSM2")
p.home()
Ts = 0.001;
# The units are SI (rad for joints 0,1,3,4,5 and m for joint 2)
# move all the axes to the start position defined below
start_pos = [3.4/180*math.pi,1.1/180*math.pi,0.08,0,0,0]
axis = 0
for item in start_pos:
    p.move_joint_one(float(item),axis,interpolate = True, blocking = True)
    axis+=1

p.move_jaw(0,interpolate = False, blocking = False)
cycleNum = 1

axisIndex = [3,4,5];
motion_rng = [[-math.pi/2,math.pi/2,0],[-math.pi/3,math.pi/3,0],[-math.pi/3,math.pi/3,0]]

v_max = 0.5
a_max = 5
j_max = 100

for counter in range(len(axisIndex)):
    for cycle_counter in range(cycleNum):
        for value in motion_rng[counter]:
            q0 = [p.get_current_joint_position()[axisIndex[counter]]]
            q1 = [value]
            v0 = [0]
            v1 = [0]

            q = ScurvePlanner()
            tr = q.plan_trajectory(q0, q1, v0, v1, v_max, a_max, j_max)

            dof = tr.dof
            timesteps = int(max(tr.time) / Ts)
            time_bw = np.linspace(0, max(tr.time), timesteps)

            # NOW
            # profiles[t]           --- profiles for each DOF at time x[t]
            # profiles[t][d]        --- profile for d DOF at time x[t]
            # profiles[t][d][k]     --- accel/vel/pos profile for d DOF at time x[t]
            p_list = [tr(t) for t in time_bw]
            profiles = np.asarray(p_list)
            ls_cmd = profiles[:,0,2]
            for item in np.nditer(ls_cmd):
                p.move_joint_one(float(item),axisIndex[counter],interpolate = False, blocking = False)
                time.sleep(Ts)

jaw_pos = [0,math.pi/3,0,-10/180*math.pi,0]
for value in jaw_pos:
    q0 = [p.get_current_jaw_position()]
    q1 = [value]
    v0 = [0]
    v1 = [0]

    q = ScurvePlanner()
    tr = q.plan_trajectory(q0, q1, v0, v1, v_max, a_max, j_max)

    dof = tr.dof
    timesteps = int(max(tr.time) / Ts)
    time_bw = np.linspace(0, max(tr.time), timesteps)

    # NOW
    # profiles[t]           --- profiles for each DOF at time x[t]
    # profiles[t][d]        --- profile for d DOF at time x[t]
    # profiles[t][d][k]     --- accel/vel/pos profile for d DOF at time x[t]
    p_list = [tr(t) for t in time_bw]
    profiles = np.asarray(p_list)
    ls_cmd = profiles[:,0,2]
    for item in np.nditer(ls_cmd):
        p.move_jaw(float(item),interpolate = False, blocking = False)
        time.sleep(Ts) 

axis = 0
for item in start_pos:
    p.move_joint_one(float(item),axis,interpolate = True, blocking = True)
    axis+=1