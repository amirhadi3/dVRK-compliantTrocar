from pyscurve import ScurvePlanner, plot_trajectory
import rospy
import dvrk
import math
import numpy as np
import matplotlib.pyplot as plt
import time

p = dvrk.psm("PSM2")
p.home()
# The units are SI (rad for joints 0,1,3,4,5 and m for joint 2)
# move all the axes to the start position defined below
start_pos = [3.4/180*math.pi,1.1/180*math.pi,0.08,0,0,0]
axis = 0
for item in start_pos:
    p.move_joint_one(float(item),axis,interpolate = True, blocking = True)
    axis+=1

# ---- Motion cycles and duration: start ----
num_of_cycles = 10                  # number of cycles
duration =50                        # duration (sec)
# ---- Motion cycles and duration: end ----

# ---- construct a sinusoidal position command: start ----
Ts = 0.01 #s
ls_min = 0.08
ls_max = 0.19

# Forward motion #
q0 = [ls_min]
q1 = [ls_max]
v0 = [0]
v1 = [0]

v_max = 0.020
a_max = 0.015
j_max = 0.0100

q = ScurvePlanner()
tr = q.plan_trajectory(q0, q1, v0, v1, v_max, a_max, j_max)

dof = tr.dof
timesteps = int(max(tr.time) / Ts)
time_fw = np.linspace(0, max(tr.time), timesteps)

# NOW
# profiles[t]           --- profiles for each DOF at time x[t]
# profiles[t][d]        --- profile for d DOF at time x[t]
# profiles[t][d][k]     --- accel/vel/pos profile for d DOF at time x[t]
p_list = [tr(t) for t in time_fw]
profiles = np.asarray(p_list)
ls_cmd = profiles[:,0,2]


# Backward motion #

q0 = [ls_cmd[-1]]
q1 = [ls_min]
v0 = [0]
v1 = [0]

v_max = 0.050
a_max = 1.5
j_max = 1.00

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
ls_cmd = np.append(ls_cmd,profiles[:,0,2])
time_array = np.append(time_fw,time_fw[-1]+time_bw)
plot = False

if plot:
    plt.plot(time_array,ls_cmd)
    plt.xlabel('time (s)')
    plt.ylabel('disp. (m)')
    plt.title('Instrument displacement inside the trocar')
    plt.grid(True)
    plt.show()
# ---- construct a sinusoidal position command: end ----

for item in np.nditer(ls_cmd):
    p.move_joint_one(float(item),2,interpolate = False, blocking = False)
    time.sleep(Ts)