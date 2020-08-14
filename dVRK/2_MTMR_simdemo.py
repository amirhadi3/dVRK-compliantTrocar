import rospy
import math
import numpy as np
import matplotlib.pyplot as plt
# I accidentally came up with two solutions to move the robot arm
# Initially I started with sending motion commands using the dVRK
# builtin function move_joint_one. I then tried publishing directly 
# to the set_joint_position topic. Both approaches successfuly move 
# the robot arm.
approach = 1

# ---- Motion cycles and duration: start ----
num_of_cycles = 10					# number of cycles
duration =10 						# duration (sec)
# ---- Motion cycles and duration: end ----

# ---- construct a sinusoidal position command: start ----
Ts = 0.01 #s
pos_cmd_length = duration // Ts

time_array = np.arange(pos_cmd_length)
ls_min = -math.pi/6
ls_max = math.pi/6

ls_period = pos_cmd_length / num_of_cycles

f = lambda x: ls_min+(ls_max-ls_min)*(math.sin(x/2/ls_period*(2*math.pi))**2)

ls_cmd = np.vectorize(f)(time_array)

plot = True

if plot:
	plt.plot(time_array,ls_cmd)
	plt.xlabel('time (s)')
	plt.ylabel('disp. (m)')
	plt.title('Instrument displacement inside the trocar')
	plt.grid(True)
	plt.show()
# ---- construct a sinusoidal position command: end ----

# Approach 1: 


if approach == 1:
	import dvrk
	import time
	m = dvrk.mtm("MTMR")
	m.home()
	m.move_joint_one(ls_cmd[0],6)
	for item in np.nditer(ls_cmd):
		m.move_joint_one(float(item),6,interpolate = False, blocking = False)
		time.sleep(Ts)

else:
# Approach 2:
	from sensor_msgs.msg import JointState

	rospy.init_node('mtm_Cmd_publisher')
	pub = rospy.Publisher('/dvrk/MTMR/set_position_joint',JointState,queue_size = 1)

	rate = rospy.Rate(1/Ts)
	cmd = JointState()
	cmd.position = [0] * 7

	for item in np.nditer(ls_cmd):
		pub.publish(cmd)
		cmd.position[0] = float(item)
		rate.sleep()
