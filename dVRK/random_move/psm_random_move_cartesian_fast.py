from numpy.random import randint
import matplotlib.pyplot as plt
import rospy
import dvrk
import math
import random
import numpy as np
import time
from std_msgs.msg import Float64,String
import PyKDL
from pyscurve import ScurvePlanner, plot_trajectory

class psm_random_move_fast(object):

    def __init__(self,psm_name, totaltime_s,num_wayPoints,joints,joint_limits, Ts):
        # create the psm
	self.statSubscriber = rospy.Subscriber(psm_name+'/recState',String,self.__subCallback)
        self.recState = ''
        self.arm = dvrk.psm(psm_name)       
        self.time = totaltime_s
        self.num_wayPoints = num_wayPoints
        self.joints = joints
        self.joint_limits = joint_limits
	self.Ts = Ts

    def __subCallback(self,msg):
	self.recState = msg.data

    def home(self):
        self.arm.home()
        # move the arm to a home position
        # The units are SI (rad for joints 0,1,3,4,5 and m for joint 2)
        # move all the axes to the start position defined below
        self.home_pos = [3.4/180*math.pi,1.1/180*math.pi,0.075,0,0,0]
        axis = 0
        for item in self.home_pos:
            self.arm.move_joint_one(float(item),axis,interpolate = True, blocking = True)
            axis+=1
	jawPos = self.arm.get_desired_jaw_position()
	self.arm.move_jaw(jawPos,interpolate = True, blocking = True)

    def localhome(self):
        # move the arm to a home position
        # The units are SI (rad for joints 0,1,3,4,5 and m for joint 2)
        # move all the axes to the start position defined below
	ins = self.arm.get_current_joint_position()[2]
        self.home_pos = [3.4/180*math.pi,1.1/180*math.pi,ins,0,0,0]
        axis = 0
        for item in self.home_pos:
            self.arm.move_joint_one(float(item),axis,interpolate = True, blocking = True)
            axis+=1
	self.arm.move_jaw(math.pi*3/180,interpolate = True, blocking = True)
    
    def __piecewiseSineInterpolate(self,x0,xe,dur,Ts):
        a = (xe+x0)/2
        b = (xe-x0)/2

        t = np.arange(0,dur,Ts)+Ts
        x = (a+b*np.sin(-math.pi/2+math.pi/dur*t))
        return x

    def __genWayPoints(self,x0,xe,numWayPoints,limits):
        locmin = limits[0]
	locmax = limits[1]
        p = locmin+(np.random.rand(numWayPoints+4))*(locmax-locmin)
        p[0] = x0
	sgn = np.ones(np.shape(p))
	sgn[1::2]=np.sign(p[1])
	sgn[2::2]=-np.sign(p[1])
	p = np.abs(p)*sgn
	p = np.append(p,xe)
	return p

    def __findTimeIntervals(self,p,t):
	num = int(t/self.Ts);
        dispVec = np.abs(np.diff(p))	
	totdist = np.cumsum(dispVec)
	timeVec = np.floor(dispVec/totdist[-1]*t/self.Ts);
	while(np.cumsum(timeVec)[-1]<num):
	    intervalIndex = np.random.randint(0,len(timeVec))
	    timeVec[intervalIndex]+=1
	while(np.cumsum(timeVec)[-1]>num):
	    intervalIndex = np.random.randint(0,len(timeVec))
	    timeVec[intervalIndex]-=1
	#print(np.cumsum(timeVec)[-1])
	return timeVec*self.Ts

    def __singleAxisInterp(self,wayPoints,dur):
	joint_traj = np.array([0])
	for index in range(1,len(wayPoints)):
		joint_traj = np.append(joint_traj,
				self.__piecewiseSineInterpolate(joint_traj[-1],
				wayPoints[index],dur[index-1],
				self.Ts)
				)
	return joint_traj

    def singleAxisTrajGen(self,xs,xe,numWayPoints,jointLimits):
	wayPoints = self.__genWayPoints(xs,xe,numWayPoints,jointLimits)
        dur = self.__findTimeIntervals(wayPoints,self.time)
	t = self.__singleAxisInterp(wayPoints,dur)
	return t

    def traj_generate(self):
	counter = 0
	for joint in self.joints:
 	    joint_traj = self.singleAxisTrajGen(0,0,self.num_wayPoints,self.joint_limits[counter])
	    numEl = len(joint_traj)
	    if counter == 0:
		traj = np.zeros((len(joint_traj),len(self.joints)))
		numRow = np.shape(traj)[0]
            
	    if numEl > numRow:
		traj[:,counter] = joint_traj[range(numRow)]
	    elif numEl < numRow:
            	traj[range(numEl),counter] = joint_traj
		padList = -(np.arange(numRow-numEl)+1)		
		traj[padList,counter]=joint_traj[-1]
	    else:
		traj[:,counter] = joint_traj
	    #plt.plot(traj[:,joint])
	    #plt.show()
	    counter+=1

	return traj

    def insertion(self):
	q0 = [0]
	q1 = [-0.102]
	v0 = [0]
	v1 = [0]

	v_max = 0.003
	a_max = 0.5
	j_max = 1.00

	q = ScurvePlanner()
	tr = q.plan_trajectory(q0, q1, v0, v1, v_max, a_max, j_max)

	dof = tr.dof
	timesteps = int(max(tr.time) / self.Ts)
	time_fw = np.linspace(0, max(tr.time), timesteps)

	# NOW
	# profiles[t]           --- profiles for each DOF at time x[t]
	# profiles[t][d]        --- profile for d DOF at time x[t]
	# profiles[t][d][k]     --- accel/vel/pos profile for d DOF at time x[t]
	p_list = [tr(t) for t in time_fw]
	profiles = np.asarray(p_list)
	self.z = profiles[:,0,2]
	self.z = np.append(self.z,np.flip(self.z))
	self.time = 2*time_fw[-1]

if __name__ == "__main__":
	armname = 'PSM2';
	rospy.init_node(str(armname+'_insertionPublisher'))
	pub = rospy.Publisher('PSM2'+'/insertion',Float64, latch = True, queue_size = 1)

	#axisIndex = np.array([2,3,4,5],dtype=np.int64);
	#motion_rng = [[-0.03,0.03],[-math.pi/3,math.pi/3],[-math.pi/4,math.pi/4],[-math.pi/4,math.pi/4]]
	axisIndex = [0,1,2,3,4,5]
	motion_rng = [[-0.02,0.02],[-0.02,0.02],[-0.02,0.02],[-math.pi/3.2,math.pi/3.2],[-math.pi/3.2,math.pi/3.2],[-math.pi/3.2,math.pi/3.2]]
	p = psm_random_move_fast(armname,45,16,axisIndex,motion_rng,0.001)
	p.home()
	p.insertion()

	if True:
	    plt.plot(p.z)
	    plt.show()
	    
            #for axis in range(len(p.joints)):
	    #    plt.plot(traj[:,axis])
	    #    plt.show()
	for i in range(1):
            curPos = p.arm.get_current_position()
	    jp = p.arm.get_desired_jaw_position()
            jgrip = p.singleAxisTrajGen(0,0,p.num_wayPoints,[-math.pi*5/180,math.pi*5/180]) + jp
	    traj = p.traj_generate()
	    dpoints = min(len(p.z),len(traj[:,1]),len(jgrip))
	    pub.publish(0.075)
            #while p.recState != 'Go':
            #    pass
	    if i==0:
	        time.sleep(10)
	    for index in range(dpoints):
		zaber_axis = 0.075-p.z[index]
		pub.publish(zaber_axis)
		print('insertion = ',zaber_axis)
		desPos = PyKDL.Frame(curPos)
		desPos.M.DoRotZ(traj[index,5])
		desPos.M.DoRotX(traj[index,3])
		desPos.M.DoRotY(traj[index,4])
		desPos.p[0] = desPos.p[0]+traj[index,0]
		desPos.p[1] = desPos.p[1]+traj[index,1]
		desPos.p[2] = desPos.p[2]+p.z[index]+traj[index,2]
		p.arm.move(desPos,interpolate = False, blocking = False)
		p.arm.move_jaw(jgrip[index],interpolate = False, blocking = False)
        	time.sleep(p.Ts)
	    p.home()
	time.sleep(3)
	pub.publish(0)
	    # generate jaw motion
 	    


