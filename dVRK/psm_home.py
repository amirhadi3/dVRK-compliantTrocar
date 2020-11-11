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

p.close_jaw()
