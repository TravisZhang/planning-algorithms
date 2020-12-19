import numpy as np
import math


START_LINE = [[0, 3], [0, 4], [0, 5], [0, 6]]
FINISH_LINE = [[34, 11], [33, 11], [32, 11]]
# we can apply 1 m/s^2 acc as input to current vx/vy
ACTION_SPACE = [[1, 1], [0, 1], [1, 0], [0, 0], [-1, 0], [0, -1], [1, -1], [-1, 1], [-1, -1]]
action_assert_list = [-1, 0, 1]
FINISH = 3
START = 2
FREE = 0
OCCUPIED = 1
OUTBOUND = -1

race_track = np.array([
    [1, 1, 1, 2, 2, 2, 2, 1, 1, 1, 1, 1],
    [1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [1, 1, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [1, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3],
    [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3]
], dtype=np.int32)



class VehicleActionSpace:
    KtoRad = 1.0/180.0*math.pi
    KtoDeg = 1.0/math.pi*180.0

    curv_spd_table = [
    [0.005, 150.0], 
    [0.0075, 130.0], 
    [0.001, 120.0], 
    [0.00143, 110.0], 
    [0.00167, 100.0], 
    [0.002, 90.0], 
    [0.0033, 80.0], 
    [0.004, 70.0], 
    [0.005, 60.0],
    [0.0067, 55.0],
    [0.01, 47.0],
    [0.0167, 38.0],
    [0.02, 35.0],
    [0.0333, 25.0],
    [0.2, 10.0]]

    # vehicle params & constraints
    axis_dist = 2.5
    max_acc = 2.0
    max_spd = 3.0
    steer_ratio = 17.0
    max_steering_angle = 480.0/steer_ratio*KtoRad
    max_steering_spd = 30.0/steer_ratio*KtoRad
    max_steering_spd_low_spd = 720.0/steer_ratio*KtoRad

    def __init__(self, v, delta, dt, sample_num):
        max_curv = self.curv_spd_table[-1][0]
        for ele in self.curv_spd_table:
            if ele[1] <= v:
                max_curv = ele[0]
        # get max front wheel steering angle
        max_delta = min(self.max_steering_angle, math.atan(max_curv*self.axis_dist))

        half = int(sample_num/2) if sample_num%2 == 1
        sample_indices = np.linspace(-int(sample_num/2), int(sample_num/2), sample_num, endpoint=True)
        for idx in sample_indices:
            next_delta = delta+idx*dt*max_steering_spd
