import math
from racetracks import *
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
    max_jerk = 1.0
    max_spd = 3.0
    steer_ratio = 17.0
    max_steering_angle = 480.0/steer_ratio*KtoRad
    max_steering_spd = 30.0/steer_ratio*KtoRad
    max_steering_spd_low_spd = 720.0/steer_ratio*KtoRad

    def __init__(self, v, a, delta, dt, lon_sample_num, lat_sample_num):
        max_curv = self.curv_spd_table[-1][0]
        for ele in self.curv_spd_table:
            if ele[1] <= v:
                max_curv = ele[0]
        # get max front wheel steering angle
        max_delta = min(self.max_steering_angle, math.atan(max_curv*self.axis_dist))

        # sample lat actions (steering angles)
        self.delta_set = self.get_samples(delta, max_delta, self.max_steering_spd, 0.005, lat_sample_num, dt)

        # sample lon actions (acc)
        self.acc_set = self.get_samples(a, self.max_acc, self.max_jerk, 0.1, lon_sample_num, dt)

    def get_samples(self, input, max_input, max_input_spd, min_delta_input_speed, sample_num, dt):
        sample_set = list()
        half_num = max(int(sample_num/2),1)
        sample_indices = np.linspace(-half_num, half_num, sample_num, endpoint=True)
        delta_input_spd = max(max_input_spd / max((len(sample_indices)-1),1), min_delta_input_speed)
        for idx in sample_indices:
            next_input = input+idx*dt*delta_input_spd
            if next_input >= max_input and next_input <= max_input:
                sample_set.append(next_input)
        if len(sample_set) == 0:
            sample_set.append(input)
        return sample_set

    def delta_set(self):
        return self.delta_set

    def acc_set(self):
        return self.acc_set

    def max_acc(self):
        return self.max_acc

    def max_jerk(self):
        return self.max_jerk

    def max_spd(self):
        return self.max_spd

    def axis_dist(self):
        return self.axis_dist


# state defined bicycle model
# px,py: vehicle position
# yaw: vehicle heading
# delta: vehicle front wheel steering angle
# v: vehicle longitudinal speed
# a: vehicle longitudinal acc
class VehicleNode:
    # state id(static)
    next_id = 1

    def __init__(self, px, py, yaw, delta, v, a, dt, lon_act_num, lat_act_num):
        # state
        self.px = px
        self.py = py
        self.yaw = yaw
        self.delta = delta
        self.v = v
        self.a = a
        self.dt = dt
        self.id = VehicleNode.next_id
        VehicleNode.next_id += 1
        # value
        self.g_value = 0.0
        # successor
        # contains keys for next states(with different actions(acc and steering) as input, stored with the same sequence as actions)
        self.next_prob_9 = []
        self.next_prob_1 = []
        # key
        self.key = self.get_key()
        self.is_goal = False
        # vehicle actions
        vas = VehicleActionSpace(self.v, self.a, self.delta, self.dt, lon_act_num, lat_act_num)
        self.delta_set = vas.delta_set()
        self.acc_set = vas.acc_set()
        self.max_spd = vas.max_spd()
        self.axis_dist = vas.axis_dist()

    def get_id(self):
        return self.id

    def get_key(self):
        return self.generate_key(self.px, self.py, self.v, self.yaw, self.delta, self.v, self.a)


    def connect_to_graph(self, grid):
        for delta in self.delta_set:
            for acc in self.acc_set:
                self.next_prob_9.append(self.control(delta, acc, grid, success=True))
                self.next_prob_1.append(self.control(u[0], u[1], grid, success=False))


    def velocity_constraints(self, v):
        # max vel currently at 4 m/s
        return np.sign(v) * min(abs(v), self.max_spd)

    def safety_constraints(self, px2, py2, grid):
        assert  0 <= self.px < grid.shape[0]
        assert  0 <= self.py < grid.shape[1]

        x_dist = np.abs(px2 - self.px)
        y_dist = np.abs(py2 - self.py)
        step = max(x_dist, y_dist)
        x_way_points = np.linspace(self.px, px2, step + 1, endpoint=True)
        y_way_points = np.linspace(self.py, py2, step + 1, endpoint=True)
        way_points = np.stack([np.ceil(x_way_points), np.ceil(y_way_points)], axis=1).astype(np.int)

        for idx in range(way_points.shape[0]):
            point = way_points[idx]
            if (0 <= point[0] < grid.shape[0]) and (0 <= point[1] < grid.shape[1]):
                if grid[point[0], point[1]] == FINISH:
                    return FINISH, point
                elif grid[point[0], point[1]] == OCCUPIED:
                    return OCCUPIED, point
                # else:
                #  free and start: continue
            else:
                return OUTBOUND, point

        if grid[way_points[-1][0], way_points[-1][1]] == START:
            return START, way_points[-1]
        else:
            return FREE, way_points[-1]
    # end definition


    def control(self, delta, a, grid, success):
        # success with probability of 0.9
        # if not success, stay put(zero acc input)
        if not success:
            delta = self.delta
            a = self.a

        # bicycle model
        v = self.v + a*self.dt
        v = self.velocity_constraints(v)
        vx = v*math.cos(self.yaw+delta)
        vy = v*math.sin(self.yaw+delta)
        px = self.px + vx*self.dt
        py = self.py + vy*self.dt
        yaw = self.yaw + v/self.axis_dist*math.sin(delta)*self.dt

        # check collision
        status, point = self.safety_constraints(px, py, grid)
        if status == FREE:
            assert px == point[0] and py == point[1]
            return self.get_id()
        elif status == START:
            assert grid[point[0], point[1]] == START
            assert px == point[0] and py == point[1]
            return self.generate_key(point[0], point[1], 0, 0)
        elif status == FINISH:
            assert grid[point[0], point[1]] == FINISH
            return self.generate_key(point[0], point[1], 0, 0)
        else: # out of bound or occupied
            assert status == OUTBOUND or status == OCCUPIED
            rand_start = START_LINE[np.random.randint(low=0, high=3, size=1)[0]]
            return self.generate_key(rand_start[0], rand_start[1], 0, 0)