import pickle
from racetracks import *
from graph_node import Node
import matplotlib.pyplot as plt
import math

seed = np.random.seed(1234)
# a dictionary contains all [key, state] pairs
# each state has a g value and connects to other states through kinametic equations
graph = {}


# save generated graph to a .dat file by pickle
def build_up_graph(grid, save_path):
    max_vel = 5

    # velocity dimension
    vel_list = []
    for i_vel in range(-max_vel+1, max_vel):
        for j_vel in range(-max_vel+1, max_vel):
            vel_list.append([i_vel, j_vel])

    # position dimension
    x_idx, y_idx = np.where(grid == FREE)
    coord = np.stack([x_idx, y_idx], axis=1)
    for p_idx in range(coord.shape[0]):
        pnt = coord[p_idx]
        for vel in vel_list:
            state = Node(pnt[0], pnt[1], vel[0], vel[1])
            state.connect_to_graph(grid)
            graph[state.key] = state

    for pnt in START_LINE:
        state = Node(pnt[0], pnt[1], 0, 0)
        state.connect_to_graph(grid)
        graph[state.key] = state

    for pnt in FINISH_LINE:
        state = Node(pnt[0], pnt[1], 0, 0)
        state.is_goal = True
        graph[state.key] = state

    output = open(save_path, 'wb')
    pickle.dump(graph, output)



def check_graph(grid):
    plt.figure(figsize=(4.5, 16))
    plt.pcolor(grid, edgecolors='k', linewidths=1)
    for key in graph.keys():
        for child_idx, child_key in enumerate(graph[key].next_prob_1): # or next_prob_1
            ux, uy = ACTION_SPACE[child_idx]
            vx, vy = graph[key].vx + ux,  graph[key].vy + uy
            child = graph[child_key]
            # check a specific connection
            # plt.title(str(vy) + '_' + str(vx))
            # plt.show()
            if [child.px, child.py] in START_LINE:
                # print('found')
                continue
            plt.arrow(graph[key].py + 0.5, graph[key].px + 0.5,
                      child.py - graph[key].py, child.px - graph[key].px,
                      color='r', head_width=0.3, head_length=0.1)
            # print(key, child_idx)
        # end for
    # end for
    plt.show()


def track_the_best_plan(idx = 0):
    start_node = Node(START_LINE[idx][0], START_LINE[idx][1], 0, 0)
    start_key = start_node.key
    state = graph[start_key]
    trajectory = [state]
    # for i in range(grid.shape[0]+grid.shape[1]) a safer condition
    while not state.is_goal:
        value_uk = []
        for child_idx in range(len(ACTION_SPACE)):
            child_key_9 = state.next_prob_9[child_idx]
            child_9 = graph[child_key_9]
            value_uk.append(child_9.g_value)
        child_key = state.next_prob_9[np.argmin(value_uk)]
        state = graph[child_key]
        trajectory.append(state)
        print(state.px, state.py)
    return trajectory



def visualize_the_best_plan(plan, grid_para):
    assert isinstance(plan, list)
    plt.figure(figsize=(4.5, 16))
    plt.pcolor(grid_para, edgecolors='k', linewidths=1)
    plan_len = len(plan)
    plan.append(plan[-1])
    for i in range(plan_len):
        plt.arrow(plan[i].py + 0.5, plan[i].px + 0.5,
                  plan[i+1].py - plan[i].py, plan[i+1].px - plan[i].px,
                  color='r', head_width=0.3, head_length=0.1)
    plt.show()



def dynamic_programming():
    itr_num = 0
    bellman_error = np.inf
    bellman_error_list = []
    while bellman_error > 0.0001:
        itr_num += 1
        bellman_error = 0.0
        for key in graph.keys():
            state = graph[key]
            if state.is_goal:
                state.g_value = 0
            else:
                value_uk = []
                for child_idx in range(len(ACTION_SPACE)):
                    child_key_9 = state.next_prob_9[child_idx]
                    child_9 = graph[child_key_9]
                    child_key_1 = state.next_prob_1[child_idx]
                    child_1 = graph[child_key_1]
                    # at first each action's cost are the same(= 1)
                    expected_cost_uk = 0.9 * (1 + child_9.g_value) + 0.1 * (1 + child_1.g_value)
                    value_uk.append(expected_cost_uk)
                current_value = min(value_uk)
                bellman_error += np.linalg.norm(state.g_value - current_value)
                state.g_value = min(value_uk)
            # end if
        # end for
        bellman_error_list.append(bellman_error)
        print("{}th iteration: {}".format(itr_num, bellman_error))
    # end while

    plt.figure()
    x_axis = range(len(bellman_error_list))
    plt.plot(x_axis, bellman_error_list)
    plt.show()

def get_heuristic(px0, py0, px1, py1, mode = 0):
    dx = abs(px0 - px1)
    dy = abs(py0 - py1)
    if mode == 0:
        # L1 norm (manhatton norm)
        dist = dx + dy
    elif mode == 1:
        # L2 nrom (euclidean norm)
        dist = math.sqrt(dx*dx + dy*dy)
    elif mode == 2:
        # Linf norm (max element norm)
        dist = max(dx, dy)
    elif mode == 3:
        # diagonal norm
        d_max = max(dx, dy)
        d_min = min(dx, dy)
        dist = d_max - d_min + math.sqrt(2*d_min*d_min)
    elif mode == 4:
        # Linf norm / max spd
        dist = max(dx, dy) / 4
    elif mode == 5:
        # L2 norm / max spd
        dist = math.sqrt(dx*dx + dy*dy) / 4
    elif mode == 6:
        # diagonal norm / max spd
        d_max = max(dx, dy)
        d_min = min(dx, dy)
        dist = d_max - d_min + math.sqrt(2*d_min*d_min)
        dist = dist / 4
    elif mode == 7:
        # L1 norm (manhatton norm) /  max spd
        dist = dx + dy
        dist = dist / 4
    return dist

def start_check(input_key):
    for pnt in START_LINE:
        state = Node(pnt[0], pnt[1], 0, 0)
        if input_key == state.key:
            return 1
    return 0

def real_time_dynamic_programming(mode):
    iter_num = 0
    bellman_error_list = []
    random_prob = 0.0
    # step 1: Initialize G values of all states to admissible values
    for key in graph.keys():
        state = graph[key]
        if state.is_goal:
            state.g_value = 0.0
        else:
            g_value_list = []
            for goal in FINISH_LINE:
                g_value_temp = get_heuristic(state.px, state.py, goal[0], goal[1], mode)
                g_value_list.append(g_value_temp)
            state.g_value = min(g_value_list)
        # print('[INIT] g_value: ', state.g_value)

    while 1:
        # step 2: Follow greedy policy picking outcomes at random until goal is reached 
        iter_num += 1
        bellman_error_start_list = []
        bellman_error_start = np.inf
        for start_pnt in START_LINE:
            bellman_error = 0.0
            start_state = Node(start_pnt[0], start_pnt[1], 0, 0)
            current_state = graph[start_state.key]
            while not current_state.is_goal:
                # find next state by bellman expectation equation
                min_value = np.inf
                min_idx = 0
                first_flag = 0
                value_list = []
                idx_list = []
                for child_idx in range(len(ACTION_SPACE)):
                    child_key_9 = current_state.next_prob_9[child_idx]
                    child_9 = graph[child_key_9]
                    # skip start node(might be occupied or out of bound)
                    if start_check(child_key_9) == 1 and first_flag == 1:
                        continue
                    child_key_1 = current_state.next_prob_1[child_idx]
                    child_1 = graph[child_key_1]
                    expected_cost_uk = 0.9 * (1 + child_9.g_value) + 0.1 * (1 + child_1.g_value)
                    value_list.append(expected_cost_uk)
                    idx_list.append(child_idx)
                    if first_flag == 0:
                        first_flag = 1
                        min_value = expected_cost_uk
                        min_idx = child_idx
                    elif min_value > expected_cost_uk:
                        min_value = expected_cost_uk
                        min_idx = child_idx  
                current_value = 0.0                  
                if np.random.rand() < random_prob:
                    rand_idx_inner = np.random.randint(low=0, high=len(idx_list)-1, size=1)[0]
                    rand_idx = idx_list[rand_idx_inner]
                    next_key = current_state.next_prob_9[rand_idx]
                    current_value = value_list[rand_idx_inner]
                    # print('[ITER] next_key(rand): ', next_key)
                else:
                    # go to next state
                    next_key = current_state.next_prob_9[min_idx]
                    current_value = min_value
                    # print('[ITER] next_key(min): ', next_key)
                # update bellman error
                bellman_error += np.linalg.norm(current_state.g_value - current_value)
                # step 3: Backup all states visited on the way
                current_state.g_value = current_value
                current_state = graph[next_key]
                # print('[ITER] bellman error: ', bellman_error)
            # collect bellman error for each start point
            bellman_error_start_list.append(bellman_error)

        bellman_error_start = max(bellman_error_start_list)
        bellman_error_list.append(bellman_error_start)
        print(iter_num, 'th iteration, error: ', bellman_error_start)
        if bellman_error_start <= 0.0001:
            print('Bellman error <= 0.0001, finished iteration !!!')
            break    
 
    plt.figure()
    x_axis = range(len(bellman_error_list))
    plt.plot(x_axis, bellman_error_list)
    plt.show()

if __name__ == '__main__':
    path = './solution/graph_dp.dat'
    track_map = race_track
    # build_up_graph(track_map, path)
    graph = pickle.load(open(path, 'rb'))

    # solve
    # dynamic_programming()
    real_time_dynamic_programming(4)
    for i in range(len(START_LINE)):
        plan = track_the_best_plan(i)
        visualize_the_best_plan(plan, track_map)
    # check_graph(track_map)
