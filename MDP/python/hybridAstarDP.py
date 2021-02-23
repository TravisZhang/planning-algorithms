import pickle
from racetracks import *
from vehicle_graph_node import VehicleNode
import matplotlib.pyplot as plt
import math

seed = np.random.seed(1234)
# a graph indicates its node's dist to nearest occupied node
graph = np.zeros(race_track.shape, dtype = float) 


# save generated graph to a .dat file by pickle
# build up gradient graph by calculate dist to nearest occupied Node
def build_up_graph(grid, save_path):
    # first set free nodes value to inf
    x_idx_free, y_idx_free = np.where(grid == FREE)
    coord_free = np.stack([x_idx_free, y_idx_free], axis=1)
    for type_name in [START, FINISH]:
        x_idx, y_idx = np.where(grid == type_name)
        coord = np.stack([x_idx, y_idx], axis=1)
        coord_free = np.concatenate((coord_free,coord),axis=0)
    for p_idx in range(coord_free.shape[0]):
        pnt = coord_free[p_idx]
        graph[pnt[0], pnt[1]] = math.inf

    print(graph)

    # update dist to nearest occupied node for free nodes
    # may use kd-tree to store points to increase find speed
    x_idx_occ, y_idx_occ = np.where(grid == OCCUPIED)
    coord_occ = np.stack([x_idx_occ, y_idx_occ], axis=1)
    for p_idx0 in range(coord_occ.shape[0]):
        pnt0 = coord_occ[p_idx0]
        graph[pnt0[0], pnt0[1]] = -1.0
        for p_idx1 in range(coord_free.shape[0]): 
            pnt1 = coord_free[p_idx1]
            dist_x = pnt1[0] - pnt0[0]
            dist_y = pnt1[1] - pnt0[1]
            dist = math.sqrt(dist_x*dist_x+dist_y*dist_y)
            if graph[pnt1[0], pnt1[1]] > dist:
                graph[pnt1[0], pnt1[1]] = dist

    output = open(save_path, 'wb')
    pickle.dump(graph, output)



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

if __name__ == '__main__':
    path = './solution/graph_dp.dat'
    track_map = race_track
    build_up_graph(track_map, path)
    print(graph)
    # graph = pickle.load(open(path, 'rb'))

