import scipy.io
import numpy as np
import yaml

from dependentcode import Environment, nextstate, NetworkxNeighbours
from FieldView import FieldView
from main import read_csv_from_dir, interpolate_points_between
from shapely.geometry import Point
from collections import defaultdict
import matplotlib.pyplot as plt
from functools import partial
from tqdm import tqdm


def detection(start_point, goal_point, history, k, fov, targetStates, VMAX):
    trajectory = interpolate_points_between(start_point, goal_point, VMAX)
    for l, goalcord in enumerate(trajectory):

        if (l == 0) or (l == len(trajectory) - 1):
            continue

        robotcord = trajectory[l - 1]
        trajData = {"start": robotcord, "goal": goalcord}
        history[k].append(trajData)

        try:
            FOV = fov.clippedFOV(robotcord, goalcord)

            for t in targetStates:
                tcord = nx_points[t]
                tpoint = Point(tcord[0], tcord[1])
                if FOV.contains(tpoint):

                    # print(k, "catches the target ", t)
                    return t
        except:
            pass

    return None

def patrolling_simulation(initial_robot_states, divs, DMCs, number_iterations, targetStates, nx_points, targets):
    obstacles = list(read_csv_from_dir('obstacles'))
    fov = FieldView(obstacles=obstacles, fov=120)


    catchStates = set()
    robot_states = initial_robot_states.reshape(-1, 1)
    robot = np.squeeze(robot_states.copy())
    VMAX = 5 # m/s
    history = defaultdict(list)
    # print('robot',robot, robot_states)
    for j in range(number_iterations):
        nextstatevalue = []
        for k, state in enumerate(robot):

            nstate = nextstate(divs[k].index(state), DMCs[k])
            goalstate = divs[k][nstate]
            start_point = nx_points[state]
            goal_point = nx_points[goalstate]
            nextstatevalue.append(goalstate)

            if len(catchStates) > 0:
                stationaryRobotIndexes = [kIndex for (_, kIndex) in catchStates]
                if k in stationaryRobotIndexes:
                    continue

            t = detection(start_point, goal_point, history, k, fov, targetStates, VMAX)
            if t is not None:
                catchStates.add((t, k))

        if len(catchStates) == len(targetStates):
            break

    return history, catchStates

def display_animation(history):
    obstacles = list(read_csv_from_dir('obstacles'))
    fov = FieldView(obstacles=obstacles, fov=120)
    colors = ['r', 'g', 'b', 'cyan', 'brown']
    def getPath(i):
        path = []

        for data in history[i]:
            # print(data)
            FOV = fov.clippedFOV(data['start'], data['start'])
            plt.fill(*FOV.exterior.xy, alpha=0.4, color=colors[i])
            path.append(data['start'])
        return np.array(path)

    paths = list(map(getPath, history))
    print(len(paths))
    for ob in obstacles:
        ob.plot()
    for i, path in enumerate(paths):
        plt.plot(path[:, 0], path[:, 1])



def compute_path_length(path):
    total_length = 0
    for i in range(len(path) - 1):
        total_length += np.linalg.norm(path[i - 1] - path[i])
    return total_length

def eval(num_trials, targets, eval_func):
    def getPath(pathHistory):
        path = []
        for data in pathHistory:
            path.append(data['start'])
        return np.array(path)

    count = 0
    totalPathLen = 0
    for _ in tqdm(range(num_trials)):
        # print(f'epoch {i + 1} of {num_trials}')
        history, catchStates = eval_func()
        for (_, k) in catchStates:
            path = getPath(history[k])
            # print(path)
            totalPathLen += compute_path_length(path)
            # print(totalPathLen)

        if len(catchStates) == len(targets):
            count += 1

    return {'avgPathLength': totalPathLen/count, 'detectionRate': count/num_trials}


def main(param):
    env_file_name = param['env_file_name']
    new_environment = Environment(env_file_name)
    new_environment.parse_obstacles()
    new_environment.find_vertical_lines()
    nx_points = new_environment.get_env_graph_points()

    divs = [param['MTSP'][f'robot{i + 1}'] for i in range(param['num_robots'])]
    DMC_FILEs = [param['DMC'][f'robot{i + 1}'] for i in range(param['num_robots'])]
    DMCs = list(map(lambda fname: scipy.io.loadmat(fname)['Ser'], DMC_FILEs))
    initial_robot_states = np.array(param['initial_robot_states'])
    target_states = param['target_states']
    number_trails = param['number_trails']
    number_iterations = param['number_iterations']

    targets = []
    for t in target_states:
        targets.append(nx_points[t].tolist())
    eval_func = partial(patrolling_simulation, initial_robot_states, divs, DMCs, number_iterations, target_states,
                        nx_points, np.array(targets))

    metric = eval(number_trails, targets, eval_func)
    print(metric)

if __name__ == "__main__":
    with open('param.yaml', 'r') as f:
        param = yaml.safe_load(f)
    # main(param)

    input_file_name = 'resources/input_file'
    new_environment = Environment(input_file_name)
    new_environment.parse_obstacles()
    new_environment.find_vertical_lines()
    Linelist = new_environment.LineList
    nx_points = new_environment.get_env_graph_points()
    nx_obstacles = new_environment.get_env_new_obstacles()

    div1 = [0, 2, 14, 10, 12, 26, 29, 25, 27, 68, 71, 70, 72, 74, 67, 65, 66, 23, 18, 17, 15, 16]
    div2 = [34, 32, 30, 31, 33, 69, 22, 20, 21, 24, 19, 7, 8, 5, 9, 6, 84, 81, 80, 82, 4]
    div3 = [39, 35, 36, 38, 37, 51, 53, 50, 54, 52, 59, 56, 55, 57, 64, 79, 86, 85, 87, 89, 47, 49, 45, 48, 46]
    div4 = [1, 83, 3, 11, 13, 28, 73, 41, 44, 40, 43, 42, 78, 88, 77, 75, 76, 63, 62, 60, 61, 58]
    mat = scipy.io.loadmat('resources/MCR1.mat')
    DMC1 = mat['Ser']
    mat = scipy.io.loadmat('resources/MCR2.mat')
    DMC2 = mat['Ser']
    mat = scipy.io.loadmat('resources/MCR3.mat')
    DMC3 = mat['Ser']
    mat = scipy.io.loadmat('resources/MCR4.mat')
    DMC4 = mat['Ser']

    DMCs = [DMC1, DMC2, DMC3, DMC4]
    divs = [div1, div2, div3, div4]
    initial_robot_states = np.array([14, 5, 89, 40])
    target_states = [15, 25]
    number_iterations = 100
    number_trails = 100
    # Create the Networx_x Graph
    number_of_neigbours = 9
    env_nx_graph = NetworkxNeighbours(nx_points, nx_obstacles, number_of_neigbours)
    G = env_nx_graph.get_neighbours()
    #eval=MCEvaluation( divs, DMCs, nx_points, G, target_states, initial_robot_states, number_iterations, number_trails)
    #print(eval)


    targets = []
    for t in target_states:
        targets.append(nx_points[t].tolist())
    eval_func = partial(patrolling_simulation, initial_robot_states, divs, DMCs, number_iterations, target_states, nx_points, np.array(targets))

    metric = eval(number_trails, targets, eval_func)
    print(metric)


    plt.figure()
    vizTargets = np.array(targets)
    history, catchStates = patrolling_simulation(initial_robot_states, divs, DMCs, number_iterations, target_states,
                                                 nx_points, np.array(targets))
    display_animation(history)
    plt.scatter(vizTargets[:, 0], vizTargets[:, 1], 50, c='k')
    plt.show()

