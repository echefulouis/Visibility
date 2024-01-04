import numpy as np
from shapely.geometry import Point
from threading import Thread
import yaml
import scipy.io
from dependentcode import Environment, nextstate, NetworkxNeighbours
from FieldView import FieldView
from ObstaclePolygon import ObstaclePolygon
import pathlib
from collections import Counter


def read_csv_from_dir(dirname):
    def read_csv(filename):
        return np.loadtxt(filename, delimiter=",")
    path = pathlib.Path(dirname)
    for filename in path.glob("*.csv"):
        xy = read_csv(filename)
        yield ObstaclePolygon(xy[:, 0].tolist(),xy[:, 1].tolist())


def interpolate_points_between(start_point, goal_point, fixed_distance):
    # Convert points to NumPy arrays
    start_point = np.array(start_point)
    goal_point = np.array(goal_point)

    # Calculate total distance
    total_distance = np.linalg.norm(goal_point - start_point)

    # Calculate the number of segments
    num_segments = int(total_distance // fixed_distance)

    # Calculate the fraction for each segment
    fractions = np.linspace(0, 1, num_segments + 1, endpoint=True)

    # Interpolate coordinates for each fraction
    interpolated_points = [start_point + frac * (goal_point - start_point) for frac in fractions]

    return interpolated_points

class PatrolRobot(Thread):
    VMAX = 5
    def __init__(self, initState, div, dmc, fov, transition, coordMap, targetStates, maxSteps, CATCHES, verbose=True):
        self.state = initState
        self.div = div
        self.dmc = dmc
        self.fov = fov
        self.transition = transition
        self.coordMap = coordMap
        self.targetStates = targetStates
        self.maxSteps = maxSteps
        self.success = False
        self.trajectory = []
        self.CATCHES = CATCHES
        self.verbose = verbose
        Thread.__init__(self)


    def toPoint(self, state):
        return self.coordMap[state]

    @property
    def pathLen(self):
        return self.compute_path_length(self.trajectory)

    @staticmethod
    def compute_path_length(path):
        total_length = 0
        for i in range(len(path) - 1):
            total_length += np.linalg.norm(path[i - 1] - path[i])
        return total_length

    def __repr__(self):
        return f"[{self.name}] pathLen = {self.pathLen:.3f}"

    def detectAlongTraj(self, trajectory):
        for l, goalcord in enumerate(trajectory):

            if (l == 0) or (l == len(trajectory) - 1):
                continue
            robotcord = trajectory[l - 1]
            self.trajectory.append(robotcord)
            try:
                FOV = self.fov.clippedFOV(robotcord, goalcord)
                for tstate in self.targetStates:
                    tcord = self.toPoint(tstate)
                    tpoint = Point(tcord[0], tcord[1])
                    if FOV.contains(tpoint):
                        if self.verbose:
                            print(self.name, f"catches {tstate} target")
                        self.CATCHES.add(self.name)
                        return True
            except:
                pass

        return False

    def run(self):
        step = 0
        while step < self.maxSteps:
            nstate = self.transition(self.div.index(self.state), self.dmc)
            goalState = self.div[nstate]
            start_point = self.toPoint(self.state)
            goal_point = self.toPoint(goalState)
            trajectory = interpolate_points_between(start_point, goal_point, self.VMAX)
            if self.detectAlongTraj(trajectory):
                self.success = True
                break
            if len(self.CATCHES) == len(self.targetStates):
                break
            self.state = goalState
            step += 1

def search(param, verbose=True):
    CATCHES = set()
    env_file_name = param['env_file_name']
    new_environment = Environment(env_file_name)
    new_environment.parse_obstacles()
    new_environment.find_vertical_lines()
    coordMap = new_environment.get_env_graph_points()

    divs = [param['MTSP'][f'robot{i + 1}'] for i in range(param['num_robots'])]
    DMC_FILEs = [param['DMC'][f'robot{i + 1}'] for i in range(param['num_robots'])]
    DMCs = list(map(lambda fname: scipy.io.loadmat(fname)['Ser'], DMC_FILEs))
    initial_robot_states = param['initial_robot_states']
    target_states = param['target_states']
    number_iterations = param['number_iterations']

    obstacles = list(read_csv_from_dir('obstacles'))
    fov = FieldView(obstacles=obstacles, fov=120)

    patrols = {}
    for initState, div, dmc in zip(initial_robot_states, divs, DMCs):
        pAgent = PatrolRobot(initState, div, dmc, fov, nextstate, coordMap, target_states, number_iterations, CATCHES, verbose)
        pAgent.start()
        patrols[pAgent.name] = pAgent

    for agent in patrols.values():
        agent.join()
        # print(agent)
    return patrols, CATCHES


if __name__ == '__main__':
    with open('param.yaml') as file:
        param = yaml.safe_load(file)
    patrols, catchers = search(param)
    numSuccess = Counter([p.success for p in patrols.values()])[True]
    print(patrols, catchers)
    print(numSuccess)

    import matplotlib.pyplot as plt
    plt.cla()
    colors = ['r', 'g', 'b', 'cyan', 'brown']
    obstacles = list(read_csv_from_dir('obstacles'))

    N = max(map(lambda a: len(a.trajectory), patrols.values()))
    for t in range(N):
        if t == 0:
            continue
        # plt.cla()
        for ob in obstacles:
            ob.plot()

        for i, agent in enumerate(patrols.values()):
            if t < len(agent.trajectory):
                fov = agent.fov
                try:
                    FOV = fov.clippedFOV(agent.trajectory[t-1], agent.trajectory[t])
                    plt.fill(*FOV.exterior.xy, alpha=0.4, color=colors[i])
                except Exception as error:
                    print(error)
        plt.axis([-10, 300, -10, 200])
        plt.pause(0.1)
    plt.show()


