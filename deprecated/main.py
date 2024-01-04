from ObstaclePolygon import ObstaclePolygon
import matplotlib.pyplot as plt
from FieldView import FieldView
import pathlib
import numpy as np
import yaml
from collections import defaultdict
show_animation = True
def read_csv(filename):
    return np.loadtxt(filename, delimiter=",")

def read_csv_from_dir(dirname):
    path = pathlib.Path(dirname)
    for filename in path.glob("*.csv"):
        xy = read_csv(filename)
        yield ObstaclePolygon(xy[:, 0].tolist(),xy[:, 1].tolist())


def display_animation(fov, obstacles, vizData, targets):
    simTime = max(map(len, vizData.values()))
    print(simTime)
    for t in range(simTime):
        plt.cla()
        plt.plot(targets[:,0], targets[:,1], "ro")
        # display obstacles
        for ob in obstacles:
            ob.plot()
        for robotName, value in vizData.items():
            k = min(t, len(value) - 1)
            robot = value[k]['robot']
            goal = value[k]['goal']
            color = value[k]['color']
            try:
                FOV = fov.clippedFOV(robot, goal)
                plt.fill(*FOV.exterior.xy, alpha=0.4, color=color)
            except:
                pass



            plt.scatter(robot[0], robot[1], c=color)
        plt.axis([-10, 300, -10, 200])
        plt.pause(0.1)
    plt.show()
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

def main():
    print(__file__ + " start!!")
    result = "result.yaml"
    xy_coord = read_csv('pvis.csv')
    with open(result) as file:
        data = yaml.safe_load(file)

    obstacles = list(read_csv_from_dir('obstacles'))
    targets = np.array(data["Targets"])
    fov = FieldView(obstacles=obstacles, fov=120)


    colors = ['r', 'g', 'b', 'cyan', 'brown']
    vizData = defaultdict(list)
    VMAX = 5.0
    for i in range(data['NUM_ROBOTS']):
        name = f'robot{i+1}'
        indexes = data[name]
        for j in range(1, len(indexes)):
            robot = xy_coord[indexes[j-1], :2]
            goal = xy_coord[indexes[j], :2]
            color = colors[i]
            vizData[name].append({'robot': robot, 'goal': goal, 'color': color})
            # path = interpolate_points_between(robot, goal, VMAX)
            # for k, g in enumerate(path):
            #     if k == 0 or k == len(path) - 1:
            #         continue
            #     p = path[k-1]
            #     vizData[name].append({'robot': p, 'goal':g, 'color':color})

    display_animation(fov, obstacles, vizData, targets)



# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
