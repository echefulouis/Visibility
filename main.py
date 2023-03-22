from ObstaclePolygon import ObstaclePolygon
from VisibilityRoadmap import VisibilityRoadMap
import matplotlib.pyplot as plt
from FieldView import FieldView
show_animation = True

def main():
    print(__file__ + " start!!")

    # start and goal position
    sx, sy = 10.0, 10.0  # [m]
    gx, gy = 50.0, 50.0  # [m]

    expand_distance = 5.0  # [m]

    obstacles = [
        ObstaclePolygon(
            [20.0, 30.0, 15.0],
            [20.0, 20.0, 30.0],
        ),
        ObstaclePolygon(
            [40.0, 45.0, 50.0, 40.0],
            [50.0, 40.0, 20.0, 40.0],
        ),
        ObstaclePolygon(
            [20.0, 30.0, 30.0, 20.0],
            [40.0, 45.0, 60.0, 50.0],
        )
    ]

    if show_animation:  # pragma: no cover
        plt.plot(sx, sy, "or")
        plt.plot(gx, gy, "ob")
        for ob in obstacles:
            ob.plot()
        plt.axis("equal")
        plt.pause(1.0)

    rx, ry = VisibilityRoadMap(expand_distance, do_plot=show_animation) \
        .planning(sx, sy, gx, gy, obstacles)

    fov = FieldView(fov=30)
    for i, (x, y) in enumerate(zip(rx, ry)):
        if i == 0: continue
        robot = [rx[i-1], ry[i-1]]
        goal = [x, y]
        FOV = fov.clippedFOV(robot, goal, obstacles)
        plt.fill(*FOV.exterior.xy, alpha=0.4)


    if show_animation:  # pragma: no cover
        plt.plot(rx, ry, "-r")
        plt.pause(0.1)
        plt.show()


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    main()
