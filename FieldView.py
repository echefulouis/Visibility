import math
import numpy as np
from shapely.geometry import Polygon, Point, LineString

class FieldView:
    def __init__(self, obstacles, fov):
        self.obstacles = obstacles
        self.fov = fov
    @staticmethod
    def sector(center, start_angle, end_angle, radius, steps=200):
        '''
        # https://stackoverflow.com/questions/54284984/sectors-representing-and-intersections-in-shapely
        :param center: camera location
        :param start_angle: starting angle from x axis
        :param end_angle: ending angle from y axis
        :param radius: the range of FOV
        :param steps: how many points on FOV
        :return: Polygon (shapely.geometry.Polygon)
        '''
        def polar_point(origin_point, angle, distance):
            return [origin_point.x + math.sin(math.radians(angle)) * distance,
                    origin_point.y + math.cos(math.radians(angle)) * distance]

        if start_angle > end_angle:
            start_angle = start_angle - 360
        else:
            pass
        step_angle_width = (end_angle - start_angle) / steps
        sector_width = (end_angle - start_angle)
        segment_vertices = []

        segment_vertices.append(polar_point(center, 0, 0))
        segment_vertices.append(polar_point(center, start_angle, radius))

        for z in range(1, steps):
            segment_vertices.append((polar_point(center, start_angle + z * step_angle_width, radius)))
        segment_vertices.append(polar_point(center, start_angle + sector_width, radius))
        segment_vertices.append(polar_point(center, 0, 0))
        return Polygon(segment_vertices)


    def getFOV(self, x, y, theta, radius):
        '''
        parameters
        ==========
        @param x: x coordinate
        @param y: y coordinate
        @param theta: heading angle
        @param radius: sensing range
        @return Visibility Polygon
        '''

        center = Point(x, y)
        PHI_OFFSET = self.fov / 2
        THETA_OFFSET = 90
        theta = THETA_OFFSET - theta
        phi1 = theta - PHI_OFFSET
        phi2 = theta + PHI_OFFSET
        sect = FieldView.sector(center, phi1, phi2, radius)
        return sect


    def sectorView(self, robot, goal):
        '''
        :param robot:
        :param goal:
        :return:
        '''

        robot = np.array(robot)
        goal = np.array(goal)
        diff = goal - robot
        theta = np.arctan2(diff[1], diff[0])
        theta = np.rad2deg(theta)
        # radius = np.linalg.norm(diff)
        radius = 50
        x, y = robot

        return self.getFOV(x, y, theta, radius)

    def clippedFOV(self, robot, goal):
        currentFOV = self.sectorView(robot, goal)
        # convert obstacles to polygon
        def toPolygons(obstacles):
            result = []
            for obs in obstacles:
                points = []
                for x, y in zip(obs.x_list, obs.y_list):
                    p = Point(x, y)
                    points.append(p)
                result.append(Polygon(points))
            return result

        obsPolys = toPolygons(self.obstacles)

        def ray_trace(fov, poly):
            coords = fov.exterior.coords
            N = len(coords)
            start = Point(coords[0])
            result = [start]
            for i in range(1, N - 1):
                line = LineString([start, coords[i]])
                if line.intersects(poly):
                    dist = np.inf
                    chosenPoint = None
                    intersections = line.intersection(poly)
                    # print(intersections)
                    for p in intersections.coords:
                        # compute minimum distance
                        point = Point(p[0], p[1])
                        newdist = start.distance(point)
                        if newdist < dist:
                            chosenPoint = point
                            dist = newdist
                    result.append(chosenPoint)
                else:
                    result.append(Point(coords[i][0], coords[i][1]))
            result.append(start)
            return Polygon(result)

        for poly in obsPolys:
            if currentFOV.intersects(poly):
                currentFOV = ray_trace(currentFOV, poly)

        return currentFOV