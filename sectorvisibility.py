import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon
import math


def polar_point(origin_point, angle, distance):
    return [origin_point.x + math.sin(math.radians(angle)) * distance,
            origin_point.y + math.cos(math.radians(angle)) * distance]


def sector(center, start_angle, end_angle, radius, steps=200):
    ANGLE_LIMIT = 360
    if start_angle > end_angle:
        start_angle = start_angle - ANGLE_LIMIT

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


def getFOV(x, y, theta, radius):
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
    PHI_OFFSET = 45
    THETA_OFFSET = 90
    theta = THETA_OFFSET - theta
    phi1 = theta - PHI_OFFSET
    phi2 = theta + PHI_OFFSET
    sect = sector(center, phi1, phi2, radius)
    return sect






if __name__ == '__main__':
    radius = 30
    for theta in range(0, 360, 15):
        plt.cla()
        sect = getFOV(0, 0, theta, radius)
        plt.fill(*sect.exterior.xy, color='tab:cyan')
        plt.axis([-radius - 5, radius + 5, -radius - 5, radius + 5])
        plt.pause(1)
