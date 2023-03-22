from pylab import figure, show, rand
from matplotlib.patches import Ellipse, Circle, Rectangle
from matplotlib.patches import Polygon as MPolygon
import geopandas as gpd
import alphashape
from descartes import PolygonPatch
import matplotlib as mpl
import numpy as np
import math
from numpy import *
import itertools
import subprocess
import time
import pylab
import os
from shapely.geometry import Polygon
from collections import defaultdict
import itertools
from shapely.ops import cascaded_union
from itertools import combinations
from shapely.geometry import MultiPolygon
from shapely.geometry import Point

#Imports from decomposition file
import sys
from helpers.graph import *
from helpers.geometry import *
import matplotlib.pyplot as plt
from shapely.geometry import *
from sklearn.neighbors import KNeighborsClassifier ,NearestNeighbors
import networkx as nx
from matplotlib.patches import Wedge
from matplotlib.collections import PatchCollection


from matplotlib.path import Path
from matplotlib.patches import PathPatch

def plot_polygon(ax, poly, **kwargs):
    path = Path.make_compound_path(
        Path(np.asarray(poly.exterior.coords)[:, :2]),
        *[Path(np.asarray(ring.coords)[:, :2]) for ring in poly.interiors])

    patch = PathPatch(path, **kwargs)
    collection = PatchCollection([patch], **kwargs)

    ax.add_collection(collection, autolim=True)
    ax.autoscale_view()
    return collection


def extract_poly_coords(geom):
    if geom.type == 'Polygon':
        exterior_coords = geom.exterior.coords[:]
        interior_coords = []
        for interior in geom.interiors:
            interior_coords.append(list(interior.coords))
    elif geom.type == 'MultiPolygon':
        exterior_coords = []
        interior_coords = []
        for part in geom:
            epc = extract_poly_coords(part)  # Recursive call
            exterior_coords += epc['exterior_coords']
            interior_coords += epc['interior_coords']
    else:
        raise ValueError('Unhandled geometry type: ' + repr(geom.type))
    return {'exterior_coords': exterior_coords,
            'interior_coords': interior_coords}

def writedat(filename, x, y, xprecision=1, yprecision=1):
    with open(filename, 'w') as f:
        for a, b in itertools.zip(x, y):
            print >> f, "%.*g\t%.*g" % (xprecision, a, yprecision, b)


def genVisibMatrix(guardno):
    a = str(os.popen("./main gridnewenvironment2 gridnewguards " + str(guardno)).read())
    l = a.split()
    mat = []
    # print(a)
    for i in range(0, len(l) - 1, 2):
        mat.append([float(l[i]), float(l[i + 1])])
    return mat


def pointOnPolygon(point, poly):
    for i in range(len(poly)):
        a, b = poly[i - 1], poly[i]
        if abs(dist(a, point) + dist(b, point) - dist(a, b)) < EPSILON:
            return true
    return false


def point_in_polyen(x, y, poly):
    # check if point is a vertex
    if (x, y) in poly: return True

    # check if point is on a boundary
    for i in range(len(poly)):
        p1 = None
        p2 = None
        if i == 0:
            p1 = poly[0]
            p2 = poly[1]
        else:
            p1 = poly[i - 1]
            p2 = poly[i]
        if p1[1] == p2[1] and p1[1] == y and x > min(p1[0], p2[0]) and x < max(p1[0], p2[0]):
            return True

    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y

    if inside:
        return inside
    else:
        return inside


def point_in_poly(x, y, poly):

    n = len(poly)
    inside = False

    p1x, p1y = poly[0]
    for i in range(n + 1):
        p2x, p2y = poly[i % n]
        if y > min(p1y, p2y):
            if y <= max(p1y, p2y):
                if x <= max(p1x, p2x):
                    if p1y != p2y:
                        xints = (y - p1y) * (p2x - p1x) / (p2y - p1y) + p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x, p1y = p2x, p2y

    return inside


# https://stackoverflow.com/questions/54284984/sectors-representing-and-intersections-in-shapely
def sector(center, start_angle, end_angle, radius, steps=200):
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


def insidenvironment(x, y):
    if point_in_poly(x, y, env):
        return True
    else:
        return False


def insideobstacle(x, y):
    if point_in_poly(x, y, obs1) or point_in_poly(x, y, obs2) or point_in_poly(x, y, obs3) :
        return True
    else:
        return False


def cellcenter(x, y):
    x_center = (floor(x) + floor(x) + 1) * 0.5
    y_center = (floor(y) + floor(y) + 1) * 0.5
    return x_center, y_center


def cellposition(cell):
    cc = cell % (maxc_x * maxr_y)
    x = cc % maxc_x
    y = floor(cc / maxc_x)
    return (x + (x + 1)) * 0.5, (y + (y + 1)) * 0.5


# mpl.rcParams.update({'font.size': 18})

def midcellnumber(y, x):
    cell = floor(y) * maxc_x + floor(x)

    return int(cell)


def cellnumber(y, x):
    cell = y * maxc_x + x

    return cell


#Function imports from
def parse_input_line(line):
    temp2 = []
    line = [i.strip() for i in line.split(",")]
    vertex = []
    for index,i in enumerate(line):
        if(i[0] == "("):
            i = i[1:]
        if(i[len(i)-1] == ")"):
            i= i[:-1]
        vertex.append(int(i))
        if(index%2 != 0):
            temp2.append(vertex)
            vertex = []
    return temp2

def draw_problem():
    bnd_x = [i.x for i in boundary]
    bnd_x.append(boundary[0].x)
    bnd_y = [i.y for i in boundary]
    bnd_y.append(boundary[0].y)
    poly_x = []
    poly_y = []

    # Draw the boundary
    plt.plot(bnd_x, bnd_y)

    for index, i in enumerate(obstacles):
        poly_x.append([p[0] for p in i])
        poly_y.append([p[1] for p in i])

        plt.fill( poly_x[index], poly_y[index], color="#512DA8")

    plt.plot(source.x, source.y, marker="o")
    plt.plot(dest.x, dest.y, marker="o")
    plt.annotate('Source', xy=(source.x, source.y), xytext=(source.x+5, source.y-6) )
    plt.annotate('Destination', xy=(dest.x, dest.y), xytext=(dest.x-4, dest.y-10) )


def plot_polygon_decomposition(x, y):
    # Create a polygon with the coordinates of the vertices
    polycoords = [(x[0], y[0]), (x[1], y[1]), (x[2], y[2]), (x[3], y[3])]
    # coords = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
    polygon = Polygon(polycoords)

    # distance(shortest)
    line1 = LineString([(x[0], y[0]), (x[1], y[1])])
    line2 = LineString([(x[2], y[2]), (x[3], y[3])])

    distance1 = line1.distance(line2)
    # print(distance1)

    # distance(shortest)
    line3 = LineString([(x[0], y[0]), (x[3], y[3])])
    line4 = LineString([(x[1], y[1]), (x[2], y[2])])

    distance2 = line3.distance(line4)
    # print(distance2)

    # Get the minimum bounding rectangle of the polygon
    min_x, min_y, max_x, max_y = polygon.bounds

    # Find the centroid of the polygon
    # polycentroid = list(polygon.centroid.coords)
    polycentroid = polygon.centroid

    if (distance2 == max_x - min_x):
        # Calculate the height of the polygon along x-axis
        height = max_x - min_x

        # find area and median
        area = polygon.area
        median = area / height
        # print(height)
        # print(median)

        # create x and y coordinate for horizontal line
        y_horizontal = polycentroid.y
        x_horizontal_left = polycentroid.x - height / 4
        x_horizontal_right = polycentroid.x + height / 4

        # calculate x and y coordinate for vertical line
        x_vertical = polycentroid.x
        y_vertical_up = polycentroid.y + median / 4
        y_vertical_down = polycentroid.y - median / 4

        point_x_cord.append(polycentroid.x)

        point_y_cord.append(polycentroid.y)

        point_x_cord.extend([x_horizontal_left, x_horizontal_right, x_vertical, x_vertical])

        point_y_cord.extend([y_horizontal, y_horizontal, y_vertical_up, y_vertical_down])

        # Draw horizontal line through the centriod
        plt.plot([x_horizontal_left, x_horizontal_right], [y_horizontal, y_horizontal], 'x')

        # Draw vertical line through the centriod
        plt.plot([x_vertical, x_vertical], [y_vertical_up, y_vertical_down], 'x')

    elif (distance1 == max_y - min_y):
        # Calculate the height of the polygon
        height = max_y - min_y

        # find area and median
        area = polygon.area
        median = area / height
        # print(height)
        # print(median)

        # create x and y coordinate for vertical line
        x_vertical = polycentroid.x
        y_vertical_up = polycentroid.y - height / 4
        y_vertical_down = polycentroid.y + height / 4

        # create x and y coordinate for horizontal line
        y_horizontal = polycentroid.y
        x_horizontal_left = polycentroid.x - median / 4
        x_horizontal_right = polycentroid.x + median / 4

        point_x_cord.append(polycentroid.x)

        point_y_cord.append(polycentroid.y)

        point_x_cord.extend([x_horizontal_left, x_horizontal_right, x_vertical, x_vertical])

        point_y_cord.extend([y_horizontal, y_horizontal, y_vertical_up, y_vertical_down])

        # Draw horizontal line through the centriod
        plt.plot([x_horizontal_left, x_horizontal_right], [y_horizontal, y_horizontal], 'x')

        # Draw horizontal line through the centriod
        plt.plot([x_vertical, x_vertical], [y_vertical_up, y_vertical_down], 'x')

    # print((polycentroid.x,polycentroid.y))
    plt.plot(x, y)

def main():
    global env, obs1, obs2, obs3, maxc_x, maxr_y, boundary, obstacles,source , dest,point_x_cord ,point_y_cord
    # #LineList = [[0, 0], [8, 0], [8, 8], [0, 8]]
    # #env = [(0, 0), (8, 0), (8, 8), (0, 8)]
    # hole1 = [[2, 5], [5, 5], [5, 6], [6, 6], [6, 3], [3, 3], [3, 2], [2, 2]]
    # # hole2=[[18,11],[26,11],[26,5],[18,5]]
    #
    # #obs1 = [(2, 5), (5, 5), (5, 6), (6, 6), (6, 3), (3, 3), (3, 2), (2, 2)]
    # # obs2=[(18,5),(26,5),(26,11),(18,11)]
    # maxc_x = 8  # number of colmuns in x
    # maxr_y = 8  # number of rows in y
    # freecell = []
    # ltable = {}
    # lvispoly = {}
    # Check for empty lines
    file_handler = open("input_file", "r")
    raw_data = file_handler.read()
    raw_data = raw_data.split("\n")

    #------------------------------------------
    temp = parse_input_line(raw_data[0])
    LineList=temp
    env=[tuple(i) for i in LineList]
    #_-----------------------------------------



    #Point is from the geometry library and it takes the values of temp
    boundary=[point(i[0], i[1]) for i in temp]

    temp = parse_input_line(raw_data[len(raw_data)-1])
    source = point(temp[0][0], temp[0][1])
    dest = point(temp[1][0], temp[1][1])

    # Extract obstacles which is made up of three lists
    obstacles = []
    for i in raw_data[1:len(raw_data) - 1]:
        obstacles.append(parse_input_line(i))

    #-------------------------------------------------
    hole1=parse_input_line(raw_data[1])
    obs1= [tuple(i) for i in hole1]
    print(hole1)

    hole2=parse_input_line(raw_data[2])
    obs2= [tuple(i) for i in hole2]
    print(hole2)

    hole3=parse_input_line(raw_data[3])
    obs3= [tuple(i) for i in hole3]
    print(hole3)
    #--------------------------------------------------

    sorted_vertices = []
    for index, i in enumerate(obstacles):
        for j in i:
            j.append(index)
            sorted_vertices.append(j)


    sorted_vertices.sort(key=lambda x: x[0])

    draw_problem()
    plt.show()

    new_sorted_vertices = []
    for i in sorted_vertices:
        temp = point(i[0], i[1], i[2])
        new_sorted_vertices.append(temp)

    new_obstacles = []
    for index, i in enumerate(obstacles):
        temp_obs = []
        for j in i:
            temp = point(j[0], j[1], index)
            # print(temp)
            temp_obs.append(temp)
        new_obstacles.append(temp_obs)

    # Find vertical lines
    open_line_segments = []

    y_limit_lower = boundary[0].y
    y_limit_upper = boundary[2].y

    for pt in new_sorted_vertices:
        curr_line_segment = [point(pt.x, y_limit_lower), point(pt.x, y_limit_upper)]
        lower_obs_pt = curr_line_segment[0]
        upper_obs_pt = curr_line_segment[1]
        upper_gone = False
        lower_gone = False
        break_now = False

        # Find intersection points with the vertical proposed lines. the intersection function returns false if segments are same, so no need to worry about same segment checking
        for index, obs in enumerate(new_obstacles):
            # Add the first point again for the last line segment of a polygon.

            obs.append(obs[0])
            for vertex_index in range(len(obs) - 1):
                res = segment_intersection(curr_line_segment[0], curr_line_segment[1], obs[vertex_index],
                                           obs[vertex_index + 1])
                if (res != -1):
                    if (index == pt.obstacle):
                        if pt.equals(res) == False:
                            if (res.y > pt.y):
                                upper_gone = True
                            elif (res.y < pt.y):
                                lower_gone = True
                    else:
                        if pt.equals(res) == False:
                            if (upper_gone is False):
                                if ((res.y > pt.y) and res.y < (upper_obs_pt.y)):
                                    upper_obs_pt = res
                            if (lower_gone is False):
                                if ((res.y < pt.y) and (res.y > lower_obs_pt.y)):
                                    lower_obs_pt = res
                if (upper_gone is True and lower_gone is True):
                    break_now = True

            # No need to check for current point anymore...completely blocked
            if (break_now is True):
                break

        # Draw the vertical cell lines
        if (lower_gone is False):
            plt.plot([lower_obs_pt.x, pt.x], [lower_obs_pt.y, pt.y])

        if (upper_gone is False):
            plt.plot([pt.x, upper_obs_pt.x], [pt.y, upper_obs_pt.y])

        # Add to the global segment list
        if (lower_gone and upper_gone):
            open_line_segments.append([None, None])
        elif (lower_gone):
            open_line_segments.append([None, upper_obs_pt])
        elif (upper_gone):
            open_line_segments.append([lower_obs_pt, None])
        else:
            open_line_segments.append([lower_obs_pt, upper_obs_pt])

    cells = []

    for index1 in range(len(open_line_segments)):
        curr_segment = open_line_segments[index1]
        curr_vertex = new_sorted_vertices[index1]
        break_now = False
        done = [False, False, True]
        if (curr_segment[0] is None):
            done[0] = True
        if (curr_segment[1] is None):
            done[1] = True
        if (curr_segment[1] is None and open_line_segments[index1][0] is None):
            done[2] = False

        for index2 in range(index1 + 1, len(open_line_segments)):
            next_segment = open_line_segments[index2]
            next_vertex = new_sorted_vertices[index2]

            double_index1 = -2
            double_index2 = -2
            lines_to_check = []
            trapezoids = []
            double_check = False

            if (next_segment[0] is not None and next_segment[1] is not None):
                double_check = True

            if (done[0] is False):
                if (double_check):
                    double_index1 = len(lines_to_check)
                    lines_to_check.append(
                        [centroidver([curr_segment[0], curr_vertex]), centroidver([next_segment[0], next_vertex]), 0])
                    lines_to_check.append(
                        [centroidver([curr_segment[0], curr_vertex]), centroidver([next_segment[1], next_vertex]), 0])
                    trapezoids.append([curr_segment[0], next_segment[0], next_vertex, curr_vertex])
                    trapezoids.append([curr_segment[0], next_vertex, next_segment[1], curr_vertex])
                elif (next_segment[0] is not None):
                    lines_to_check.append(
                        [centroidver([curr_segment[0], curr_vertex]), centroidver([next_segment[0], next_vertex]), 0])
                    trapezoids.append([curr_segment[0], next_segment[0], next_vertex, curr_vertex])
                elif (next_segment[1] is not None):
                    lines_to_check.append(
                        [centroidver([curr_segment[0], curr_vertex]), centroidver([next_segment[1], next_vertex]), 0])
                    trapezoids.append([curr_segment[0], next_vertex, next_segment[1], curr_vertex])
                else:
                    lines_to_check.append([centroidver([curr_segment[0], curr_vertex]), next_vertex, 0])
                    trapezoids.append([curr_segment[0], next_vertex, curr_vertex])

            if (done[1] is False):
                if (double_check):
                    double_index2 = len(lines_to_check)
                    lines_to_check.append(
                        [centroidver([curr_segment[1], curr_vertex]), centroidver([next_segment[0], next_vertex]), 1])
                    lines_to_check.append(
                        [centroidver([curr_segment[1], curr_vertex]), centroidver([next_segment[1], next_vertex]), 1])
                    trapezoids.append([curr_vertex, next_segment[0], next_vertex,
                                       point(curr_segment[1].x, curr_segment[1].y, curr_segment[1].obstacle, 34)])
                    trapezoids.append([curr_vertex, next_vertex, next_segment[1], curr_segment[1]])
                elif (next_segment[1] is not None):
                    lines_to_check.append(
                        [centroidver([curr_segment[1], curr_vertex]), centroidver([next_segment[1], next_vertex]), 1])
                    trapezoids.append([curr_vertex, next_vertex, next_segment[1], curr_segment[1]])
                elif (next_segment[0] is not None):
                    lines_to_check.append(
                        [centroidver([curr_segment[1], curr_vertex]), centroidver([next_segment[0], next_vertex]), 1])
                    trapezoids.append([curr_vertex, next_segment[0], next_vertex, curr_segment[1]])
                else:
                    lines_to_check.append([centroidver([curr_segment[1], curr_vertex]), next_vertex, 1])
                    trapezoids.append([curr_vertex, next_vertex, curr_segment[1]])

            if (done[2] is False):
                if (double_check):
                    double_index = len(lines_to_check)
                    lines_to_check.append([curr_vertex, centroidver([next_segment[0], next_vertex]), 2])
                    trapezoids.append([curr_vertex, next_segment[0], next_vertex])
                    lines_to_check.append([curr_vertex, centroidver([next_segment[1], next_vertex]), 2])
                    trapezoids.append([curr_vertex, next_vertex, next_segment[1]])
                elif (next_segment[0] is not None):
                    lines_to_check.append([curr_vertex, centroidver([next_segment[0], next_vertex]), 2])
                    trapezoids.append([curr_vertex, next_segment[0], next_vertex])
                elif (next_segment[1] is not None):
                    lines_to_check.append([curr_vertex, centroidver([next_segment[1], next_vertex]), 2])
                    trapezoids.append([curr_vertex, next_vertex, next_segment[1]])
                # Will this ever occur though??
                else:
                    lines_to_check.append([curr_vertex, next_vertex, 2])
                    trapezoids.append([curr_vertex, next_vertex])

            temp_to_remove = []
            for index5, q in enumerate(lines_to_check):
                ok = [True, True, True]
                for index3, obs in enumerate(new_obstacles):
                    # Add the last line to make closed polygon
                    obs.append(obs[0])
                    for index4 in range(len(obs) - 1):
                        if (segment_intersection(q[0], q[1], obs[index4], obs[index4 + 1]) != -1):
                            ok[q[2]] = False
                            if (index5 not in temp_to_remove):
                                temp_to_remove.append(index5)

                if (ok[q[2]] is True):
                    done[q[2]] = True

            for i in range(len(lines_to_check)):
                if i not in temp_to_remove:
                    cells.append(trapezoids[i])

            if (done[0] == True and done[1] == True and done[2] == True):
                break

    to_draw = []
    for i in cells:
        i.append(i[0])
        to_draw.append(i)

    # Merge overlapping Polygons
    quad_cells = [i for i in cells if len(i) > 3]
    tri_cells = [i for i in cells if len(i) == 3]
    others = [i for i in cells if len(i) < 3]
    quads_to_remove = []
    quads_to_add = []

    quads_to_remove = []
    quads_to_add = []
    for index_cell in range(len(quad_cells)):
        for index_cell2, cell in enumerate(quad_cells):
            if (index_cell != index_cell2):
                if (quad_cells[index_cell][0].x == cell[0].x and quad_cells[index_cell][1].x == cell[1].x):
                    temp1 = list(quad_cells[index_cell])
                    temp1.append(temp1[0])
                    temp2 = list(cell)
                    temp2.append(temp2[0])
                    area1 = polygon_area(temp1, 4)
                    area2 = polygon_area(temp2, 4)
                    new_quad = []

                    new_quad.append(point(temp1[0].x, min(temp1[0].y, temp2[0].y)))
                    new_quad.append(point(temp1[1].x, min(temp1[1].y, temp2[1].y)))
                    new_quad.append(point(temp1[1].x, max(temp1[2].y, temp2[2].y)))
                    new_quad.append(point(temp1[0].x, max(temp1[3].y, temp2[3].y)))
                    new_quad.append(point(temp1[0].x, min(temp1[0].y, temp2[0].y)))
                    area3 = polygon_area(new_quad, 4)
                    if (area1 + area2 >= area3):
                        # merge
                        quads_to_remove.append(index_cell)
                        quads_to_remove.append(index_cell2)

                        quads_to_add.append(new_quad)

    quads_to_remove = list(set(quads_to_remove))
    for index in sorted(quads_to_remove, reverse=True):
        del quad_cells[index]

    for i in quads_to_add:
        quad_cells.append(i)

    # Remove duplicates
    to_remove = []
    for index1 in range(len(quad_cells)):
        for index2 in range(index1 + 1, len(quad_cells)):
            duplicate = True
            for k, m in zip(quad_cells[index1], quad_cells[index2]):
                if k.equals(m) is False:
                    duplicate = False
                    break
            if (duplicate is True):
                if index2 not in to_remove:
                    to_remove.append(index2)

    for index in sorted(to_remove, reverse=True):
        del quad_cells[index]

    # One more pass to remove extra quads generated because of cross - segments
    quads_to_remove = []
    for index1 in range(len(quad_cells)):
        for index2 in range(len(quad_cells)):
            if (index1 != index2 and quad_cells[index1][0].x == quad_cells[index2][0].x and quad_cells[index1][1].x ==
                    quad_cells[index2][1].x):

                if ((quad_cells[index1][0].y <= quad_cells[index2][0].y) and (
                        quad_cells[index1][1].y <= quad_cells[index2][1].y)
                        and (quad_cells[index1][2].y >= quad_cells[index2][2].y) and (
                                quad_cells[index1][3].y >= quad_cells[index2][3].y)):
                    quads_to_remove.append(index2)

    quads_to_remove = list(set(quads_to_remove))
    for index in sorted(quads_to_remove, reverse=True):
        del quad_cells[index]

    # Add boundary lines
    if (boundary[0].x != new_sorted_vertices[0].x):
        quad_cells.append([boundary[0], point(new_sorted_vertices[0].x, y_limit_lower),
                           point(new_sorted_vertices[0].x, y_limit_upper), boundary[3]])
    if (boundary[1].x != new_sorted_vertices[len(new_sorted_vertices) - 1].x):
        quad_cells.append(
            [point(new_sorted_vertices[len(new_sorted_vertices) - 1].x, y_limit_lower), boundary[1], boundary[2],
             point(new_sorted_vertices[len(new_sorted_vertices) - 1].x, y_limit_upper)])

    point_x_cord = []

    point_y_cord = []

    # Plot final cells

    to_draw = quad_cells + tri_cells + others

    polycentroidlist = []
    polyheight = []

    for i in quad_cells:
        x = [j.x for j in i]
        y = [j.y for j in i]
        plot_polygon_decomposition(x, y)

    # Dictionary for multipolygon
    m_polygons = {}

    # Loop throgh Obstacles and create polygons
    for index, i in enumerate(new_obstacles):
        x = [j.x for j in i]
        y = [j.y for j in i]

        m_polygons[index] = Polygon(np.array([x, y]).T)
        #print(index, len(x), len(y))

    # loop through the dictionaly and create a multipolygon
    m_polygon_obstacles = MultiPolygon([m_polygons[i] for i in m_polygons])

    maxc_x = 200  # number of colmuns in x
    maxr_y = 300 # number of rows in y

    freecell = []
    ltable = {}
    lvispoly = {}

    graph_points=np.array([point_x_cord,point_y_cord]).T


    directions = {31545: [], 45135: [], 135225: [], 225315: []}
    directions_area = {31545: 0, 45135: 0, 135225: 0, 225315: 0}
    #student_id = [1, 2, 3, 4, 5]
    visibility_polygons = {}
    visibility_polygons_area={}
    graph_points_list = graph_points.tolist()
    # print(graph_points_list)area



    for i in range(0, len(graph_points_list)):
        visibility_polygons.update({graph_points_list.index(graph_points_list[i]): dict(directions)})
        visibility_polygons_area.update({graph_points_list.index(graph_points_list[i]): dict(directions_area)})

    graph_vertex_index=[]
    for index, position in enumerate(graph_points.tolist()):
        x,y=position
        print ("x,y", x,y)
        graph_vertex_index.append((index))
        if not insidenvironment(x, y) or insideobstacle(x, y) :
            continue
        else:
            freecell.append(i)
            np.savetxt('gridnewguards', (x, y), fmt='%5.1f')
            unitA = Circle((x, y), .2, facecolor='none', fill=True, color='blue', edgecolor=(0, 0.8, 0.8), linewidth=2,
                           alpha=0.5)
            ccellno = int(midcellnumber(y, x))

            circlepoly = Point(x, y).buffer(2)
            circleview = Polygon(env).intersection(circlepoly)
            #sectorpoly = sector(Point(x, y), 225, 315, 40)

            #315-45
            sectorpolyd1 = sector(Point(x, y), 45, 135, 65)
            #45-135
            sectorpolyd2 = sector(Point(x, y), 315, 45, 65)
            #135-225
            sectorpolyd3 = sector(Point(x, y), 225, 315, 65)
            #225-315
            sectorpolyd4 = sector(Point(x, y), 135, 225, 65)

            circlecoords = list(circleview.exterior.coords)
            polygon1 = genVisibMatrix(0)
            vispolygon = genVisibMatrix(0)

            points = []
            for i in range(len(polygon1)):
                points.append(Point(polygon1[i][0], polygon1[i][1]))
            polygon1 = Polygon([[p.x, p.y] for p in points])

            intersectingpoly = polygon1.intersection(circlepoly)
            circleview = polygon1.intersection(circlepoly)
            #sectorview = polygon1.intersection(sectorpoly)

            sectorviewd1 = polygon1.intersection(sectorpolyd1)
            sectorviewd2 = polygon1.intersection(sectorpolyd2)
            sectorviewd3 = polygon1.intersection(sectorpolyd3)
            sectorviewd4 = polygon1.intersection(sectorpolyd4)

            circlecoords = list(circleview.exterior.coords)
            #sectorcoords = list(sectorview.exterior.coords)

            sectorcoordsd1 = list(sectorviewd1.exterior.coords)
            sectorcoordsd2 = list(sectorviewd2.exterior.coords)
            sectorcoordsd3 = list(sectorviewd3.exterior.coords)
            sectorcoordsd4 = list(sectorviewd4.exterior.coords)

            intspoints = np.array(intersectingpoly)
            #print('Intersect', intspoints)
            ix, iy = intersectingpoly.exterior.xy

            intspoly = []
            for i in range(len(ix)):
                intspoly.append([int(ix[i]), int(iy[i])])
            lvispoly.update({ccellno: intspoly})

            visibility_polygons[index][31545]= sectorcoordsd1
            visibility_polygons[index][45135] = sectorcoordsd2
            visibility_polygons[index][135225] = sectorcoordsd3
            visibility_polygons[index][225315] = sectorcoordsd4

            #print(visibility_polygons.values())
            fig = figure(figsize=(18, 16))
            ax = fig.add_subplot(111, aspect='equal', xlabel="S", ylabel="t")

            if (len(intspoly) > 1):
                # ax.add_patch(MPolygon(intspoly, closed=True, fill=True, color='g', linewidth=0))
                # ax.add_patch(MPolygon(vispolygon, closed=True, fill=True, color='g', linewidth=0))
                #ax.add_patch(MPolygon(sectorcoords, closed=True, fill=True, color='r', linewidth=0))

                ax.add_patch(MPolygon(sectorcoordsd1, closed=True, fill=True, color='r', linewidth=0))
                ax.add_patch(MPolygon(sectorcoordsd2, closed=True, fill=True, color='r', linewidth=0))
                ax.add_patch(MPolygon(sectorcoordsd3, closed=True, fill=True, color='r', linewidth=0))
                ax.add_patch(MPolygon(sectorcoordsd4, closed=True, fill=True, color='r', linewidth=0))


            ax.add_patch(MPolygon(LineList, closed=True, fill=False, color='black', label='line 1', linewidth=3))
            ax.add_patch(MPolygon(hole1, closed=True, fill=True, color='black', label='line 1', linewidth=2))
            ax.add_patch(MPolygon(hole2, closed=True, fill=True, color='black', label='line 1', linewidth=2))
            ax.add_patch(MPolygon(hole3, closed=True, fill=True, color='black', label='line 1', linewidth=2))

            ax.add_patch(unitA)
            ax.set_xlim(-2, 300)
            ax.set_ylim(-2, 200)

            #timestr = time.strftime("%Y%m%d-%H%M%S")
            #pylab.savefig('output/'+timestr+'No_'+str(l)+".jpeg", bbox_inches='tight')
            show()

            maxpolyx = -10
            maxpolyy = -10
            minpolyx = 1000
            minpolyy = 1000

            for j in range(0, len(intspoly)):
                if (intspoly[j][0] > maxpolyx):
                    maxpolyx = intspoly[j][0]
                if (intspoly[j][1] > maxpolyy):
                    maxpolyy = intspoly[j][1]
                if (intspoly[j][0] < minpolyx):
                    minpolyx = intspoly[j][0]
                if (intspoly[j][1] < minpolyy):
                    minpolyy = intspoly[j][1]
                    # print floor(maxpolyx), ceil(int(minpolyx)),floor(maxpolyy),ceil(int(minpolyy))
            maxpolyx = floor(maxpolyx)
            minpolyx = ceil(int(minpolyx))
            maxpolyy = floor(maxpolyy)
            minpolyy = ceil(int(minpolyy))

            # # instpoly=Polygon(intspoly)
            # # print instpoly
            # ltable.setdefault(ccellno, [])
            # for j in range(int(minpolyx), int(maxpolyx) + 1):
            #     for k in range(int(minpolyy), int(maxpolyy) + 1):
            #         # print j, k
            #         # print intersectingpoly.intersects(Point(1,1))
            #         # print point_in_polyen(0,0, intspoly), point_in_polyen(1,0, intspoly), point_in_polyen(1,1, intspoly), pointOnPolygon([0,1], intspoly)
            #         if intersectingpoly.intersects(Point(j, k)) and intersectingpoly.intersects(
            #                 Point(j + 1, k)) and intersectingpoly.intersects(
            #                 Point(j + 1, k + 1)) and intersectingpoly.intersects(Point(j, k + 1)):
            #             # print cellnumber(j, k), j, k
            #             viscell = midcellnumber(k, j)
            #             # if viscell==204:
            #             #    print viscell, cellposition(viscell),j,k
            #             ltable[ccellno].append(viscell)


    #print(graph_vertex_index)
    F = graph_vertex_index
    lvalues = list(visibility_polygons.values())

    #visibility_polygons_area ={}


    visibility_polygons_area_list = []
    for index, value in enumerate(lvalues):
        #print(index, Polygon(value[31545]).area)
        visibility_polygons_area[index][31545] = Polygon(value[31545]).area
        visibility_polygons_area_list.append(Polygon(value[31545]).area)
        visibility_polygons_area[index][45135] = Polygon(value[45135]).area
        visibility_polygons_area_list.append(Polygon(value[45135]).area)
        visibility_polygons_area[index][135225] = Polygon(value[135225]).area
        visibility_polygons_area_list.append(Polygon(value[135225]).area)
        visibility_polygons_area[index][225315] = Polygon(value[225315]).area
        visibility_polygons_area_list.append(Polygon(value[225315]).area)
    # print(len(lvalues))
    # print(lvalues[0][31545])
    # print(lvalues[0][45135])
    # print(lvalues[0][135225])
    # print(lvalues[0][225315])

    visibility_polygons_area_sorted_list = sorted(visibility_polygons_area_list,reverse=True)
    print(visibility_polygons_area_sorted_list)
    # free_env_space=Polygon(env).exterior.difference(m_polygon_obstacles.exterior)
    print(Polygon(obs1).area)
    print(Polygon(env).area)
    free_env_space = Polygon(env).difference(Polygon(obs1))
    free_env_space2 = free_env_space.difference(Polygon(obs2))
    free_env_space3 = free_env_space2.difference(Polygon(obs3))

    free = Polygon(env).difference(Polygon(env))
    if free.is_empty:
        print("Polygon is Empty")
    else:
        print("Polygon is not empty")
    print(free_env_space2)
    print(free_env_space3)
    print(extract_poly_coords(free_env_space3)['exterior_coords'],extract_poly_coords(free_env_space3)['interior_coords'])

    fig, ax = plt.subplots()

    # Create Polygon
    exterior = extract_poly_coords(free_env_space3)['exterior_coords']
    #exterior_int = [(int(x),int(y)) for x,y in exterior]
    #print(exterior_int)
    interiors = extract_poly_coords(free_env_space3)['interior_coords']
    #print(interiors)
    #interiors_int = [[(int(x),int(y)) for x,y in interiors]]
    #print((interiors_int))
    # exterior = [(20, 20), (80, 20), (50, 70)]
    # interiors = [[(40, 25), (40, 30), (45, 30), (45, 25)],
    #              [(50, 35), (50, 40), (55, 40)]]
    poly = Polygon(exterior, holes=interiors)

    plot_polygon(ax, poly, facecolor="blue", edgecolor="black")

    ax.axis([0, 300, 0, 200])
    plt.show()

    #free_space_poly = gpd.GeoSeries([free_env_space])
    #print(free_env_space)

    #print(free_space_poly)
    U = set()

    # x, y = free_env_space.exterior.xy
    # plt.plot(x, y, color='black', linewidth=2)
    # for interior in free_space_poly.interiors:
    #     x,y=list(interior).xy
    #     plt.plot(x, y, color='red', linewidth=2)

    plt.show()
    # C = set()
    # cellsequence = []
    # while U is not None:
    #     maxl = -10001
    #     for k in range(len(F)):
    #         S = set(ltable.get(int(F[k])))
    #         # print S
    #         L = U.intersection(S)
    #         # print U.intersection(S)
    #         # print len(L)
    #         if len(L) > maxl:
    #             maxl = len(L)
    #             maxkey = int(F[k])
    #
    #     print('mmaxl', maxl, maxkey)
    #     MS = set(ltable.get(int(maxkey)))
    #     U = U - MS
    #     print(len(U), len(F), MS)
    #
    #     if maxl == 0:
    #         break
    #     F.remove(int(maxkey))
    #
    #     # print len(ltable.get(int(maxkey)))
    #     # i=i+1
    #     C.update(list(MS))
    #     print(len(C))
    #     cellsequence.append(maxkey)
    # print(len(U))
    # print(C)
    # print(len(C))
    # print(U - C)
    # print('cell sequence', cellsequence)
    # for g in range(len(lvalues)):
    #     U.update(lvalues[g])
    # G = set()

    # fig = figure(figsize=(18, 16))
    # ax = fig.add_subplot(111, aspect='equal', xlabel="S", ylabel="t")
    #
    # ax.add_patch(MPolygon(LineList, closed=True, fill=False, color='black', label='line 1', linewidth=3))
    # ax.add_patch(MPolygon(hole1, closed=True, fill=True, color='black', label='line 1', linewidth=2))
    # ax.add_patch(MPolygon(hole2, closed=True, fill=True, color='black', label='line 1', linewidth=2))
    # ax.add_patch(MPolygon(hole3, closed=True, fill=True, color='black', label='line 1', linewidth=2))
    # #ax.add_patch(poly)

    # ax.add_patch(unitA)
    # ax.set_xlim(-2, 300)
    # ax.set_ylim(-2, 200)

    # timestr = time.strftime("%Y%m%d-%H%M%S")
    # pylab.savefig('output/'+timestr+'No_'+str(l)+".jpeg", bbox_inches='tight')
    #show()

if __name__ == "__main__":main()
