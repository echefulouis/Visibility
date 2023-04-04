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
from scipy.spatial.distance import pdist, squareform
from networkx.algorithms import approximation as approx

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


def get_hamiltonian_path_with_path_legnth_as_weights(vertex_cords_dict,adj_matrix_distance):
    vertex_cord = list(vertex_cords_dict.values())
    vertex_node = list(vertex_cords_dict.keys())

    F = nx.Graph()
    for i in range(len(vertex_node)):
        F.add_node(vertex_node[i], pos=vertex_cord[i])

    for i in range(len(vertex_node)):
        for j in range(len(vertex_node)):
            if i !=j:
                F.add_edge(vertex_node[i], vertex_node[j],
                           weight=((adj_matrix_distance[i][j])))

    hamiltonian_path = nx.approximation.traveling_salesman_problem(F,cycle=False)
    hamiltonian_path = approx.greedy_tsp(F, source=hamiltonian_path[0])
    #print(hamiltonian_path)
    # Q = nx.Graph()
    # for i in hamiltonian_path:
    #     Q.add_node(i, pos=vertex_cords_dict[i])
    #
    # for i in range(len(hamiltonian_path)-1):
    #     Q.add_edge(hamiltonian_path[i],hamiltonian_path[i+1],weight = F[hamiltonian_path[i]][hamiltonian_path[i+1]]["weight"])

    #hamiltonian_path = nx.approximation.traveling_salesman_problem(F,nodes=first_path,weight='weight',cycle=False)
    #hamiltonian_path = approx.threshold_accepting_tsp(F , "greedy")
    #hamiltonian_path = approx.simulated_annealing_tsp(F, "greedy")
    #hamiltonian_path  = nx.algorithms.tournament.hamiltonian_path(F)

    #print("Hamiltonian Path is: {}".format(hamiltonian_path))


    # Calculate total distance
    total_distance = 0
    hamiltonian_edge_distance = {}
    for i in range(len(hamiltonian_path) - 1):
        #print("distance from node {} to {} is {}".format(hamiltonian_path[i],hamiltonian_path[i+1],F[hamiltonian_path[i]][hamiltonian_path[i + 1]]['weight']))
        hamiltonian_edge_distance[hamiltonian_path[i]]=F[hamiltonian_path[i]][hamiltonian_path[i + 1]]['weight']
        total_distance += F[hamiltonian_path[i]][hamiltonian_path[i + 1]]['weight']

    return F.copy(),hamiltonian_path ,hamiltonian_edge_distance

def get_hamiltonian_path_with_ratio_as_weights(vertex_cords_dict,adj_matrix_distance,max_sector_area_dict):
    vertex_cord = list(vertex_cords_dict.values())
    vertex_node = list(vertex_cords_dict.keys())

    F = nx.Graph()
    for i in range(len(vertex_node)):
        F.add_node(vertex_node[i], pos=vertex_cord[i])

    for i in range(len(vertex_node)):
        for j in range(len(vertex_node)):
            if i !=j:
                F.add_edge(vertex_node[i], vertex_node[j], weight=((adj_matrix_distance[i][j] /max_sector_area_dict[vertex_node[i]])*100))


    hamiltonian_path = nx.approximation.traveling_salesman_problem(F,cycle=False)
    hamiltonian_path = approx.greedy_tsp(F, source=hamiltonian_path[0])
    #print(hamiltonian_path)

    total_distance = 0
    hamiltonian_edge_distance = {}
    for i in range(len(hamiltonian_path) - 1):
        #print("distance from node {} to {} is {}".format(hamiltonian_path[i],hamiltonian_path[i+1],F[hamiltonian_path[i]][hamiltonian_path[i + 1]]['weight']))
        hamiltonian_edge_distance[hamiltonian_path[i]]=F[hamiltonian_path[i]][hamiltonian_path[i + 1]]['weight']
        total_distance += F[hamiltonian_path[i]][hamiltonian_path[i + 1]]['weight']

    return F.copy(),hamiltonian_path ,hamiltonian_edge_distance



def draw_hamiltonian_cycle(h_path, tsp_cord):
    Z = nx.Graph()
    vertexsequence_edges = []
    vertexsequence = h_path
    vertexsequence_cords = tsp_cord

    for i in range(len(vertexsequence) - 1):
        vertexsequence_edges.append(tuple([vertexsequence[i], vertexsequence[i + 1]]))

    for i in vertexsequence:
        Z.add_node(i, pos=vertexsequence_cords[i])

    Z.add_edges_from(vertexsequence_edges)
    pos = nx.get_node_attributes(Z, 'pos')
    nx.draw(Z, pos, with_labels=True)
    show()


def draw_hamiltonian_circles_on_main_graph(nx_graph,nx_coords,h_path_coords,h_path,max_sec_coords):
    for index,num in enumerate(h_path):
        if index+1 == len(h_path):
            pass
        else:
            x, y = h_path_coords[num][0],h_path_coords[num][1]
            a, b = h_path_coords[h_path[index+1]][0], h_path_coords[h_path[index+1]][1]
            unitA = Circle((x, y), 5, facecolor='none', fill=True, color='blue', edgecolor=(0, 0.8, 0.8), linewidth=2,
                           alpha=0.5)
            unitB = Circle((a, b), 5, facecolor='none', fill=True, color='magenta', edgecolor=(0, 0.8, 0.8), linewidth=2,
                           alpha=0.5)
            fig = figure(figsize=(18, 16))
            ax = fig.add_subplot(111, aspect='equal', xlabel="S", ylabel="t")
            ax.add_patch(MPolygon(max_sec_coords[num], closed=True, fill=True, color='r', linewidth=0))
            ax.add_patch(MPolygon(LineList, closed=True, fill=False, color='black', label='line 1', linewidth=3))
            ax.add_patch(MPolygon(hole1, closed=True, fill=True, color='gray', label='line 1', linewidth=2))
            ax.add_patch(MPolygon(hole2, closed=True, fill=True, color='gray', label='line 1', linewidth=2))
            ax.add_patch(MPolygon(hole3, closed=True, fill=True, color='gray', label='line 1', linewidth=2))

            ax.add_patch(unitA)
            ax.set_xlim(-2, 300)
            ax.set_ylim(-2, 200)

            nx.draw(nx_graph, nx_coords, with_labels=True)
            # Find the shortest path between the nodes

            path = nx.shortest_path(nx_graph, source=h_path[index], target=h_path[index+1])
            path_edges = list(zip(path, path[1:]))

            nx.draw(nx_graph, nx_coords,edgelist = path_edges, edge_color = 'y', width = 2.0,with_labels=True)
            ax.add_patch(unitB)
            show()

def assign_path_to_robots(edges_dist,no_of_robots):
    #Average distance
    W = sum(list(edges_dist.values())) / no_of_robots
    #print("average distance "+str(W))

    path_for_robot = {}

    #Store the index of node the previous robot stopped and start of new robot
    start = 0
    for robot in range(1, no_of_robots+1):
        #Store the total distance a robot has
        path_list = []
        total_distance = 0
        for index_node in range(start, len(edges_dist.keys())):
            if  total_distance <W:
                path_list.extend([list(edges_dist.keys())[index_node]])
                total_distance = total_distance + list(edges_dist.values())[index_node]
                start = start +1

        path_for_robot[robot]=path_list

    return path_for_robot





def get_total_path_length(h_edge_distance):
    return sum(list(h_edge_distance.values()))

def plot_polygon_decomposition(x, y):
    # Create a polygon with the coordinates of the vertices
    polycoords = [(x[0], y[0]), (x[1], y[1]), (x[2], y[2]), (x[3], y[3])]
    # coords = [(0, 0), (0, 1), (1, 1), (1, 0), (0, 0)]
    polygon = Polygon(polycoords)

    # distance(shortest)
    line1 = LineString([(x[0], y[0]), (x[1], y[1])])
    line2 = LineString([(x[2], y[2]), (x[3], y[3])])

    distance1 = line1.distance(line2)

    # distance(shortest)
    line3 = LineString([(x[0], y[0]), (x[3], y[3])])
    line4 = LineString([(x[1], y[1]), (x[2], y[2])])

    distance2 = line3.distance(line4)

    # Get the minimum bounding rectangle of the polygon
    min_x, min_y, max_x, max_y = polygon.bounds

    # Find the centroid of the polygon
    polycentroid = polygon.centroid

    if (distance2 == max_x - min_x):
        # Calculate the height of the polygon along x-axis
        height = max_x - min_x

        # find area and median
        area = polygon.area
        median = area / height

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

    plt.plot(x, y)

def get_networkx_neighbours(networkx_points,new_obstacles,num_neigbours):
    from sklearn.metrics.pairwise import pairwise_distances
    neigh = NearestNeighbors(n_neighbors=num_neigbours)
    # Fit the model
    neigh.fit(networkx_points)

    # Calculate pairwise distances between points
    distances = pairwise_distances(networkx_points,networkx_points,metric="euclidean")
    # Create weighted adjacency matrix from distances
    adj_matrix = distances

    # Predict the nearest neighbors
    neighbors = neigh.kneighbors_graph(networkx_points, mode='distance')
    G = nx.from_numpy_matrix(neighbors.toarray())
    # Remove self loop

    # Add weights to edges
    for i, j in G.edges():
        G[i][j]['weight'] = adj_matrix[i, j]

    G.remove_edges_from(nx.selfloop_edges(G))
    # print(G.nodes())
    ax = plt.gca()
    for index, i in enumerate(new_obstacles):
        g = [j.x for j in i]
        h = [j.y for j in i]

        polygon_patch = MPolygon(np.array([g, h]).T, closed=True, color="black")
        ax.add_patch(polygon_patch)
    # Dictionary for multipolygon
    polygons = {}

    for index, i in enumerate(new_obstacles):
        x = [j.x for j in i]
        y = [j.y for j in i]

        polygons[index] = Polygon(np.array([x, y]).T)

    # loop through the dictionaly and create a multipolygon
    mp = MultiPolygon([polygons[i] for i in polygons])

    # Convert networkx node lines into shapely lines
    line_nodes = G.edges()
    # print(len(line_nodes))

    # Check for collision between each line and collision object
    for node in line_nodes:
        # print(list(points[node[0]])[0], list(points[node[0]])[1])
        point1 = (list(networkx_points[node[0]])[0], list(networkx_points[node[0]])[1])

        point2 = (list(networkx_points[node[1]])[0], list(networkx_points[node[1]])[1])
        lines = LineString([point1, point2])
        collision = lines.intersects(mp)

        # Check Colllision with Multi-polygon
        if collision:
            G.remove_edge(node[0], node[1])

    nx.draw(G,networkx_points, with_labels=True, node_size=20)
    plt.show()

    return G.copy()

def create_environment():
    global env, obs1, obs2, obs3, maxc_x, maxr_y, boundary, obstacles, source, dest, point_x_cord, point_y_cord, new_obstacles, hole1,hole2,hole3,LineList
    maxc_x = 200  # number of colmuns in x
    maxr_y = 300  # number of rows in y
    #-------------------------------------------
    file_handler = open("input_file", "r")
    raw_data = file_handler.read()
    raw_data = raw_data.split("\n")
    # ------------------------------------------
    temp = parse_input_line(raw_data[0])
    LineList = temp
    env = [tuple(i) for i in LineList]
    # _-----------------------------------------
    # Point is from the geometry library and it takes the values of temp
    boundary = [point(i[0], i[1]) for i in temp]

    temp = parse_input_line(raw_data[len(raw_data) - 1])
    source = point(temp[0][0], temp[0][1])
    dest = point(temp[1][0], temp[1][1])

    # Extract obstacles which is made up of three lists
    obstacles = []
    for i in raw_data[1:len(raw_data) - 1]:
        obstacles.append(parse_input_line(i))
    # -------------------------------------------------
    hole1 = parse_input_line(raw_data[1])
    obs1 = [tuple(i) for i in hole1]

    hole2 = parse_input_line(raw_data[2])
    obs2 = [tuple(i) for i in hole2]

    hole3 = parse_input_line(raw_data[3])
    obs3 = [tuple(i) for i in hole3]
    # --------------------------------------------------
    sorted_vertices = []
    for index, i in enumerate(obstacles):
        for j in i:
            j.append(index)
            sorted_vertices.append(j)

    sorted_vertices.sort(key=lambda x: x[0])

    # draw_problem()
    # plt.show()

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

    to_draw = quad_cells + tri_cells + others

    for i in quad_cells:
        x = [j.x for j in i]
        y = [j.y for j in i]
        plot_polygon_decomposition(x, y)
    #print(len(point_y_cord), len(point_x_cord))
    # Dictionary for multipolygon
    m_polygons = {}

    # Loop throgh Obstacles and create polygons
    for index, i in enumerate(new_obstacles):
        x = [j.x for j in i]
        y = [j.y for j in i]

        m_polygons[index] = Polygon(np.array([x, y]).T)

def main():
    create_environment()
    freecell = []
    ltable = {}
    lvispoly = {}

    graph_points = np.array([point_x_cord, point_y_cord]).T

    directions = {31545: [], 45135: [], 135225: [], 225315: []}
    directions_area = {31545: 0, 45135: 0, 135225: 0, 225315: 0}
    visibility_polygons = {}
    visibility_polygons_area = {}
    graph_points_list = graph_points.tolist()

    for i in range(0, len(graph_points_list)):
        visibility_polygons.update({graph_points_list.index(graph_points_list[i]): dict(directions)})
        visibility_polygons_area.update({graph_points_list.index(graph_points_list[i]): dict(directions_area)})

    #print(visibility_polygons,visibility_polygons_area)

    networkx_points = np.array([point_x_cord, point_y_cord]).T
    R=get_networkx_neighbours(networkx_points,new_obstacles,9)

    print(graph_points_list)
    r=150
    #Try
    #try 500
    for index, position in enumerate(graph_points.tolist()):
        x,y=position
        if not insidenvironment(x, y) or insideobstacle(x, y) :
            continue
        else:
            freecell.append(i)
            np.savetxt('gridnewguards', (x, y), fmt='%5.1f')
            unitA = Circle((x, y), 5, facecolor='none', fill=True, color='blue', edgecolor=(0, 0.8, 0.8), linewidth=2,
                           alpha=0.5)
            ccellno = int(midcellnumber(y, x))
            # 315-45
            sectorpolyd1 = sector(Point(x, y), 45, 135, r)
            # 45-135
            sectorpolyd2 = sector(Point(x, y), 315, 45, r)
            # 135-225
            sectorpolyd3 = sector(Point(x, y), 225, 315, r)
            # 225-315
            sectorpolyd4 = sector(Point(x, y), 135, 225, r)
            polygon1 = genVisibMatrix(0)
            points = []
            for i in range(len(polygon1)):
                points.append(Point(polygon1[i][0], polygon1[i][1]))
            polygon1 = Polygon([[p.x, p.y] for p in points])

            sectorviewd1 = polygon1.intersection(sectorpolyd1)
            sectorviewd2 = polygon1.intersection(sectorpolyd2)
            sectorviewd3 = polygon1.intersection(sectorpolyd3)
            sectorviewd4 = polygon1.intersection(sectorpolyd4)

            sectorcoordsd1 = list(sectorviewd1.exterior.coords)
            sectorcoordsd2 = list(sectorviewd2.exterior.coords)
            sectorcoordsd3 = list(sectorviewd3.exterior.coords)
            sectorcoordsd4 = list(sectorviewd4.exterior.coords)

            visibility_polygons[index][31545] = sectorcoordsd1
            visibility_polygons[index][45135] = sectorcoordsd2
            visibility_polygons[index][135225] = sectorcoordsd3
            visibility_polygons[index][225315] = sectorcoordsd4

            visibility_polygons_area[index][31545] = sectorviewd1.area
            visibility_polygons_area[index][45135] = sectorviewd2.area
            visibility_polygons_area[index][135225] = sectorviewd3.area
            visibility_polygons_area[index][225315] = sectorviewd4.area

            # fig = figure(figsize=(18, 16))
            # ax = fig.add_subplot(111, aspect='equal', xlabel="S", ylabel="t")


            # ax.add_patch(MPolygon(sectorcoordsd1, closed=True, fill=True, color='r', linewidth=0))
            # ax.add_patch(MPolygon(sectorcoordsd2, closed=True, fill=True, color='r', linewidth=0))
            # ax.add_patch(MPolygon(sectorcoordsd3, closed=True, fill=True, color='r', linewidth=0))
            # ax.add_patch(MPolygon(sectorcoordsd4, closed=True, fill=True, color='r', linewidth=0))
            #
            # ax.add_patch(MPolygon(LineList, closed=True, fill=False, color='black', label='line 1', linewidth=3))
            # ax.add_patch(MPolygon(hole1, closed=True, fill=True, color='black', label='line 1', linewidth=2))
            # ax.add_patch(MPolygon(hole2, closed=True, fill=True, color='black', label='line 1', linewidth=2))
            # ax.add_patch(MPolygon(hole3, closed=True, fill=True, color='black', label='line 1', linewidth=2))
            #
            # ax.add_patch(unitA)
            # ax.set_xlim(-2, 300)
            # ax.set_ylim(-2, 200)

            # timestr = time.strftime("%Y%m%d-%H%M%S")
            # pylab.savefig('output/'+timestr+'No_'+str(l)+".jpeg", bbox_inches='tight')
            #show()

            maxpolyx = -10
            maxpolyy = -10
            minpolyx = 1000
            minpolyy = 1000


#--------------------Loop through the environment (free space) and explore it using the visibility---------
    free_env_space1 = Polygon(env).difference(Polygon(obs1))
    free_env_space2 = free_env_space1.difference(Polygon(obs2))
    free_env_space = free_env_space2.difference(Polygon(obs3))

    vertexsequence = []
    directions_keys = list(directions.keys())

    tmp_visibility_polygons = visibility_polygons.copy()
    free_space = free_env_space
    count = 0
    max_sector_visibility_polygon_area={}
    max_sectorcoords_dict = {}
    while not free_space.is_empty:
        count = count + 1
        max_visibility_area = 0
        max_vertex_index=0
        max_visibilty_polygon = 0
        for i in tmp_visibility_polygons.keys():
            value_count = 0
            for value in tmp_visibility_polygons[i].values():
                #max_sector = free_space.intersection(Polygon(value))
                try:
                    max_sector = free_space.intersection(Polygon(value))
                except Exception as e:
                    # If a TopologyException is raised, try buffering and simplifying the polygon
                    if "TopologyException" in str(e):
                        free_space = free_space.buffer(0.001).simplify(0.001)
                        max_sector = free_space.intersection(Polygon(value))
                if (max_sector.area > max_visibility_area):
                    max_visibility_area = max_sector.area
                    ##max intersecting  visibility
                    max_visibilty_polygon = max_sector
                    #max sector visibility
                    max_visibilty_polygon_for_view = value

                    # print(starting_adjacent_vertex_index_for_i + value_count, i)
                    #max_adjacent_vertex = starting_adjacent_vertex_index_for_i + value_count
                    max_value_count= value_count
                    max_vertex_index = i
                value_count = value_count + 1
        max_sector_visibility_polygon_area[max_vertex_index]=max_visibility_area

        max_sectorcoords_dict[max_vertex_index] = max_visibilty_polygon_for_view
        # print(max_visibility_area, max_vertex_index) #Store the max_vis and Index in a dictionary
        # print("-------",free_space.difference(max_visibilty_polygon))
        # #sectorcoords = list(max_visibilty_polygon.exterior.coords)
        # sectorcoords = max_visibilty_polygon_for_view
        # x,y=networkx_points.tolist()[max_vertex_index][0], networkx_points.tolist()[max_vertex_index][1]
        # unitA = Circle((x, y), 5, facecolor='none', fill=True, color='blue', edgecolor=(0, 0.8, 0.8), linewidth=2,  alpha=0.5)
        # fig = figure(figsize=(18, 16))
        # ax = fig.add_subplot(111, aspect='equal', xlabel="S", ylabel="t")
        # ax.add_patch(MPolygon(sectorcoords, closed=True, fill=True, color='r', linewidth=0))
        # ax.add_patch(MPolygon(LineList, closed=True, fill=False, color='black', label='line 1', linewidth=3))
        # ax.add_patch(MPolygon(hole1, closed=True, fill=True, color='gray', label='line 1', linewidth=2))
        # ax.add_patch(MPolygon(hole2, closed=True, fill=True, color='gray', label='line 1', linewidth=2))
        # ax.add_patch(MPolygon(hole3, closed=True, fill=True, color='gray', label='line 1', linewidth=2))
        #
        # ax.add_patch(unitA)
        # ax.set_xlim(-2, 300)
        # ax.set_ylim(-2, 200)
        # show()

        #if  max_centroid_index  not in cellsequence:
        vertexsequence.append(max_vertex_index)

        tmp_visibility_polygons[max_vertex_index][directions_keys[max_value_count]] = []
        visibility_polygons_area[max_vertex_index][directions_keys[max_value_count]] = 0
        #print (visibility_polygons_area)
        #print("The number of iteration is ",count)
        #print("The previous area of the free space is ", free_space.area)
        free_space= free_space.difference(max_visibilty_polygon)
        #print(free_space)
        #print("The maximum area deducted is ",max_visibility_area)
        #print("The current area of the free space is ",free_space.area)

    print(vertexsequence)



    #------------------------------------------Create a adjacent matrix of the shortest edge distance using main graph between two points-------------------------
    vertexsequence_cords = {}
    for i in vertexsequence:
        vertexsequence_cords[i]=tuple(graph_points.tolist()[i])

    mat = np.zeros((len(vertexsequence), len(vertexsequence)))
    mat_dist = np.zeros((len(vertexsequence), len(vertexsequence)))

    for i in range(len(vertexsequence)):
        for j in range(len(vertexsequence)):
            if i != j:
                path = nx.shortest_path(R, source=vertexsequence[i], target=vertexsequence[j])
                path_length = 0
                for k in range(len(path) - 1):
                    path_length = path_length + R[path[k]][path[k + 1]]['weight']
                if not path:
                    pass
                else:
                    mat[i][j] =1
                    mat_dist[i][j] = path_length


    #--------------------------------------------Hamiltonianpath and edge weights--------------------
    ham_graph1,ham_path1, edge_distance1 = get_hamiltonian_path_with_ratio_as_weights(vertexsequence_cords,mat_dist,max_sector_visibility_polygon_area)
    draw_hamiltonian_cycle(ham_path1, vertexsequence_cords)

    ham_graph2,ham_path2, edge_distance2 = get_hamiltonian_path_with_path_legnth_as_weights(vertexsequence_cords, mat_dist)
    draw_hamiltonian_cycle(ham_path2, vertexsequence_cords)

    print(ham_path1)
    #draw_hamiltonian_circles_on_main_graph(R,networkx_points,vertexsequence_cords,ham_path1,max_sectorcoords_dict)

    print(edge_distance1)
    print(edge_distance2)

    #total_distance = get_total_path_length(edge_distance1)

    # nuumber_of_robots = 4
    #
    # path_for_each_robot = assign_path_to_robots(edge_distance,nuumber_of_robots)
    #
    # print(path_for_each_robot)
    #
    # #print(max_sector_visibility_polygon_area)
    #
    # #print(mat)
    # #print(mat_dist)

if __name__ == "__main__":
    main()