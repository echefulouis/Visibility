from .helpers.geometry import *
import matplotlib.pyplot as plt
from shapely.geometry import *
import numpy as np

class Point:
    def __init__(self, x, y, obstacle=None):
        self.x = x
        self.y = y
        self.obstacle = obstacle

    def equals(self, other):
        return self.x == other.x and self.y == other.y and self.obstacle == other.obstacle


class Obstacle:
    def __init__(self, vertices):
        self.vertices = [Point(x, y) for x, y in vertices]


class Environment:
    def __init__(self, filename):
        self.maxc_x = 200  # number of columns in x
        self.maxr_y = 300  # number of rows in y
        self.boundary = []
        self.obstacles = []
        self.LineList,self.raw_data=self.load_environment(filename)
        # print(self.LineList)
        self.source = None
        self.dest = None
        self.point_x_cord = []
        self.point_y_cord = []
        self.new_obstacles = None
        self.open_line_segments = None
        self.new_sorted_vertices = None
        self.sorted_vertices = None

        self.env = self.parse_input_line(self.raw_data[0])
        self.hole1 = self.parse_input_line(self.raw_data[1])
        self.hole2 = self.parse_input_line(self.raw_data[2])
        self.hole3 = self.parse_input_line(self.raw_data[3])



    def load_environment(self, filename):
        with open(filename, "r") as file_handler:
            raw_data = file_handler.read()
            raw_data = raw_data.split("\n")
            Linelist = self.parse_input_line(raw_data[0])
            temp = self.parse_input_line(raw_data[0])
            self.boundary = [Point(x, y) for x, y in temp]
            temp = self.parse_input_line(raw_data[len(raw_data) - 1])
            self.source = Point(temp[0][0], temp[0][1])
            self.dest = Point(temp[1][0], temp[1][1])
            for i in raw_data[1:len(raw_data) - 1]:
                self.obstacles.append(Obstacle(self.parse_input_line(i)))
            self.parse_obstacles()
        return Linelist,raw_data



    def parse_obstacles(self):
        self.sorted_vertices = []
        for index, obs in enumerate(self.obstacles):
            for j in obs.vertices:
                j.obstacle = index
                self.sorted_vertices.append(j)

        self.sorted_vertices.sort(key=lambda x: x.x)
        self.draw_problem()

        self.new_sorted_vertices = []
        for i in self.sorted_vertices:
            self.new_sorted_vertices.append(Point(i.x, i.y, i.obstacle))

        self.new_obstacles = []
        for index, obs in enumerate(self.obstacles):
            temp_obs = []
            for j in obs.vertices:
                temp_obs.append(Point(j.x, j.y, index))
            self.new_obstacles.append(temp_obs)

    def parse_input_line(self, line):
        temp2 = []
        line = [i.strip() for i in line.split(",")]
        vertex = []
        for index, i in enumerate(line):
            if (i[0] == "("):
                i = i[1:]
            if (i[len(i) - 1] == ")"):
                i = i[:-1]
            vertex.append(int(i))
            if (index % 2 != 0):
                temp2.append(vertex)
                vertex = []
        return temp2

    def segment_intersection(self, a, b, c, d):
        if (intersect(a, b, c, d) == True):
            return line_intersection([a, b], [c, d])
        else:
            return -1


    def find_vertical_lines(self):
        self.open_line_segments = []
        self.cells=[]
        self.to_draw = []

        y_limit_lower = self.boundary[0].y
        y_limit_upper = self.boundary[2].y

        for pt in self.new_sorted_vertices:
            curr_line_segment = [Point(pt.x, y_limit_lower), Point(pt.x, y_limit_upper)]
            lower_obs_pt = curr_line_segment[0]
            upper_obs_pt = curr_line_segment[1]
            upper_gone = False
            lower_gone = False
            break_now = False

            for index, obs in enumerate(self.new_obstacles):
                obs.append(obs[0])
                for vertex_index in range(len(obs) - 1):
                    res = self.segment_intersection(curr_line_segment[0], curr_line_segment[1], obs[vertex_index], obs[vertex_index + 1])
                    if res != -1:
                        if index == pt.obstacle:
                            if not pt.equals(res):
                                if res.y > pt.y:
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
            # if (lower_gone is False):
            #     plt.plot([lower_obs_pt.x, pt.x], [lower_obs_pt.y, pt.y])
            #
            # if (upper_gone is False):
            #     plt.plot([pt.x, upper_obs_pt.x], [pt.y, upper_obs_pt.y])

            # Add to the global segment list
            if (lower_gone and upper_gone):
                self.open_line_segments.append([None, None])
            elif (lower_gone):
                self.open_line_segments.append([None, upper_obs_pt])
            elif (upper_gone):
                self.open_line_segments.append([lower_obs_pt, None])
            else:
                self.open_line_segments.append([lower_obs_pt, upper_obs_pt])

        for index1 in range(len(self.open_line_segments)):
            curr_segment = self.open_line_segments[index1]
            curr_vertex = self.new_sorted_vertices[index1]
            break_now = False
            done = [False, False, True]
            if (curr_segment[0] is None):
                done[0] = True
            if (curr_segment[1] is None):
                done[1] = True
            if (curr_segment[1] is None and self.open_line_segments[index1][0] is None):
                done[2] = False

            for index2 in range(index1 + 1, len(self.open_line_segments)):
                next_segment = self.open_line_segments[index2]
                next_vertex = self.new_sorted_vertices[index2]

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
                            [centroidver([curr_segment[0], curr_vertex]), centroidver([next_segment[0], next_vertex]),
                             0])
                        lines_to_check.append(
                            [centroidver([curr_segment[0], curr_vertex]), centroidver([next_segment[1], next_vertex]),
                             0])
                        trapezoids.append([curr_segment[0], next_segment[0], next_vertex, curr_vertex])
                        trapezoids.append([curr_segment[0], next_vertex, next_segment[1], curr_vertex])
                    elif (next_segment[0] is not None):
                        lines_to_check.append(
                            [centroidver([curr_segment[0], curr_vertex]), centroidver([next_segment[0], next_vertex]),
                             0])
                        trapezoids.append([curr_segment[0], next_segment[0], next_vertex, curr_vertex])
                    elif (next_segment[1] is not None):
                        lines_to_check.append(
                            [centroidver([curr_segment[0], curr_vertex]), centroidver([next_segment[1], next_vertex]),
                             0])
                        trapezoids.append([curr_segment[0], next_vertex, next_segment[1], curr_vertex])
                    else:
                        lines_to_check.append([centroidver([curr_segment[0], curr_vertex]), next_vertex, 0])
                        trapezoids.append([curr_segment[0], next_vertex, curr_vertex])

                if (done[1] is False):
                    if (double_check):
                        double_index2 = len(lines_to_check)
                        lines_to_check.append(
                            [centroidver([curr_segment[1], curr_vertex]), centroidver([next_segment[0], next_vertex]),
                             1])
                        lines_to_check.append(
                            [centroidver([curr_segment[1], curr_vertex]), centroidver([next_segment[1], next_vertex]),
                             1])
                        trapezoids.append([curr_vertex, next_segment[0], next_vertex,
                                           point(curr_segment[1].x, curr_segment[1].y, curr_segment[1].obstacle, 34)])
                        trapezoids.append([curr_vertex, next_vertex, next_segment[1], curr_segment[1]])
                    elif (next_segment[1] is not None):
                        lines_to_check.append(
                            [centroidver([curr_segment[1], curr_vertex]), centroidver([next_segment[1], next_vertex]),
                             1])
                        trapezoids.append([curr_vertex, next_vertex, next_segment[1], curr_segment[1]])
                    elif (next_segment[0] is not None):
                        lines_to_check.append(
                            [centroidver([curr_segment[1], curr_vertex]), centroidver([next_segment[0], next_vertex]),
                             1])
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
                    for index3, obs in enumerate(self.new_obstacles):
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
                        self.cells.append(trapezoids[i])

                if (done[0] == True and done[1] == True and done[2] == True):
                    break

        for i in self.cells:
            i.append(i[0])
            self.to_draw.append(i)


        quad_cells = [i for i in self.cells if len(i) > 3]
        tri_cells = [i for i in self.cells if len(i) == 3]
        others = [i for i in self.cells if len(i) < 3]



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

        quads_to_remove = []
        for index1 in range(len(quad_cells)):
            for index2 in range(len(quad_cells)):
                if (index1 != index2 and quad_cells[index1][0].x == quad_cells[index2][0].x and quad_cells[index1][
                    1].x ==
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
        if (self.boundary[0].x != self.new_sorted_vertices[0].x):
            quad_cells.append([self.boundary[0], point(self.new_sorted_vertices[0].x, y_limit_lower),
                               point(self.new_sorted_vertices[0].x, y_limit_upper), self.boundary[3]])
        if (self.boundary[1].x != self.new_sorted_vertices[len(self.new_sorted_vertices) - 1].x):
            quad_cells.append(
                [point(self.new_sorted_vertices[len(self.new_sorted_vertices) - 1].x, y_limit_lower), self.boundary[1], self.boundary[2],
                 point(self.new_sorted_vertices[len(self.new_sorted_vertices) - 1].x, y_limit_upper)])

        for i in quad_cells:
            x = [j.x for j in i]
            y = [j.y for j in i]
            self.plot_polygon_decomposition(x,y)


        m_polygons = {}

        for index, i in enumerate(self.new_obstacles):
            x = [j.x for j in i]
            y = [j.y for j in i]

            m_polygons[index] = Polygon(np.array([x, y]).T)


    def draw_problem(self):
        bnd_x = [i.x for i in self.boundary]
        bnd_x.append(self.boundary[0].x)
        bnd_y = [i.y for i in self.boundary]
        bnd_y.append(self.boundary[0].y)
        poly_x = []
        poly_y = []

        # Draw the boundary
        plt.plot(bnd_x, bnd_y)



    def plot_polygon_decomposition(self,x ,y):
        number_of_sides_of_poly=4
        self.polycoords=[]

        # distance(shortest)
        line1 = LineString([(x[0], y[0]), (x[1], y[1])])
        line2 = LineString([(x[2], y[2]), (x[3], y[3])])
        line3 = LineString([(x[0], y[0]), (x[3], y[3])])
        line4 = LineString([(x[1], y[1]), (x[2], y[2])])

        for side in range(number_of_sides_of_poly):
            self.polycoords.append((x[side], y[side]))

        polygon = Polygon(self.polycoords)

        distance1 = line1.distance(line2)
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


            self.point_x_cord.append(polycentroid.x)

            self.point_y_cord.append(polycentroid.y)

            self.point_x_cord.extend([x_horizontal_left, x_horizontal_right, x_vertical, x_vertical])

            self.point_y_cord.extend([y_horizontal, y_horizontal, y_vertical_up, y_vertical_down])

            # Draw horizontal line through the centriod
            plt.plot([x_horizontal_left, x_horizontal_right], [y_horizontal, y_horizontal], 'x')

            # Draw vertical line through the centriod
            plt.plot([x_vertical, x_vertical], [y_vertical_up, y_vertical_down], 'x')


    def get_env_graph_points(self):
        graph_points = np.array([self.point_x_cord, self.point_y_cord]).T
        return graph_points

    def get_env_new_obstacles(self):
        return self.new_obstacles

    def get_env(self):
        return self.env




