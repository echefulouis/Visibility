import math
import numpy as np
import matplotlib.pyplot as plt
from dijkstra_search import DijkstraSearch
from geometry import Geometry


class VisibilityRoadMap:

    def __init__(self, expand_distance, do_plot=False, show_animation=True):
        self.expand_distance = expand_distance
        self.do_plot = do_plot
        self.show_animation = show_animation

    def planning(self, start_x, start_y, goal_x, goal_y, obstacles):

        nodes = self.generate_visibility_nodes(start_x, start_y,
                                               goal_x, goal_y, obstacles)

        road_map_info = self.generate_road_map_info(nodes, obstacles)

        if self.do_plot:
            self.plot_road_map(nodes, road_map_info)
            plt.pause(1.0)

        rx, ry = DijkstraSearch(self.show_animation).search(
            start_x, start_y,
            goal_x, goal_y,
            [node.x for node in nodes],
            [node.y for node in nodes],
            road_map_info
        )

        return rx, ry

    def generate_visibility_nodes(self, start_x, start_y, goal_x, goal_y,
                                  obstacles):

        # add start and goal as nodes
        nodes = [DijkstraSearch.Node(start_x, start_y),
                 DijkstraSearch.Node(goal_x, goal_y, 0, None)]

        # add vertexes in configuration space as nodes
        for obstacle in obstacles:

            cvx_list, cvy_list = self.calc_vertexes_in_configuration_space(
                obstacle.x_list, obstacle.y_list)

            for (vx, vy) in zip(cvx_list, cvy_list):
                nodes.append(DijkstraSearch.Node(vx, vy))

        if self.do_plot:
            for node in nodes:
                plt.plot(node.x, node.y, "xr")

        return nodes

    def calc_vertexes_in_configuration_space(self, x_list, y_list):
        x_list = x_list[0:-1]
        y_list = y_list[0:-1]
        cvx_list, cvy_list = [], []

        n_data = len(x_list)

        for index in range(n_data):
            offset_x, offset_y = self.calc_offset_xy(
                x_list[index - 1], y_list[index - 1],
                x_list[index], y_list[index],
                x_list[(index + 1) % n_data], y_list[(index + 1) % n_data],
            )
            cvx_list.append(offset_x)
            cvy_list.append(offset_y)

        return cvx_list, cvy_list

    def generate_road_map_info(self, nodes, obstacles):

        road_map_info_list = []

        for target_node in nodes:
            road_map_info = []
            for node_id, node in enumerate(nodes):
                if np.hypot(target_node.x - node.x,
                            target_node.y - node.y) <= 0.1:
                    continue

                is_valid = True
                for obstacle in obstacles:
                    if not self.is_edge_valid(target_node, node, obstacle):
                        is_valid = False
                        break
                if is_valid:
                    road_map_info.append(node_id)

            road_map_info_list.append(road_map_info)

        return road_map_info_list

    @staticmethod
    def is_edge_valid(target_node, node, obstacle):

        for i in range(len(obstacle.x_list) - 1):
            p1 = Geometry.Point(target_node.x, target_node.y)
            p2 = Geometry.Point(node.x, node.y)
            p3 = Geometry.Point(obstacle.x_list[i], obstacle.y_list[i])
            p4 = Geometry.Point(obstacle.x_list[i + 1], obstacle.y_list[i + 1])

            if Geometry.is_seg_intersect(p1, p2, p3, p4):
                return False

        return True

    def calc_offset_xy(self, px, py, x, y, nx, ny):
        p_vec = math.atan2(y - py, x - px)
        n_vec = math.atan2(ny - y, nx - x)
        offset_vec = math.atan2(math.sin(p_vec) + math.sin(n_vec),
                                math.cos(p_vec) + math.cos(
                                    n_vec)) + math.pi / 2.0
        offset_x = x + self.expand_distance * math.cos(offset_vec)
        offset_y = y + self.expand_distance * math.sin(offset_vec)
        return offset_x, offset_y

    @staticmethod
    def plot_road_map(nodes, road_map_info_list):
        for i, node in enumerate(nodes):
            for index in road_map_info_list[i]:
                plt.plot([node.x, nodes[index].x],
                         [node.y, nodes[index].y], "-b")
