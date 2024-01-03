import networkx as nx
import matplotlib.pyplot as plt
from sklearn.neighbors import NearestNeighbors
from sklearn.metrics.pairwise import pairwise_distances
from shapely.geometry import LineString, MultiPolygon, Polygon
from matplotlib.patches import Polygon as MPolygon
from .create_env import Environment
from pylab import figure
import numpy as np
import csv


class NetworkxNeighbours:
    def __init__(self, networkx_points, new_obstacles, num_neigbours):
        self.networkx_points = networkx_points
        self.new_obstacles = new_obstacles
        self.num_neigbours = num_neigbours
        self.graph = nx.Graph()

    def fit_model(self):
        neigh = NearestNeighbors(n_neighbors=self.num_neigbours)
        # Fit the model
        neigh.fit(self.networkx_points)

        # Calculate pairwise distances between points
        distances = pairwise_distances(self.networkx_points, self.networkx_points, metric="euclidean")
        # Create weighted adjacency matrix from distances
        adj_matrix = distances

        # Predict the nearest neighbors
        neighbors = neigh.kneighbors_graph(self.networkx_points, mode='distance')
        self.graph = nx.from_numpy_array(neighbors.toarray())
        # Remove self loop
        self.graph.remove_edges_from(nx.selfloop_edges(self.graph))

        # Add weights to edges
        for i, j in self.graph.edges():
            self.graph[i][j]['weight'] = adj_matrix[i, j]

    def get_adj_matrix(self):
        return pairwise_distances(self.networkx_points, self.networkx_points, metric="euclidean")

    def plot_obstacles(self):
        ax = plt.gca()
        for index, i in enumerate(self.new_obstacles):
            g = [j.x for j in i]
            h = [j.y for j in i]

            polygon_patch = MPolygon(np.array([g, h]).T, closed=True, color="black")
            ax.add_patch(polygon_patch)
        # Dictionary for multipolygon
        polygons = {}

        for index, i in enumerate(self.new_obstacles):
            x = [j.x for j in i]
            y = [j.y for j in i]

            polygons[index] = Polygon(np.array([x, y]).T)

        # loop through the dictionary and create a multipolygon
        mp = MultiPolygon([polygons[i] for i in polygons])

        return mp

    def check_collision(self, point1, point2, mp):
        lines = LineString([point1, point2])
        return lines.intersects(mp)

    def remove_collisions(self):
        mp = self.plot_obstacles()
        line_nodes = self.graph.edges()

        # Check for collision between each line and collision object
        for node in line_nodes:
            point1 = (list(self.networkx_points[node[0]])[0], list(self.networkx_points[node[0]])[1])
            point2 = (list(self.networkx_points[node[1]])[0], list(self.networkx_points[node[1]])[1])
            collision = self.check_collision(point1, point2, mp)

            # Check Colllision with Multi-polygon
            if collision:
                self.graph.remove_edge(node[0], node[1])



    def draw_graph(self):
        nx.draw(self.graph, self.networkx_points, with_labels=True, node_size=20)
        #plt.show()

    def create_adjacent_matrix_of_nx_graph(self):

        #adj = [[0] * 90 for _ in range(90)]
        #np.savetxt('adjmatrix.csv', adj, delimiter=',')
        G = self.graph
        MG_Edges = G.edges()
        MG_Edges = list(MG_Edges)

        MG_Adj = np.zeros((G.number_of_nodes(), G.number_of_nodes()))
        for i in range(G.number_of_edges()):
            MG_Adj[MG_Edges[i][0]][MG_Edges[i][1]] = 1
            MG_Adj[MG_Edges[i][1]][MG_Edges[i][0]] = 1
        np.savetxt('adjmatrix.csv', MG_Adj, delimiter=',')

        #div1 = [0, 87, 85, 86, 79, 89, 47, 49, 45, 48, 46, 34, 31, 30, 32, 33, 69, 22, 20, 21, 23, 18]
        #Adj1 = np.zeros((len(div1), len(div1)))
        #for i in range(G.number_of_edges()):
        #    MG_Adj[MG_Edges[i][0]][MG_Edges[i][1]] = 1
        #    MG_Adj[MG_Edges[i][1]][MG_Edges[i][0]] = 1
        #np.savetxt('adjmatrix.csv', MG_Adj, delimiter=',')


        #div2 = [72, 70, 71, 74, 66, 65, 67, 36, 39, 35, 38, 37, 51, 53, 50, 54, 52, 59, 56, 55, 57, 58]
        #div3 = [1, 83, 3, 11, 10, 13, 16, 15, 17, 24, 19, 7, 8, 5, 9, 6, 84, 81, 80, 82, 4]
        #div4 = [64, 62, 60, 61, 63, 76, 75, 77, 88, 78, 42, 43, 40, 44, 41, 73, 68, 27, 28, 25, 29, 26, 12, 14, 2]
        return MG_Adj
    def divisionsubgraph(self, division):
        self.networkx_points
        SG = nx.Graph()
        n = len(division)
        SubgraphAdj = np.zeros((n,n))
        OG = self.graph

        OGAdj= self.create_adjacent_matrix_of_nx_graph()
        print(OGAdj.shape)
        for i in range(n):
            for j in range (n):
                if i !=j:
                    a = int(division[i])
                    b = int(division[j])
                    if OGAdj[a, b] == 1:
                        print(OGAdj[a, b])
                        SubgraphAdj[i,j]= 1
        # Create Environment
        input_file_name = 'input_file'
        new_environment = Environment(input_file_name)
        new_environment.parse_obstacles()
        new_environment.find_vertical_lines()
        Linelist = new_environment.LineList

        fig = figure(figsize=(18, 16))
        ax = fig.add_subplot(111, aspect='equal', xlabel="S", ylabel="t")
        ax.set_xlim(-2, 300)
        ax.set_ylim(-2, 200)

        ax.add_patch(MPolygon(Linelist, closed=True, fill=False, color='black', label='line 1', linewidth=3))
        ax.add_patch(MPolygon(new_environment.hole1, closed=True, fill=True, color='gray', label='line 1', linewidth=2))
        ax.add_patch(MPolygon(new_environment.hole2, closed=True, fill=True, color='gray', label='line 1', linewidth=2))
        ax.add_patch(MPolygon(new_environment.hole3, closed=True, fill=True, color='gray', label='line 1', linewidth=2))

        nx_coords = self.networkx_points
        #ax.plot(nx_coords[division,0], nx_coords[division,1], 'r')
        for i in range(n):
            for j in range(n):
                if SubgraphAdj[i,j] == 1:
                    ax.plot([nx_coords[division[i],0],nx_coords[division[j],0]],[nx_coords[division[i],1],nx_coords[division[j],1]],'k')
        #plt.show()

        # neighbours = []
        # for node in OG.nodes():
        #     #if node in division:
        #     #    neighbours.append(node)
        #     if node in division:
        #         neighbours.append(node)
        #         SG.add_node(node)
        #         neighbor_list = [n for n in OG.neighbors(node)]
        #         for neigh in neighbor_list:
        #             if neigh in division:
        #                 neighbours.append(neigh)
        #                 SG.add_edge(node, neigh)
        #                 #print (node, "--> Neighlist", neighbor_list, 'neigh', neigh,'div', division )
        #
        # # Make the list of neigbours unique
        # neighbours = list(set(neighbours))
        #
        # # create a subgraph of the unique list of nodes
        # #nx_subgraph = OG.subgraph(neighbours).copy()
        # print(SG.nodes())
        # print(SG.edges())
        #
        #
        # index_mapping = {node: i for i , node in enumerate(SG.nodes())}
        # new_nx_cords = [list(nx_coords[node]) for i, node in enumerate(SG.nodes())]
        #
        #
        #
        # print (new_nx_cords)
        # print (len(new_nx_cords))
        # I = nx.relabel_nodes(SG,index_mapping,copy=True)
        # print(len(I.nodes()))
        # # Create Coord Matrix and Store in a file
        # #self.get_adj_matrix(I,robot_index)
        # #self.get_coord_matrix(new_nx_cords,robot_index)
        #
        # #new_edge_list = [(index_mapping[u], index_mapping[v]) for u, v in division_edges]
        # #print(list(set(node for edge in new_edge_list for node in edge)))
        #
        # # Create Environment
        # input_file_name = 'input_file'
        # new_environment = Environment(input_file_name)
        # new_environment.parse_obstacles()
        # new_environment.find_vertical_lines()
        # Linelist = new_environment.LineList
        #
        #
        # fig = figure(figsize=(18, 16))
        # ax = fig.add_subplot(111, aspect='equal', xlabel="S", ylabel="t")
        # ax.set_xlim(-2, 300)
        # ax.set_ylim(-2, 200)
        #
        # ax.add_patch(MPolygon(Linelist, closed=True, fill=False, color='black', label='line 1', linewidth=3))
        # ax.add_patch(MPolygon(new_environment.hole1, closed=True, fill=True, color='gray', label='line 1', linewidth=2))
        # ax.add_patch(MPolygon(new_environment.hole2, closed=True, fill=True, color='gray', label='line 1', linewidth=2))
        # ax.add_patch(MPolygon(new_environment.hole3, closed=True, fill=True, color='gray', label='line 1', linewidth=2))
        # print (len(SG.nodes()), len(new_nx_cords))
        # nx.draw(SG, new_nx_cords, with_labels=True)
        # #nx.draw(nx_subgraph, new_nx_cords, edgelist=nx_subgraph.edges(), edge_color='b', width=5.0)
        #

        #return SG

    def create_adjacent_matrix_of_nx_graph_from_csv(self, filename):
        division =[]
        with open(filename, 'r') as file:
            csvreader = csv.reader(file)
            #csvreader = csv.reader("pointsvisibility.")
            for row in csvreader:
                #print(row, row[0], row[1])
                k=0
                for i in range(len(self.networkx_points)):
                    #print(round(self.networkx_points[i][0], 5), row[0], round(self.networkx_points[i][1], 5), row[1])
                    if (str(round(self.networkx_points[i][0], 5)) == str(row[0])) and (str(round(self.networkx_points[i][1], 5)) ==  str(row[1])):
                        #print (i, round(self.networkx_points[i][0], 5), row[0], round(self.networkx_points[i][1], 5), row[1])
                        division.append(i)
                        k=i;

                print(row, row[0], row[1], k)
        print ("Division", division)
        Adj1 = np.zeros((len(division), len(division)))
        OG = self.graph

        neighbours = []
        for node in OG.nodes():
            #if node in division:
            #    neighbours.append(node)
            if node in division:
                neighbours.append(node)
                neighbor_list = [n for n in OG.neighbors(node)]
                for neigh in neighbor_list:
                    if neigh in division:
                        neighbours.append(neigh)
                print (node, "--> Neigh", neighbor_list, "all neigh", neighbours)

        # Make the list of neigbours unique
        neighbours = list(set(neighbours))

        # create a subgraph of the unique list of nodes
        nx_subgraph = OG.subgraph(neighbours).copy()

        div = []
        k = 0
        with open("pvis.csv", 'r') as file:
            csvreader = csv.reader(file)
            for row in csvreader:
                # print(row, row[0], row[1], k)
                # print (i, round(self.networkx_points[i][0], 5), row[0], round(self.networkx_points[i][1], 5), row[1])
                if k in SG.nodes():
                    m = float(row[0])
                    n = float(row[1])
                    div.append([m, n])
                k = k + 1;

                # print(row, row[0], row[1], k)
        print("Division", div)
        #adj = [[0] * 90 for _ in range(90)]
        #np.savetxt('adjmatrix.csv', adj, delimiter=',')
        #G = self.graph
        #MG_Edges = G.edges()
        #MG_Edges = list(MG_Edges)

        #MG_Adj = np.zeros((G.number_of_nodes(), G.number_of_nodes()))
        #for i in range(G.number_of_edges()):
        #    MG_Adj[MG_Edges[i][0]][MG_Edges[i][1]] = 1
        #    MG_Adj[MG_Edges[i][1]][MG_Edges[i][0]] = 1
        #np.savetxt('adjmatrix.csv', MG_Adj, delimiter=',')
        #return MG_Adj
    def get_neighbours(self):
        self.fit_model()
        self.remove_collisions()
        #self.draw_graph()
        #self.create_adjacent_matrix_of_nx_graph()
        #div=[0, 87, 85, 86, 79, 89, 47, 49, 45, 48, 46, 34, 31, 30, 32, 33, 69, 22, 20, 21, 23, 18,]
        #div1 = [0, 2, 14, 10, 12, 26, 29, 25, 27, 68, 71, 70, 72, 74, 67, 65, 66, 23, 18, 17, 15, 16]

        #div2 = [34, 32, 30, 31, 33, 69, 22, 20, 21, 24, 19, 7, 8, 5, 9, 6, 84, 81, 80, 82, 4]
        #div3 = [39, 35, 36, 38, 37, 51, 53, 50, 54, 52, 59, 56, 55, 57, 64, 79, 86, 85, 87, 89, 47, 49, 45, 48, 46]
        #div4 = [1, 83, 3, 11, 13, 28, 73, 41, 44, 40, 43, 42, 78, 88, 77, 75, 76, 63, 62, 60, 61, 58]
        #self.divisionsubgraph(div)
        #self.create_adjacent_matrix_of_nx_graph_from_csv("1.csv")

        return self.graph.copy()
