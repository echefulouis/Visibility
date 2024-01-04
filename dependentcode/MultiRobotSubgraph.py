import networkx as nx
import matplotlib.pyplot as plt
from .create_nx_graph import  NetworkxNeighbours
from shapely.geometry import LineString, MultiPolygon, Polygon
from matplotlib.patches import Polygon as MPolygon
from .create_env import Environment
import numpy as np
from pylab import figure
import csv
import scipy.io
import random

input_file_name = 'resources/input_file'
new_environment = Environment(input_file_name)
new_environment.parse_obstacles()
new_environment.find_vertical_lines()
Linelist = new_environment.LineList

# Get parameters for Networkx graph of the environment
nx_points = new_environment.get_env_graph_points()
nx_obstacles = new_environment.get_env_new_obstacles()

# Create the Networx_x Graph
number_of_neigbours = 9
div1 = [0, 2, 14, 10, 12, 26, 29, 25, 27, 68, 71, 70, 72, 74, 67, 65, 66, 23, 18, 17, 15, 16]
div2 = [34, 32, 30, 31, 33, 69, 22, 20, 21, 24, 19, 7, 8, 5, 9, 6, 84, 81, 80, 82, 4]
div3 = [39, 35, 36, 38, 37, 51, 53, 50, 54, 52, 59, 56, 55, 57, 64, 79, 86, 85, 87, 89, 47, 49, 45, 48, 46]
div4 = [1, 83, 3, 11, 13, 28, 73, 41, 44, 40, 43, 42, 78, 88, 77, 75, 76, 63, 62, 60, 61, 58]


def create_adjacent_matrix_of_nx_graph():
    #G = nx.graph()
    env_nx_graph = NetworkxNeighbours(nx_points, nx_obstacles, number_of_neigbours)
    G = env_nx_graph.get_neighbours()
    MG_Edges = G.edges()
    MG_Edges = list(MG_Edges)
    MG_Adj = np.zeros((G.number_of_nodes(), G.number_of_nodes()))
    for i in range(G.number_of_edges()):
            MG_Adj[MG_Edges[i][0]][MG_Edges[i][1]] = 1
            MG_Adj[MG_Edges[i][1]][MG_Edges[i][0]] = 1
            #np.savetxt('adjmatrix.csv', MG_Adj, delimiter=',')
    return MG_Adj

#def create_adjacent_matrix_of_nx_subgraph():

def multirobotsubgraph(division):
    n = len(division)
    SubgraphAdj = np.zeros((n, n))
    OGAdj = create_adjacent_matrix_of_nx_graph()
    #print(OGAdj.shape)
    for i in range(n):
            for j in range(n):
                    if i != j:
                            a = int(division[i])
                            b = int(division[j])
                            if OGAdj[a, b] == 1:
                                    #print(OGAdj[a, b])
                                    SubgraphAdj[i, j] = 1
    return SubgraphAdj

def draw_subgraph(division):
    SubgraphAdj = multirobotsubgraph(division)
    fig2 = figure(figsize=(18, 16))
    ax2 = fig2.add_subplot(111, aspect='equal', xlabel="S", ylabel="t")
    ax2.set_xlim(-2, 300)
    ax2.set_ylim(-2, 200)

    ax2.add_patch(MPolygon(Linelist, closed=True, fill=False, color='black', label='line 1', linewidth=3))
    ax2.add_patch(MPolygon(new_environment.hole1, closed=True, fill=True, color='gray', label='line 1', linewidth=2))
    ax2.add_patch(MPolygon(new_environment.hole2, closed=True, fill=True, color='gray', label='line 1', linewidth=2))
    ax2.add_patch(MPolygon(new_environment.hole3, closed=True, fill=True, color='gray', label='line 1', linewidth=2))
    nx_coords = nx_points

    #ax2.plot(nx_coords[division, 0], nx_coords[division, 1], 'r')
    # new_nx_cords =[]
    # for i in range(n):
    #        for j in range(n):
    #                if SubgraphAdj[i, j] == 1:
    #                        ax2.plot([nx_coords[division[i], 0], nx_coords[division[j], 0]],
    #                                [nx_coords[division[i], 1], nx_coords[division[j], 1]], 'k')

    #plt.show()
    SG = nx.from_numpy_array(SubgraphAdj)
    # print(division)
    # print(SG.nodes())
    index_mapping = {node: i for i, node in enumerate(SG.nodes())}
    new_nx_cords = [list(nx_coords[node]) for i, node in enumerate(division)]
    # print(new_nx_cords)
    # nx.draw(SG, new_nx_cords, with_labels=True)
    # I = nx.relabel_nodes(nx_graph, index_mapping, copy=True)
    # plt.show()
    # Create Coord Matrix and Store in a file
    #get_adj_matrix(SG)
    #get_coord_matrix(new_nx_cords)


def get_adj_matrix(sub_graph):
    MG_Edges = sub_graph.edges()
    MG_Edges = list(MG_Edges)
    MG_Adj = np.zeros((sub_graph.number_of_nodes(), sub_graph.number_of_edges()))

    for i in range(sub_graph.number_of_edges()):
        MG_Adj[MG_Edges[i][0]][i] = 1
        MG_Adj[MG_Edges[i][1]][i] = -1

    #print(MG_Adj)

    #Open a file for writing Adjacency file
    filename=f"SubAdj.txt"
    #print(filename)

    with open(filename, 'w') as file:
        # Write each row of the matrix to the file, followed by a semicolon
        file.write("[")
        for i, row in enumerate(MG_Adj):
            if i == len(MG_Adj) - 1:
                file.write(" ".join(str(x) for x in row) + "];\n")
            else:
                file.write(" ".join(str(x) for x in row) + ";\n")

def nextstate(i, P):
    """
    A function that takes a transition matrix P, a current state i (assumingstarting at 0) and returns the next state j
    """
    r = random.random()
    cumulativerow = [P[i][0]]
    for k in P[i][1:]:  # Iterate through elements of the transition matrix
        cumulativerow.append(cumulativerow[-1] + k)  # Obtain the cumulative distribution
    for j in range(len(cumulativerow)):
        if cumulativerow[j] >= r:  # Find the next state using inverse sampling
            return j
    return j
def get_coord_matrix(coords_list):
    filename = f"SubCordsMatrix.txt"
    with open(filename, 'w') as file:
        # Write each row of the matrix to the file, followed by a semicolon
        file.write("[")
        for i, point in enumerate(coords_list):
            file.write(" " + (str(point[0])))
        file.write(";\n")
        for i, point in enumerate(coords_list):
            file.write(" " + (str(point[1])))
        file.write("]';\n")

def MCSimulation():
    mat = scipy.io.loadmat('MCR1.mat')
    DMC1 = mat['Ser']
    mat = scipy.io.loadmat('MCR2.mat')
    DMC2 = mat['Ser']
    mat = scipy.io.loadmat('MCR3.mat')
    DMC3 = mat['Ser']
    mat = scipy.io.loadmat('MCR4.mat')
    DMC4 = mat['Ser']
    #print(DMC1.shape)
    #print(DMC2.shape)
    #print(DMC3.shape)
    #print(DMC4.shape)
    print(div1)
    print(div2)
    print(div3)
    print(div4)

    initial_robot_states = [14,5,89,40]
    state1=[]
    state2 = []
    state3 = []
    state4 = []
    state1.append(initial_robot_states[0])
    state2.append(initial_robot_states[1])
    state3.append(initial_robot_states[2])
    state4.append(initial_robot_states[3])

    OGAdj = create_adjacent_matrix_of_nx_graph()
    OG = nx.from_numpy_array(OGAdj)
    SubgraphAdj = multirobotsubgraph(div1)
    SG = nx.from_numpy_array(SubgraphAdj)
    Samplingstates = list(set(OG.nodes()) - set(initial_robot_states))
    print(Samplingstates)
    targetStates = [random.choice(Samplingstates)]
    print('state',targetStates[0], state1, state2, state3, state4)
    #count=1
    #robotid = 0
    catchStates = set()
    for j in range(1000):
        states = [state1[j],state2[j],state3[j],state4[j]]
        catch= False
        for k, s in enumerate(states):
            for t in targetStates:
                if s == t:
                    catchStates.add((t,k))

            if len(catchStates)==len(targetStates):
                #robotid = k
                print("catch", j)
                catch = True
                break
        if catch:
            break
        #if state1[j] == state[0] or state2[j] == state[0] or state3[j] == state[0] or state4[j] == state[0]:
        #    print("catch", j)
        #    break
        #print(state1)
        nstate1 = nextstate(div1.index(state1[j]), DMC1)
        #print('nstateind1', nstate1, 'nstate1', div1[nstate1])
        state1.append(div1[nstate1])


        #print(state2)
        nstate2 = nextstate(div2.index(state2[j]), DMC2)
        #print('nstateind2', nstate2, 'nstate2', div2[nstate2])
        state2.append(div2[nstate2])

        #print(state3)
        nstate3 = nextstate(div3.index(state3[j]), DMC3)
        #print('nstateind3', nstate3, 'nstate3', div3[nstate3])
        state3.append(div3[nstate3])

        #print(state4)
        nstate4 = nextstate(div4.index(state4[j]), DMC4)
        #print('nstateind4', nstate4, 'nstate4', div4[nstate4])
        state4.append(div4[nstate4])

        # #if state[0] in div1:
        #     print(div1.index(state[0]))
        #     nstate1 = nextstate(div1.index(state[0]), DMC1)
        #     print('nstateind', nstate1, 'nstate', div1[nstate1])
        # elif state[0] in div2:
        #     print(div2.index(state[0]))
        #     nstate1 = nextstate(div2.index(state[0]), DMC2)
        #     print('nstateind',nstate1, 'nstate', div2[nstate1])
        # elif state[0] in div3:
        #     print(div3.index(state[0]))
        #     nstate1 = nextstate(div3.index(state[0]), DMC3)
        #     print('nstateind', nstate1, 'nstate', div3[nstate1])
        # elif state[0] in div4:
        #     print(div4.index(state[0]))
        #     nstate1 = nextstate(div4.index(state[0]), DMC4)
        #     print('nstateind', nstate1, 'nstate', div4[nstate1])
        #print ('nextstate', nstate1)

        #count += 1
    #print(state1)
    #print(state2)
    #print(state3)
    #print(state4)
    maxTime=0
    allPaths = [state1,state2,state3,state4]
    for _,robotid in catchStates:
        selectedPath = allPaths[robotid]
        totalPathDist= 0
        dt=0.5
        for i in range(1, len(selectedPath)):
            pointA = nx_points[selectedPath[i-1],: ]
            pointB = nx_points[selectedPath[i], :]
            dist = np.linalg.norm(pointA-pointB)
            totalPathDist = totalPathDist + dist

        totalTime = totalPathDist/dt
        maxTime = max(maxTime,totalTime)
        #assert targetState == selectedPath[-1]
        print(totalTime,maxTime, robotid+1, allPaths[robotid])

if __name__ == "__main__":
    multirobotsubgraph(div4)
    draw_subgraph(div4)
    #MCSimulation()
