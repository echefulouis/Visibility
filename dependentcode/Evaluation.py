import scipy.io
import random
import networkx as nx
import numpy as np
from .create_env import Environment
from .create_nx_graph import  NetworkxNeighbours
from collections import defaultdict
from .MultiRobotSubgraph import create_adjacent_matrix_of_nx_graph, multirobotsubgraph, nextstate

input_file_name = 'resources/input_file'
new_environment = Environment(input_file_name)
new_environment.parse_obstacles()
new_environment.find_vertical_lines()
Linelist = new_environment.LineList
nx_points = new_environment.get_env_graph_points()
nx_obstacles = new_environment.get_env_new_obstacles()
def path_length(path):
    totalPathDist = 0
    for i in range(1, len(path)):
        pointA = nx_points[path[i - 1], :]
        pointB = nx_points[path[i], :]
        dist = np.linalg.norm(pointA - pointB)
        totalPathDist = totalPathDist + dist
    return totalPathDist

def patrol_path_simulation(initial_robot_states, divs, DMCs, number_iterations, targetStates):
    catchStates = set()
    robot_states = initial_robot_states.reshape(-1, 1)
    robot = np.squeeze(robot_states.copy())
    # print('robot',robot, robot_states)
    for j in range(number_iterations):
        # states = [state1[j],state2[j],state3[j],state4[j]]
        catch = False
        nextstatevalue = []
        for k, state in enumerate(robot):
            for t in targetStates:
                if state == t and (t, k) not in catchStates:
                    # if state == t and (t, k) not in catchStates:
                    # print (k+1, "caught", t)
                    catchStates.add((t, k))

            if len(catchStates) == len(targetStates):
                # robotid = k
                # print("catch", catchStates)
                catch = True
                break
            else:
                nstate = nextstate(divs[k].index(state), DMCs[k])
                # print(divs[k][nstate])
                nextstatevalue.append(divs[k][nstate])

        if catch:
            #count = count + 1
            # print(count)
            break
        else:
            robot = np.array(nextstatevalue)
            for t, k in catchStates:
                robot[k] = t
            var = robot.reshape(-1, 1)
            robot_states = np.hstack((robot_states, var))
    return robot_states, catchStates


def MCEvaluation(divs, DMCs,nx_points, Graph, targetStates, initial_robot_states, number_iterations, number_trials):


    robot_states = initial_robot_states.reshape(-1,1)

    OGAdj = create_adjacent_matrix_of_nx_graph()
    #print (OGAdj[0])
    OG = Graph
    SubgraphAdj = multirobotsubgraph(div1)
    SG = nx.from_numpy_array(SubgraphAdj)
    Samplingstates = list(set(OG.nodes()) - set(initial_robot_states))

    #targetStates = [random.choice(Samplingstates) for _ in range(num_targets)]
    #targetStates=[]
    targetDistance = float('inf')
    plength = 0
    for i, goal in enumerate(targetStates):
        #samplingstates = list(set(divs[i]) - set(robot_states[i]))
        #targetStates.append(random.choice(samplingstates))
        #targetStates = [15]
        #print(targetStates[i], robot_states[i][0], OG.nodes())
        #for initial in initial_robot_states:
        initial = robot_states[i][0]
        #goal = target
        path = nx.shortest_path(OG, source=initial, target=goal)
        plength += path_length(path)
    targetDistance = plength/len(targetStates)
    #targetDistance = min(targetDistance, plength)

        #print ('pathlength', path_length(path), path)
    #targetStates = [random.choice(divs[i]) for i in range(len(divs))]:
        #print(divs[i], robot_states[i])
    #print(len(divs))
    #targetStates = [25]
    #print('targetstate',targetStates)
    #count=1
    #robotid = 0

    #f = open('evaluationresults1.csv', 'w')
    #evaluationresulttitle = "InitialState" +  ", " + "TargetState"+ ", "+"TargetStatePathLength" + ", " + "DetectionRate" + ", "  + "SamplingRate_Time" + ", "+ "Number_MCSteps" + ", " + "Number_Trails"
    #f.write(evaluationresulttitle)
    #f.write("\n")
    count = 0;
    simulationtime=0
    for trail in range(number_trials):
        catchStates = set()
        robot_states = initial_robot_states.reshape(-1, 1)
        robot = np.squeeze(robot_states.copy())
        #print('robot',robot, robot_states)
        for j in range(number_iterations):
            #states = [state1[j],state2[j],state3[j],state4[j]]
            catch= False
            nextstatevalue = []
            for k, state in enumerate(robot):
                for t in targetStates:
                    if state == t and (t,k) not in catchStates:
                    #if state == t and (t, k) not in catchStates:
                        #print (k+1, "caught", t)
                        catchStates.add((t,k))

                if len(catchStates)==len(targetStates):
                    #robotid = k
                    #print("catch", catchStates)
                    catch = True
                    break
                else:
                    nstate = nextstate(divs[k].index(state), DMCs[k])
                    #print(divs[k][nstate])
                    nextstatevalue.append(divs[k][nstate])

            if catch:
                count=count+1
                #print(count)
                break
            else:
                robot = np.array(nextstatevalue)
                for t, k in catchStates:
                    robot[k] = t
                var = robot.reshape(-1, 1)
            robot_states = np.hstack((robot_states, var))

        maxTime=0
        for _,robotid in catchStates:
            selectedPath = robot_states[robotid]
            totalPathDist= 0
            dt=0.5
            for i in range(1, len(selectedPath)):
                pointA = nx_points[selectedPath[i-1],: ]
                pointB = nx_points[selectedPath[i], :]
                dist = np.linalg.norm(pointA-pointB)
                totalPathDist = totalPathDist + dist

            totalTime = totalPathDist/dt
            maxTime = max(maxTime,totalTime)
            simulationtime += maxTime
            #assert targetState == selectedPath[-1]
            #print('max', totalTime,maxTime, robotid+1, len(selectedPath), selectedPath)
    detectionrate = count / number_trials
    avgdetectiontime = simulationtime/count
    #print('count', count, number_iterations, trail, detectionrate)
    return {'DetectionRate': detectionrate, 'AverageDetectionTime': avgdetectiontime, 'TargetDistance':targetDistance}

    #evaluationresult = (str(initial) + ", " + str(goal) + ", " + str(plength) + ", " +
    #                    str(detectionrate)  + ", " + str(maxTime) + ", "+str(number_iterations)  + ", " + str(trail))
    #print(evaluationresult)
    #f.write(evaluationresult)
    #f.write("\n")

if __name__ == "__main__":
    div1 = [0, 2, 14, 10, 12, 26, 29, 25, 27, 68, 71, 70, 72, 74, 67, 65, 66, 23, 18, 17, 15, 16]
    div2 = [34, 32, 30, 31, 33, 69, 22, 20, 21, 24, 19, 7, 8, 5, 9, 6, 84, 81, 80, 82, 4]
    div3 = [39, 35, 36, 38, 37, 51, 53, 50, 54, 52, 59, 56, 55, 57, 64, 79, 86, 85, 87, 89, 47, 49, 45, 48, 46]
    div4 = [1, 83, 3, 11, 13, 28, 73, 41, 44, 40, 43, 42, 78, 88, 77, 75, 76, 63, 62, 60, 61, 58]
    mat = scipy.io.loadmat('../resources/MCR1.mat')
    DMC1 = mat['Ser']
    mat = scipy.io.loadmat('../resources/MCR2.mat')
    DMC2 = mat['Ser']
    mat = scipy.io.loadmat('../resources/MCR3.mat')
    DMC3 = mat['Ser']
    mat = scipy.io.loadmat('../resources/MCR4.mat')
    DMC4 = mat['Ser']

    DMCs = [DMC1, DMC2, DMC3, DMC4]
    divs = [div1, div2, div3, div4]
    initial_robot_states = np.array([14, 5, 89, 40])
    target_states = [15]
    number_iterations = 100
    number_trails = 1000
    # Create the Networx_x Graph
    number_of_neigbours = 9
    env_nx_graph = NetworkxNeighbours(nx_points, nx_obstacles, number_of_neigbours)
    G = env_nx_graph.get_neighbours()
    #eval=MCEvaluation( divs, DMCs, nx_points, G, target_states, initial_robot_states, number_iterations, number_trails)
    #print(eval)
    robot_states, catchStates = patrol_path_simulation(initial_robot_states, divs, DMCs, number_iterations, target_states)
    print(robot_states)