from pylab import figure, show, rand
from matplotlib.patches import Ellipse,Circle,Rectangle
from matplotlib.patches import Polygon as MPolygon
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


def writedat(filename, x, y, xprecision=1, yprecision=1):
    with open(filename,'w') as f:
        for a, b in itertools.zip(x, y):
            print >> f, "%.*g\t%.*g" % (xprecision, a, yprecision, b)
            

def genVisibMatrix(guardno):
    a = str(os.popen("./main gridsmallenvironment gridsmallguards "+str(guardno)).read())
    l = a.split() 
    mat=[]
    print (a)
    for i in range(0,len(l)-1,2):
        mat.append ([float(l[i]),float(l[i+1])])
    return mat

def pointOnPolygon(point, poly):
    for i in range(len(poly)):
        a, b = poly[i - 1], poly[i]
        if abs(dist(a, point) + dist(b, point) - dist(a, b)) < EPSILON:
            return true
    return false

def point_in_polyen(x,y,poly):

   # check if point is a vertex
   if (x,y) in poly: return True

   # check if point is on a boundary
   for i in range(len(poly)):
      p1 = None
      p2 = None
      if i==0:
         p1 = poly[0]
         p2 = poly[1]
      else:
         p1 = poly[i-1]
         p2 = poly[i]
      if p1[1] == p2[1] and p1[1] == y and x > min(p1[0], p2[0]) and x < max(p1[0], p2[0]):
         return True

   n = len(poly)
   inside = False
   print(poly[0])

   for i in range(n+1):
      p2x,p2y = poly[i % n]
      if y > min(p1y,p2y):
         if y <= max(p1y,p2y):
            if x <= max(p1x,p2x):
               if p1y != p2y:
                  xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
               if p1x == p2x or x <= xints:
                  inside = not inside
      p1x,p1y = p2x,p2y

   if inside: return inside
   else: return inside
def point_in_poly(x,y,poly):

    n = len(poly)
    inside = False

    p1x,p1y = poly[0]
    for i in range(n+1):
        p2x,p2y = poly[i % n]
        if y > min(p1y,p2y):
            if y <= max(p1y,p2y):
                if x <= max(p1x,p2x):
                    if p1y != p2y:
                        xints = (y-p1y)*(p2x-p1x)/(p2y-p1y)+p1x
                    if p1x == p2x or x <= xints:
                        inside = not inside
        p1x,p1y = p2x,p2y

    return inside
#https://stackoverflow.com/questions/54284984/sectors-representing-and-intersections-in-shapely
def sector(center, start_angle, end_angle, radius, steps=200):
    def polar_point(origin_point, angle,  distance):
        return [origin_point.x + math.sin(math.radians(angle)) * distance, origin_point.y + math.cos(math.radians(angle)) * distance]

    if start_angle > end_angle:
        start_angle = start_angle - 360
    else:
        pass
    step_angle_width = (end_angle-start_angle) / steps
    sector_width = (end_angle-start_angle)
    segment_vertices = []

    segment_vertices.append(polar_point(center, 0,0))
    segment_vertices.append(polar_point(center, start_angle,radius))

    for z in range(1, steps):
        segment_vertices.append((polar_point(center, start_angle + z * step_angle_width,radius)))
    segment_vertices.append(polar_point(center, start_angle+sector_width,radius))
    segment_vertices.append(polar_point(center, 0,0))
    return Polygon(segment_vertices)

def insidenvironment(x,y):
    if point_in_poly(x,y, env):
        return True
    else:
        return False

def insideobstacle(x,y):
    if point_in_poly(x,y, obs1):
        return True
    else:
        return False

def cellcenter(x,y):
    x_center=(floor(x) + floor(x)+1)*0.5
    y_center= (floor(y)+floor(y)+1)*0.5
    return  x_center, y_center

def cellposition(cell):
    cc=cell%(maxc_x*maxr_y)
    x=cc%maxc_x
    y=floor(cc/maxc_x)
    return (x+(x+1))*0.5, (y+(y+1))*0.5 
#mpl.rcParams.update({'font.size': 18})

def midcellnumber(y,x):
    cell = floor(y)*maxc_x+floor(x)
        
    return int(cell)
def cellnumber(y,x):
    cell = y*maxc_x+x
        
    return cell

def main():
    global env, obs1, obs2, maxc_x, maxr_y
    LineList=[ [0, 0], [8, 0], [8,8],[0, 8]]
    env=[(0,0),(8,0),(8,8),(0,8)]
    hole1=[ [2,5],[5,5],[5,6], [6,6], [6,3], [3,3],[3, 2],  [2, 2]  ]
    #hole2=[[18,11],[26,11],[26,5],[18,5]]
     
    obs1=[(2,5),(5,5),(5,6),(6,6),(6, 3),(3,3),(3,2), (2,2)]
    #obs2=[(18,5),(26,5),(26,11),(18,11)]
    maxc_x =8 #number of colmuns in x
    maxr_y =8 #number of rows in y
    
    freecell=[]
    ltable = {}
    lvispoly={}
    #maxc_x*maxr_y
    for i in range (0,maxc_x*maxr_y):
        x,y=cellposition(i)
        print(x,y,cellposition(i))
        if not insidenvironment(x, y) or insideobstacle(x, y) :
            continue
        else:
            freecell.append(i)
            np.savetxt('gridsmallguards',(x,y),fmt='%5.1f')
            #polypointlist=[(x-1,y-2),(x+1,y-2),(x+2,y-1),(x+2,y+1),(x+1,y+2),(x-1,y+2),(x-2,y+1),(x-2,y-1)]
            trianglepointlist=[(x+1,y-2),(x+2,y-1),(x+2,y+1),(x+1,y+2), (x,y)]
            unitA = Circle((x, y), .2, facecolor='none',fill=True,edgecolor=(0,0.8,0.8), linewidth=2, alpha=0.5)
            ccellno=int(midcellnumber(y,x))
            #trianglepointlist=[(x+4,y-1),(x+4,y+1), (x,y)]
            #circlepoly=Polygon(polypointlist)
            #trianglepoly = Polygon(trianglepointlist)
            circlepoly= Point(x, y).buffer(2)
            circleview = Polygon(env).intersection(circlepoly)
            sectorpoly = sector(Point(x,y), 315, 45, 2)
            #circle=list(trianglepoly.exterior.coords)
            circlecoords= list(circleview.exterior.coords)
            polygon1=genVisibMatrix(0)
            vispolygon = genVisibMatrix(0)

            points=[]
            for i in range(len(polygon1)):
                points.append(Point(polygon1[i][0],polygon1[i][1]) ) 
            polygon1=Polygon([[p.x, p.y] for p in points])
 
            #intersectingpoly=polygon1.intersection(circlepoly)
            #intersectingpoly=polygon1.intersection(trianglepoly)
            intersectingpoly = polygon1.intersection(circlepoly)
            circleview = polygon1.intersection(circlepoly)
            sectorview = polygon1.intersection(sectorpoly)
            circlecoords = list(circleview.exterior.coords)
            sectorcoords = list(sectorview.exterior.coords)
            intspoints = np.array(intersectingpoly)
            print ('Intersect', intspoints)
            ix, iy = intersectingpoly.exterior.xy

            intspoly=[]
            for i in range(len(ix)):
                intspoly.append([int(ix[i]),int(iy[i])])
            lvispoly.update({ccellno:intspoly})
 
            fig = figure( figsize=(18, 16))
            ax = fig.add_subplot(111, aspect='equal',xlabel="S",ylabel="t")
    
   
            if (len(intspoly)>1):
                #ax.add_patch(MPolygon(intspoly, closed=True, fill=True, color='g', linewidth=0))
                ax.add_patch(MPolygon(vispolygon, closed=True, fill=True, color='g', linewidth=0))
                ax.add_patch(MPolygon(sectorcoords, closed=True, fill=True, color='r', linewidth=0))

            ax.add_patch(MPolygon(LineList, closed=True, fill=False, color='black',label='line 1', linewidth=3))
            ax.add_patch(MPolygon(hole1, closed=True, fill=True, color='black',label='line 1', linewidth=2))
            ax.add_patch(unitA)
            ax.set_xlim(-2, 10)
            ax.set_ylim(-2, 10)
    
     	    #timestr = time.strftime("%Y%m%d-%H%M%S")   
            #pylab.savefig('output/'+timestr+'No_'+str(l)+".jpeg", bbox_inches='tight')
            show()
         
             
            maxpolyx=-10
            maxpolyy=-10
            minpolyx=1000
            minpolyy=1000
  
            for j in range(0,len(intspoly)):
                if(intspoly[j][0]>maxpolyx):
                    maxpolyx=intspoly[j][0]
                if(intspoly[j][1]>maxpolyy):
                    maxpolyy=intspoly[j][1]
                if(intspoly[j][0]<minpolyx):
                    minpolyx=intspoly[j][0]
                if(intspoly[j][1]<minpolyy):
                    minpolyy=intspoly[j][1]    
            #print floor(maxpolyx), ceil(int(minpolyx)),floor(maxpolyy),ceil(int(minpolyy))
            maxpolyx=floor(maxpolyx)
            minpolyx=ceil(int(minpolyx))
            maxpolyy=floor(maxpolyy)
            minpolyy=ceil(int(minpolyy))
            
            
            #instpoly=Polygon(intspoly)
            #print instpoly
            ltable.setdefault(ccellno, [])
            for j in range(int(minpolyx), int(maxpolyx)+1):
                for k in range(int(minpolyy), int(maxpolyy)+1):
                    #print j, k
                    #print intersectingpoly.intersects(Point(1,1))
                    #print point_in_polyen(0,0, intspoly), point_in_polyen(1,0, intspoly), point_in_polyen(1,1, intspoly), pointOnPolygon([0,1], intspoly)
                    if intersectingpoly.intersects(Point(j,k)) and intersectingpoly.intersects(Point(j+1,k)) and intersectingpoly.intersects(Point(j+1,k+1)) and intersectingpoly.intersects(Point(j,k+1)):
                        #print cellnumber(j, k), j, k
                        viscell=midcellnumber(k, j)
                        #if viscell==204:
                        #    print viscell, cellposition(viscell),j,k
                        ltable[ccellno].append(viscell) 
                        


if __name__ == "__main__":main()











    

