#!/usr/bin/env python3

import math
from scipy.spatial import Delaunay
import PIL.Image
import random
import unittest

## tri:        The triangulation from scipy.spatial.Delaunay
## region:     A container of point indices to extend
## dist:       An integer distance to extend
## exclusion:  A container of point indices which extension must not include
##---
## adjacents:  A dictionary of extension point indices to a list of index
##             tuples representing a path from region to adjacent
def extend_region(tri, region, dist=1, exclusion=()):
    adjacents = dict((r, [(r,)]) for r in region)
    nindptr, nindices = tri.vertex_neighbor_vertices
    exclusion = set(exclusion);
    
    while(dist > 0):
        exclusion.update(adjacents.keys());
        localadjacents = dict();
        dist -= 1;

        for rIdx in adjacents.keys():
            for aIdx in nindices[nindptr[rIdx]:nindptr[rIdx+1]]:
                if(aIdx not in region and aIdx not in exclusion):
                    localadjacents.setdefault(aIdx, list())
                    for path in adjacents[rIdx]:
                        localadjacents[aIdx].append(path + (aIdx,))

        adjacents = localadjacents
    
    return(adjacents)

def waterize(tri, dist=1):
    downto = dict()
    depth = dict();

    
    
    for face in tri.convex_hull:
        for idx in face:
            downto[idx] = None
            depth[idx] = 0

    while(True):
        possibilities = extend_region(tri, downto.keys(), dist);
        if(not possibilities):
            break

        pkeys = tuple(possibilities.keys())
        pIdx = random.choice(pkeys)
        path = random.choice(possibilities[pIdx])

        for pathIdx in range(1,len(path)):
            downto[path[pathIdx]] = path[pathIdx-1]
            depth[path[pathIdx]] = 0

    for rIdx in depth.keys():
        qIdx = rIdx
        while(True):
            depth[qIdx] += 1;
            qIdx = downto[qIdx]
            if(qIdx is None):
                break
            pass

    namingishard = 0
    while(True):
        level = extend_region(tri, downto.keys(), 1);
        if(not level):
            break
        for lIdx, paths in level.items():
            path = random.choice(paths)
            downto[lIdx] = path[0]
            depth[lIdx] = namingishard
        namingishard -= 1
    
    return(downto, depth)

def draw_tri(view, tri, downto, depth, undeepest, scale):
    for simpl in tri.simplices:
        for aIdx, bIdx in ((simpl[0], simpl[1]),
                           (simpl[0], simpl[2]),
                           (simpl[1], simpl[2])):                
            if(aIdx not in depth or bIdx not in depth):
                color = (192, 0, 0)
            elif((depth[aIdx] <=0 and depth[bIdx] <=0 and abs(depth[aIdx] - depth[bIdx]) > 1)):
                color = (255, 0, 255)
            elif((depth[aIdx] is None and depth[bIdx] is None) or
                 (depth[aIdx] is None and depth[bIdx] > 0) or
                 (depth[aIdx] > 0 and depth[bIdx] is None) or
                 (depth[aIdx] > 0 and depth[bIdx] > 0 and
                  (downto[aIdx] == bIdx or downto[bIdx] == aIdx))):
                color = (128, 192, 255)
            elif(depth[aIdx] >= 0 and depth[bIdx] >= 0):
                color = (96, 128, 255)
            elif(depth[aIdx] >= (undeepest // 5) and depth[bIdx] >= (undeepest // 5)):
                color = (192 * 3//4, 255 * 3 //4 , 128 * 3 // 4)
            elif(depth[aIdx] >= (undeepest * 2 // 5) and depth[bIdx] >= (undeepest * 2 // 5)):
                color = (192//2, 255//2, 128//2)
            elif(depth[aIdx] >= (undeepest * 3 // 5) and depth[bIdx] >= (undeepest * 3 // 5)):
                color = (255//2, 192//2, 128//2)
            elif(depth[aIdx] >= (undeepest * 4 // 5) and depth[bIdx] >= (undeepest * 4 // 5)):
                color = (255 * 3//4, 192*3//4, 128*3//4)
            else:
                #(224*3//4, 224*3//4, 128*3//4)
                color = (255, 255, 255)


            xDiff = tri.points[bIdx][0] - tri.points[aIdx][0]
            yDiff = tri.points[bIdx][1] - tri.points[aIdx][1]

            for part in range(101):
                xDot = tri.points[aIdx][0] + xDiff * (part / 100)
                yDot = tri.points[aIdx][1] + yDiff * (part / 100)
                view.putpixel((int(xDot * scale), int(yDot * scale)), color)
                pass
            pass
        pass

    for point in tri.points:
        view.putpixel((int(point[0] * scale), int(point[1] * scale)), (255,255,255))
        pass

def draw_elev(view, tri, elev, scale):
    for x in range(view.width):
        for y in range(view.height):
            
    
def tri_dist(tri, aIdx, bIdx):
    adelt = tri.points[aIdx][0] - tri.points[bIdx][0]
    bdelt = tri.points[aIdx][1] - tri.points[bIdx][1]
    return(math.sqrt(adelt * adelt + bdelt * bdelt))

## Limits on aIdx height based on the edge to bIdx, returned as (hardmin, desired, hardmax)
def limitations(aIdx, bIdx, *, tri, downto, depth, height, gradefn, rimHeight=3, minHeight=-40):
    ## a is on the outer edge
    if(downto[aIdx] is None):  
        return(minHeight, minHeight, minHeight)

    ## a and b are river beds
    if(depth[aIdx] > 0 and depth[bIdx] > 0):
        ## a flows into b
        if(downto[aIdx] == bIdx):
            return(height[bIdx],
                   height[bIdx] + tri_dist(tri, aIdx, bIdx) * gradefn(depth[aIdx]),
                   None)

        ## b flows into a
        if(depth[bIdx] > 0 and downto[bIdx] == aIdx):
            return(None,
                   height[bIdx] - tri_dist(tri, aIdx, bIdx) * gradefn(depth[bIdx]),
                   height[bIdx])

        ## neither flow into the other
        return(None, None, None)
    
    ## a is a shore and a flows into b
    if(depth[aIdx] == 0 and downto[aIdx] == bIdx):
        return(height[bIdx],
               height[bIdx] + tri_dist(tri, aIdx, bIdx) * gradefn(depth[bIdx]),
               None)
    ## b is a shore and b flows into a
    if(depth[bIdx] == 0 and downto[bIdx] == aIdx):
        return(None,
               height[bIdx] - tri_dist(tri, aIdx, bIdx) * gradefn(depth[aIdx]),
               height[bIdx])

    ## a is a rim and a flows into b
    if(depth[aIdx] == -1 and downto[aIdx] == bIdx and downto[bIdx] is not None):
        return(height[bIdx],
               max(height[bIdx] + tri_dist(tri, aIdx, bIdx) * gradefn(depth[aIdx]),
                   height[downto[bIdx]] + rimHeight),
               None)

    ## b is a rim and b flows into a
    if(depth[bIdx] == -1 and downto[bIdx] == aIdx and downto[aIdx] is not None):
        return(None,
               None,  ## Note:  This is independent because be can be set independently of a
               height[bIdx])

    ## a is above b
    if(depth[aIdx] < depth[bIdx]):
        return(height[bIdx],
               height[bIdx] + tri_dist(tri, aIdx, bIdx) * gradefn(depth[aIdx]),
               None)

    ## b is above a
    if(depth[aIdx] > depth[bIdx]):
        return(None,
               height[bIdx] - tri_dist(tri, aIdx, bIdx) * gradefn(depth[bIdx]),
               height[bIdx])

    ## Otherwise, no restrictions
    return(None, None, None)

def oneraise(tri, downto, depth, height, gradefn):
    allIdx = list(height.keys())
    random.shuffle(allIdx);
    newheight = dict()
    nindptr, nindices = tri.vertex_neighbor_vertices
    violations = 0

    for aIdx in allIdx:
        minH = None
        desired = list()
        maxH = None
        for bIdx in nindices[nindptr[aIdx]:nindptr[aIdx+1]]:
            minL, typL, maxL = limitations(aIdx, bIdx,
                                           tri=tri, downto=downto, depth=depth,
                                           height=height, gradefn=gradefn)
            if(minL is not None and
               (minH is None or minL > minH)):
                minH = minL
            if(typL is not None):
                desired.append(typL)
            if(maxL is not None and
               (maxH is None or maxL < maxH)):
                maxH = maxL
            pass

        if((not desired) or
           (minH is not None and maxH is not None and minH > maxH)):
            newheight[aIdx] = (minH + maxH / 2)
            violations += 1
        else:
            hooboy = sum(desired) / len(desired)
            if(minH is not None and minH > hooboy):
                hooboy = minH
            if(maxH is not None and maxH < hooboy):
                hooboy = maxH

            newheight[aIdx] = hooboy
        pass

    return(newheight, violations)
            
            
                

class _IlMap_ut(unittest.TestCase):
    def test_randscatter(self):
        width = 36000
        height = 36000
        scale = 1/20
        dist = 16.666 * 7
        separate = 1500 / dist  ## about 1.5km apart
        distsq = dist * dist

        points = list()
        for idx in range(int(width * height / distsq) + 3):
            newpoint = (random.uniform(0, width), random.uniform(0, height));
            if(1):
            #if((newpoint[0] - (width/2))**2 + (newpoint[1]-(height/2))**2 <= width * height / 4):
            ## for oldpoint in points:
            ##     if((newpoint[0]-oldpoint[0])**2 + (newpoint[1]-oldpoint[1])**2 < distsq):
            ##         break
            ##     pass
            ## else:
                points.append(newpoint)
                pass
            pass
        
        print(f'POINTS:  {len(points)}')
        
        tri = Delaunay(points)

        downto, depth = waterize(tri, separate)

        ordered = sorted(depth.values())
        leastdepth = ordered[0];
        ordered = [v for v in ordered if v > 0]
        depthbits = tuple(ordered[i * len(ordered) // 5] for i in range(5))

        def gradefn(d):
            return(0.01 if d >= depthbits[4] else
                   0.02 if d >= depthbits[3] else
                   0.05 if d >= depthbits[2] else
                   0.10 if d >= depthbits[1] else
                   0.20 if d >= depthbits[0] else
                   0.01 if d >= leastdepth * 1 // 5 else
                   0.05 if d >= leastdepth * 2 // 5 else
                   0.15 if d >= leastdepth * 3 // 5 else
                   0.30 if d >= leastdepth * 4 // 5 else
                   0.75)
        
        ## histo = dict()        
        ## for value in depth.values():
        ##     histo.setdefault(value, 0)
        ##     histo[value] += 1
        ## 
        ## histosort = sorted(histo.items());
        ## if(1):
        ##     for key, value in sorted(histo.items()):
        ##         print(f"DEPTH {key:2d}: {value}")

        print(f'Depth:  [{ordered[-1]}:{depthbits}:{leastdepth}]')
        
        view = PIL.Image.new('RGB', (math.ceil(width * scale),
                                     math.ceil(height * scale)))
        
        draw_tri(view, tri, downto, depth, leastdepth, scale)
        view.show()

        height = dict( (idx, 0) for idx in depth.keys())

        for idx in range(1000):
            height, violations = oneraise(tri, downto, depth, height, gradefn)
            allheights = sorted(height.values())
            print(f"HEIGHT:  [{allheights[0]}:{allheights[-1]}] with {violations} violations")
        
        pass
    
            
            
