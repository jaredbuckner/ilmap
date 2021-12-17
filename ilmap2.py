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
            
    return(downto, depth)

class _IlMap_ut(unittest.TestCase):
    def test_randscatter(self):
        width = 2000
        height = 2000
        dist = 9
        distsq = dist * dist

        points = list()
        for idx in range(int(width * height / distsq) + 3):
            newpoint = (random.uniform(0, width), random.uniform(0, height));
            if(1):
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

        downto, depth = waterize(tri, 8)

        histo = dict()
        for value in depth.values():
            histo.setdefault(value, 0)
            histo[value] += 1

        for key, value in sorted(histo.items()):
            print(f"DEPTH {key:2d}: {value}")
        
        
        view = PIL.Image.new('RGB', (math.ceil(width),
                                     math.ceil(height)))
        
        for simpl in tri.simplices:
            for aIdx, bIdx in ((simpl[0], simpl[1]),
                               (simpl[0], simpl[2]),
                               (simpl[1], simpl[2])):
                
                color = ((128, 192, 255) if((aIdx in downto and downto[aIdx] == bIdx) or
                                            (bIdx in downto and downto[bIdx] == aIdx) or
                                            (aIdx in downto and downto[aIdx] == None and bIdx in downto and downto[bIdx] == None)) else
                         (192//2, 255//2, 128//2))
                
                xDiff = points[bIdx][0] - points[aIdx][0]
                yDiff = points[bIdx][1] - points[aIdx][1]
                
                for part in range(101):
                    xDot = points[aIdx][0] + xDiff * (part / 100)
                    yDot = points[aIdx][1] + yDiff * (part / 100)
                    view.putpixel((int(xDot), int(yDot)), color)
                    pass
                pass
            pass
        
        for point in tri.points:
            view.putpixel((int(point[0]), int(point[1])), (255,255,255))
            pass

        view.show()
        pass
    
            
            
