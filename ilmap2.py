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

            span = int(math.sqrt(xDiff*xDiff + yDiff*yDiff) * 2 * scale + 1)
            for part in range(span + 1):
                xDot = tri.points[aIdx][0] + xDiff * (part / span)
                yDot = tri.points[aIdx][1] + yDiff * (part / span)
                view.putpixel((int(xDot * scale), int(yDot * scale)), color)
                pass
            pass
        pass

    for point in tri.points:
        view.putpixel((int(point[0] * scale), int(point[1] * scale)), (255,255,255))
        pass

def draw_elev(view, tri, height, scale):
    def ecolor(idx):
        return((255,  0,255) if idx not in height else
               (  0,  0,255) if height[idx] < -20 else
               (  0,128,255) if height[idx] < -10 else
               (  96,192,192) if height[idx] < 0 else
               (192,192, 96) if height[idx] < 10 else
               (128,256,128) if height[idx] < 100 else
               ( 64,224, 64) if height[idx] < 200 else
               (192,225,128) if height[idx] < 300 else
               (192,192, 64) if height[idx] < 400 else
               (224,128, 64) if height[idx] < 500 else
               (128, 64, 32) if height[idx] < 600 else
               (224, 64,  0) if height[idx] < 700 else
               (224,192, 192) if height[idx] < 800 else
               (256,256,256) if height[idx] < 984 else
               (256,  0,  0))
    
    for simpl in tri.simplices:
        for aIdx, bIdx in ((simpl[0], simpl[1]),
                           (simpl[0], simpl[2]),
                           (simpl[1], simpl[2])):
            aRed, aGreen, aBlue = ecolor(aIdx)
            bRed, bGreen, bBlue = ecolor(bIdx)
            xDiff = tri.points[bIdx][0] - tri.points[aIdx][0]
            yDiff = tri.points[bIdx][1] - tri.points[aIdx][1]
            redDiff = bRed - aRed
            greenDiff = bGreen - aGreen
            blueDiff = bBlue - aBlue
            
            span = int(math.sqrt(xDiff*xDiff + yDiff*yDiff) * 2 * scale + 1)
            for part in range(span+1):
                xDot = tri.points[aIdx][0] + xDiff * (part / span)
                yDot = tri.points[aIdx][1] + yDiff * (part / span)
                pcolor = (int(aRed + redDiff * (part / span)),
                          int(aGreen + greenDiff * (part / span)),
                          int(aBlue + blueDiff * (part / span)))
                
                view.putpixel((int(xDot * scale), int(yDot * scale)), pcolor)
                pass
            pass
        pass

    #for point in tri.points:
    #    view.putpixel((int(point[0] * scale), int(point[1] * scale)), (255,255,255))
    #    pass

def draw_heightmap(view, tri, height, scale):
    for x in range(1081):
        triX = x * scale
        for y in range(1081):
            triY = y * scale
            smplx = tri.find_simplex(((triX, triY),))
            if(smplx[0] == -1):
                color = 0
            else:
                aH, bH, cH = (height[idx] for idx in tri.simplices[smplx[0]])
                aP, bP, cP = (tri.points[idx] for idx in tri.simplices[smplx[0]])
                det = (bP[1]-cP[1])*(aP[0]-cP[0])+(cP[0]-bP[0])*(aP[1]-cP[1])
                aW = ((bP[1]-cP[1])*(triX-cP[0]) + (cP[0]-bP[0])*(triY-cP[1]))/det
                bW = ((cP[1]-aP[1])*(triX-cP[0]) + (aP[0]-cP[0])*(triY-cP[1]))/det
                cW = 1 - aW - bW
                triH = aH * aW + bH * bW + cH * cW
                color = int((triH + 40) / 1024 * 65536)
                if color < 0: color = 0
                if color > 65535: color = 65535
                view.putpixel((x, y), color)
                
def pdist(a, b):
    xdelt = a[0] - b[0]
    ydelt = b[1] - b[1]
    return(math.sqrt(xdelt * xdelt + ydelt * ydelt))
    
def tri_dist(tri, aIdx, bIdx):
    return(pdist(tri.points[aIdx], tri.points[bIdx]))

               
## Limits on aIdx height based on the edge to bIdx, returned as (hardmin, desired, hardmax)
def limitations(aIdx, bIdx, *, tri, downto, depth, height, gradefn, rimHeight=3, minHeight=-40):
    ## a is on the outer edge
    if(downto[aIdx] is None):  
        return(minHeight, minHeight, minHeight)

    ## a and b are river beds
    if(depth[aIdx] > 0 and depth[bIdx] > 0):
        ## a flows into b
        if(downto[aIdx] == bIdx):
            nomgrade = tri_dist(tri, aIdx, bIdx) * gradefn(depth[aIdx])
            return(height[bIdx] + nomgrade / 2,
                   height[bIdx] + nomgrade,
                   height[bIdx] + nomgrade * 2)

        ## b flows into a
        if(depth[bIdx] > 0 and downto[bIdx] == aIdx):
            nomgrade = tri_dist(tri, aIdx, bIdx) * gradefn(depth[bIdx])
            return(height[bIdx] - nomgrade * 2,
                   height[bIdx] - nomgrade,
                   height[bIdx] - nomgrade / 2)

        ## neither flow into the other
        return(None, None, None)
    
    ## a is a shore and a flows into b
    if(depth[aIdx] == 0 and downto[aIdx] == bIdx):
        return(height[bIdx] + tri_dist(tri, aIdx, bIdx) * 0.001,
               height[bIdx] + tri_dist(tri, aIdx, bIdx) * gradefn(depth[bIdx]),
               height[bIdx] + tri_dist(tri, aIdx, bIdx) * 0.60)
    ## b is a shore and b flows into a
    if(depth[bIdx] == 0 and downto[bIdx] == aIdx):
        return(height[bIdx] - tri_dist(tri, aIdx, bIdx) * 0.60,
               height[bIdx] - tri_dist(tri, aIdx, bIdx) * gradefn(depth[aIdx]),
               height[bIdx] - tri_dist(tri, aIdx, bIdx) * 0.001)

    ## a is a rim and a flows into b
    if(depth[aIdx] == -1 and downto[aIdx] == bIdx and downto[bIdx] is not None):
        return(max(height[downto[bIdx]] + rimHeight,
                   height[bIdx] + tri_dist(tri, aIdx, bIdx) * 0.001),
               height[bIdx] + tri_dist(tri, aIdx, bIdx) * gradefn(depth[aIdx]),
               None)

    ## b is a rim and b flows into a
    if(depth[bIdx] == -1 and downto[bIdx] == aIdx and downto[aIdx] is not None):
        return(None,
               None,  ## Note:  This is independent because b can be set independently of a
               height[bIdx] - tri_dist(tri, aIdx, bIdx) * 0.001)

    ## a is above b
    if(depth[aIdx] < depth[bIdx]):
        return(height[bIdx] + tri_dist(tri, aIdx, bIdx) * 0.001,
               height[bIdx] + tri_dist(tri, aIdx, bIdx) * gradefn(depth[aIdx]),
               height[bIdx] + tri_dist(tri, aIdx, bIdx) * 1.50)

    ## b is above a
    if(depth[aIdx] > depth[bIdx]):
        return(height[bIdx] - tri_dist(tri, aIdx, bIdx) * 1.50,
               height[bIdx] - tri_dist(tri, aIdx, bIdx) * gradefn(depth[bIdx]),
               height[bIdx] - tri_dist(tri, aIdx, bIdx) * 0.001)

    ## Otherwise, no restrictions
    return(None, None, None)

def oneraise(tri, downto, depth, height, gradefn, moveIndices, countsAsMove=0.1):
    allIdx = list(moveIndices)
    random.shuffle(allIdx);
    newheight = dict(height)
    updates = set()
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
            hooboy = (minH + maxH) / 2
            violations += 1
        else:
            hooboy = sum(desired) / len(desired)
            if(minH is not None and minH > hooboy):
                hooboy = minH
            if(maxH is not None and maxH < hooboy):
                hooboy = maxH

        newheight[aIdx] = hooboy * .95 + height[aIdx] * 0.05
        if(abs(newheight[aIdx] - height[aIdx]) > countsAsMove):
            updates.add(aIdx)
    
    return(newheight, violations, updates)

def fix_viols(tri, downto, height):
    allIdx = set(downto.keys());
    anyMove = set();
    while(allIdx):
        nextIdx = set()
        for aIdx in allIdx:
            bIdx = downto[aIdx]
            if bIdx == None:
                continue
            if height[aIdx] < height[bIdx]:
                delta = height[bIdx] - height[aIdx]
                height[aIdx] += delta / 2 + 0.005
                height[bIdx] -= delta / 2 + 0.005
                assert(height[aIdx] > height[bIdx])
                nextIdx.update(extend_region(tri, (aIdx,)).keys())
                anyMove.add(aIdx)
                anyMove.add(bIdx)
        allIdx = nextIdx
        
    return anyMove

def river_report(tri, idx, downto, depth, height):
    rlen = 0
    rfall = 0
    lastIdx = None
    while(True):
        if(lastIdx is not None):
            seglen = tri_dist(tri, lastIdx, idx)
            segfall = height[lastIdx] - height[idx]
            #print(f'  (falls {segfall}m along a {seglen}m run (slope: {segfall/seglen}))')
            rlen += seglen
            rfall += segfall

        #print(f'I({idx} H({height[idx]}) D({depth[idx]})')
        lastIdx, idx = idx, downto[idx]
        if(idx is None):
            break
    if(rlen != 0):
        print(f'=== TOTAL FALL: {rfall}m along a {rlen} run (slope: {rfall/rlen}) ===')
    else:
        print("=== SINGULAR ===")


class _IlMap_ut(unittest.TestCase):
    def test_randscatter(self):
        width = 18000
        height = 18000
        scale = 1/10
        dist = 16.666 * 4
        separate = 1900 / dist  ## about 1.9km apart
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

        viewH = PIL.Image.new('RGB', (math.ceil(width * scale),
                                      math.ceil(height * scale)))
        viewHM = PIL.Image.new('I;16', (1081, 1081));
        height = dict( (idx, -5) for idx in depth.keys())

        updates = set(height.keys())
        repeats = 100000
        for idx in range(repeats):
            height, violations, updates = oneraise(tri, downto, depth, height, gradefn, updates)
            if(idx % 7 == 0 and idx % 13 == 0 and idx % (7*13) != 0):
                updates.update(fix_viols(tri, downto, height))
            allheights = sorted(height.values())
            print(f"HEIGHT:  [{allheights[0]}:{allheights[-1]}] updated {len(updates)} with {violations} violations")

            if(len(updates) == 0 or idx % 200 == 199):
                for rIdx, rDep in depth.items():
                    if(rDep == 1):
                        river_report(tri, rIdx, downto, depth, height)

                tempHeight = dict(height)
                fix_viols(tri, downto, tempHeight)
                draw_elev(viewH, tri, tempHeight, scale)
                viewH.show()
                draw_heightmap(viewHM, tri, tempHeight, 18000 / 1081)
                viewHM.show()
                pass

            if(len(updates) == 0):
                break
            
            updates.update(extend_region(tri, updates).keys())
            updates.update(extend_region(tri, updates).keys())
            
        pass
    
            
            
