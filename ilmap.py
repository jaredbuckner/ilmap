#!/usr/bin/env python3

import heapq
import math
import PIL.Image
import random
from scipy.spatial import Delaunay
import unittest

class IlMapper:
    def __init__(self, pointSeq):
        self._grid = Delaunay(pointSeq)
        self._nindptr, self._nindices = self._grid.vertex_neighbor_vertices

    ## Call one of the next two functions to set the shoreline based on grid features
    def set_hull_shore(self):
        self._downto = dict();
        for face in self._grid.convex_hull:
            for pIdx in face:
                self._downto[pIdx] = None
        self._depth = dict(self._downto)

    def set_boundary_shore(self):
        self._downto = dict()
        for triIdx, tri in enumerate(self._grid.simplices):
            if(-1 in self._grid.neighbors[triIdx]):
                for pIdx in tri:
                    self._downto[pIdx] = None
        self._depth = dict(self._downto)
        
    ## region:     A container of point indices to extend
    ## dist:       An integer distance (in links) to extend
    ## exclusion:  A container of point indices which extension must not include
    ##---
    ## adjacents:  A dictionary of extension point indices to a list of index
    ##             tuples representing a path from region to adjacent
    def extend_region(self, region, dist=1, exclusion=()):
        adjacents = dict((r, [(r,)]) for r in region)
        exclusion = set(exclusion);
        
        while(dist > 0):
            exclusion.update(adjacents.keys());
            localadjacents = dict();
            dist -= 1;

            for rIdx in adjacents.keys():
                for aIdx in self._nindices[self._nindptr[rIdx]:self._nindptr[rIdx+1]]:
                    if(aIdx not in region and aIdx not in exclusion):
                        localadjacents.setdefault(aIdx, list())
                        for path in adjacents[rIdx]:
                            localadjacents[aIdx].append(path + (aIdx,))

            adjacents = localadjacents

        return(adjacents)

    ## Call this to create rivers separated by the given grid distance
    def create_rivers(self, dist=1, mouths=None):
       ## Create the riverine structures
        while(True):
            
            possibilities = self.extend_region(self._downto.keys(), dist);

            if(mouths is not None):
                if(mouths > 0):
                    mouths -= 1;
                else:
                    newposs = dict()
                    for rIdx, paths in possibilities.items():
                        newpaths = []                        
                        for path in paths:
                            if path[0] in self._downto and self._downto[path[0]] is not None:
                                newpaths.append(path)
                        if newpaths:
                            newposs[rIdx] = newpaths
                    possibilities = newposs
            
            if(not possibilities):
                break
            paths = random.choice(tuple(possibilities.values()));
            path = random.choice(paths)

            for pathIdx in range(1,len(path)):
                self._downto[path[pathIdx]] = path[pathIdx-1]
                self._depth[path[pathIdx]] = 0

    
        for rIdx in self._depth.keys():
            qIdx = rIdx
            while(qIdx is not None and
                  self._depth[qIdx] is not None):
                self._depth[qIdx] += 1
                qIdx = self._downto[qIdx]

    ## After creating rivers, levelize the remaining land
    def levelize_land(self):
        level = 0
        while(True):
            steppe = self.extend_region(self._downto.keys(), 1);
            if(not steppe):
                break
            for pIdx, paths in steppe.items():
                path = random.choice(paths)
                self._downto[pIdx] = path[0]
                self._depth[pIdx] = level
            level -= 1

    ## Make a slope function suitable for elevating    
    def slope_fn(self,
                 river_slopes = (0.01, 0.015,0.02, 0.03, 0.04),
                 land_slopes = (0.05, 0.09, 0.15, 0.30, 0.75)):
        ordered = sorted(v for v in self._depth.values() if v is not None)
        leastdepth = ordered[0]
        ordered = [v for v in ordered if v >0]
        depthbits = tuple(ordered[i * len(ordered) // len(river_slopes)] for i in range(len(river_slopes)))

        def _slope(d):
            for didx, bit in enumerate(reversed(depthbits)):
                if d >= bit:
                    return(river_slopes[-didx])
            for didx, slope in enumerate(land_slopes):
                if d >= leastdepth * (didx+1) // len(land_slopes):
                    return slope

            return land_slopes[-1]
        
        return _slope

    ## Get desired, min branch height function based on a slope_fn
    def branch_height_fn(self, slope_fn, bank_min=10):
        def _branch_height(pIdx):
            qIdx = self._downto.get(pIdx)
            if(qIdx is None):
                return None

            qHeight = self._height.get(qIdx)
            if(qHeight is None):
                return None

            pDepth = self._depth[pIdx]
            sDepth = pDepth
            if(pDepth == 0):
                qDepth = self._depth[qIdx]
                if qDepth is not None:
                    sDepth = qDepth

            pPoint = self._grid.points[pIdx]
            qPoint = self._grid.points[qIdx]
            xdel = pPoint[0] - qPoint[0]
            ydel = pPoint[1] - qPoint[1]
            xydist = math.sqrt(xdel * xdel + ydel * ydel)
            pHeight = qHeight + slope_fn(sDepth) * xydist

            if(pDepth == -1):
                rIdx = self._downto[qIdx]
                if rIdx is not None:
                    rHeight = self._height[rIdx]
                    if(pHeight < rHeight + bank_min):
                        pHeight = rHeight + bank_min
                        
            return(pHeight, qHeight)

        return(_branch_height)
    
    ## Initialize shore heights
    def init_shore_heights(self, variance=0.1):
        self._height = dict((idx, -40 + random.uniform(-variance, variance))
                            for idx in self._depth.keys() if self._depth[idx] is None)

    def gen_method_heights(self, branch_height_fn, slope_min=0.001, slope_max=0.8):
        futures = []
        for pIdx in self._height.keys():
            for qIdx in self._nindices[self._nindptr[pIdx]:self._nindptr[pIdx+1]]:
                if(self._downto[qIdx] == pIdx):
                    desired, minheight = branch_height_fn(qIdx)
                    futures.append((self._height[pIdx], desired, minheight, qIdx))

        futures.sort(reverse=True)

        while(futures):
            dummy, desired, minheight, pIdx = futures.pop()
            
            pPoint = self._grid.points[pIdx]

            if(self._depth[pIdx] < -1):
                minheight = -45  ## I don't care so much about non-rivers.  Let's see what happens...
            
            for qIdx in self._nindices[self._nindptr[pIdx]:self._nindptr[pIdx+1]]:
                qHeight = self._height.get(qIdx)
                if qHeight is None:
                    continue
                qDepth = self._depth[qIdx]
                qPoint = self._grid.points[qIdx]
                xdel = pPoint[0] - qPoint[0]
                ydel = pPoint[1] - qPoint[1]
                xydist = math.sqrt(xdel * xdel + ydel * ydel)
                qmin = (qHeight + slope_min * xydist if qDepth is not None and qDepth >=0 else
                        qHeight - slope_max * xydist)
                qmax = qHeight + slope_max * xydist
                if(desired > qmax):
                    desired = qmax
                if(minheight < qmin):
                    minheight = qmin

            self._height[pIdx] = (desired if desired > minheight else minheight)
            for qIdx in self._nindices[self._nindptr[pIdx]:self._nindptr[pIdx+1]]:
                if(self._downto[qIdx] == pIdx):
                    desired, minheight = branch_height_fn(qIdx)
                    futures.append((self._height[pIdx], desired, minheight, qIdx))
            
            futures.sort(reverse=True)
                       
    
    def gen_sweep_heights(self, slope_fn):
        futures = sorted((idx for idx in self._height.keys()),
                         key=lambda idx: self._height[idx],
                         reverse=True)

        while(futures):
            pIdx = futures.pop()
            pHeight = self._height[pIdx]
            pPoint = self._grid.points[pIdx]
            for qIdx in self._nindices[self._nindptr[pIdx]:self._nindptr[pIdx+1]]:
                if qIdx in self._height:
                    continue
                qPoint = self._grid.points[qIdx]
                xdel = pPoint[0] - qPoint[0]
                ydel = pPoint[1] - qPoint[1]
                xydist = math.sqrt(xdel * xdel + ydel * ydel)                
                self._height[qIdx] = pHeight + slope_fn(self._depth[qIdx]) * xydist
                futures.append(qIdx)
            
            futures.sort(key=lambda idx: self._height[idx],
                         reverse=True)

    def force_one_underwater(self, sets_of_points):
        lowest_setmax = None
        for pSet in sets_of_points:
            highest_in_set = None
            for aPoint in pSet:
                smplx = self._grid.find_simplex((aPoint,))
                if(smplx[0] != -1):
                    for pIdx in self._grid.simplices[smplx[0]]:
                        if highest_in_set is None or self._height[pIdx] > highest_in_set:
                            highest_in_set = self._height[pIdx]

            if lowest_setmax is None or highest_in_set < lowest_setmax:
                lowest_setmax = highest_in_set

        if(lowest_setmax > 0):
            for pIdx in self._height.keys():
                self._height[pIdx] -= lowest_setmax

    def force_one_tile_shore(self, view2grid_fn):
        pointsets = list()
        omicron = 1081 * 2 / 9
        omega = 1081 * 7 / 9
        for piece in range(2, 7):
            alpha = 1081 * piece / 9
            beta = 1081 * (piece + 1) / 9
            northset = []
            southset = []
            eastset = []
            westset = []
            for part in range(101):
                weight = part / 100
                unweight = 1 - weight
                kappa = weight * alpha + unweight * beta
                northset.append(view2grid_fn((omicron, kappa)))
                southset.append(view2grid_fn((omega, kappa)))
                eastset.append(view2grid_fn((kappa, omicron)))
                westset.append(view2grid_fn((kappa, omega)))

            pointsets.append(northset)
            pointsets.append(southset)
            pointsets.append(eastset)
            pointsets.append(westset)
            
        self.force_one_underwater(pointsets)
    
    def punch_rivers(self, d=6, src_d=6):
        ## NOTE:  Constructing a list only for counting the elements is bad
        ## from a memory perspective.  However, it turns out in python this is
        ## faster that sum(1 for v in <generator>) if the number of items is
        ## much smaller than the available memory footprint.  Because these
        ## maps are only a million nodes or so, this is reasonable.

        ## This gives the number of nodes involved in a river
        riverlines = sum(1 for v in self._depth.values() if v is not None and v > 0)

        flowsize = dict()
        flowdown = dict()
        for pIdx, pHeight in self._height.items():            
            steepest = None
            flowsize[pIdx] = 0
            flowdown[pIdx] = None
            for qIdx in self._nindices[self._nindptr[pIdx]:self._nindptr[pIdx+1]]:
                qHeight = self._height[qIdx]
                qDelt = pHeight - qHeight
                if(qDelt > 0 and (steepest is None or steepest < qDelt)):
                    steepest = qDelt
                    flowdown[pIdx] = qIdx

        for pIdx in flowsize.keys():
            qIdx = pIdx
            while(qIdx is not None):
                flowsize[qIdx] += 1
                qIdx = flowdown[qIdx]

        rivers = sorted(flowsize.keys(), key=lambda idx:flowsize[idx], reverse=True)
        riverbeds = set(rivers[:riverlines])
        rivers = set(riverbeds)
        rivers.update(self.extend_region(rivers))
        
        for rIdx in rivers:
            rHeight = self._height[rIdx]
            md = d if rIdx in riverbeds else src_d
            
            self._height[rIdx] -= (md if rHeight >=0 else
                                   (40 + rHeight) / 40 * md if rHeight >= -40 else
                                   0)
                

    def draw_heightmap(self, view, view2grid_fn):
        for x in range(1081):
            for y in range(1081):
                aXY = (x, y)
                aPoint = view2grid_fn((x, y))
                smplx = self._grid.find_simplex((aPoint,))
                if(smplx[0] == -1):
                    color = 0
                else:
                    aH, bH, cH = (self._height[idx] for idx in self._grid.simplices[smplx[0]])
                    aP, bP, cP = (self._grid.points[idx] for idx in self._grid.simplices[smplx[0]])
                    det = (bP[1]-cP[1])*(aP[0]-cP[0])+(cP[0]-bP[0])*(aP[1]-cP[1])
                    aW = ((bP[1]-cP[1])*(aPoint[0]-cP[0]) + (cP[0]-bP[0])*(aPoint[1]-cP[1]))/det
                    bW = ((cP[1]-aP[1])*(aPoint[0]-cP[0]) + (aP[0]-cP[0])*(aPoint[1]-cP[1]))/det
                    cW = 1 - aW - bW
                    self._gridH = aH * aW + bH * bW + cH * cW
                    color = int((self._gridH + 40) / 1024 * 65536)
                    if color < 0: color = 0
                    if color > 65535: color = 65535
                    view.putpixel((x, y), color)

    
    def draw_grid(self, view, grid2view_fn, color_fn, vertex_color=None):
        for aIdx, aPoint in enumerate(self._grid.points):
            aXY = grid2view_fn(aPoint)
            
            for bIdx in self._nindices[self._nindptr[aIdx]:self._nindptr[aIdx+1]]:
                if bIdx < aIdx:
                    continue

                bXY = grid2view_fn(self._grid.points[bIdx])
                aColor, bColor = color_fn(aIdx, bIdx)

                xspan, yspan = bXY[0]-aXY[0], bXY[1]-aXY[1]
                ispan = int(4 * max(abs(xspan), abs(yspan)) + 1)
                for part in range(ispan + 1):
                    weight = part / ispan
                    unweight = 1 - weight
                    pXY = (int(weight * aXY[0] + unweight * bXY[0]),
                           int(weight * aXY[1] + unweight * bXY[1]))

                    pColor = tuple(int(weight * a + unweight * b) for a,b in zip(aColor, bColor))
                    view.putpixel(pXY, pColor)

            if vertex_color is not None:
                view.putpixel((int(aXY[0]), int(aXY[1])), vertex_color)

    def downto_color_fn(self, pathcolor=(255,255,128), offcolor=(128, 128, 255)):
        def _color(aIdx, bIdx):
            if((aIdx in self._downto and self._downto[aIdx] == bIdx) or
               (bIdx in self._downto and self._downto[bIdx] == aIdx)):
                return(pathcolor, pathcolor)
            else:
                return(offcolor, offcolor)

        return _color        
    
    def plan_color_fn(self, undeepest=None):
        if(undeepest is None):
            undeepest = min((v for v in self._depth.values() if v is not None),
                            default=None)

        def _pcolor(aIdx, bIdx):
            if(aIdx not in self._depth or bIdx not in self._depth):
                return (192, 0, 0)
            elif(self._depth[aIdx] is None and self._depth[bIdx] is None):
                return (0, 0, 255)
            elif(self._depth[aIdx] is None or self._depth[bIdx] is None):
                return(64, 64, 255)
            elif(self._depth[aIdx] <=0 and self._depth[bIdx] <=0 and
                  abs(self._depth[aIdx] - self._depth[bIdx]) > 1):
                return (255, 0, 255)
            elif(self._depth[aIdx] > 0 and self._depth[bIdx] > 0 and
                  (self._downto[aIdx] == bIdx or self._downto[bIdx] == aIdx)):
                return (128, 192, 255)
            elif(self._depth[aIdx] >= 0 and self._depth[bIdx] >= 0):
                return (96, 128, 255)
            elif(self._depth[aIdx] >= (undeepest // 5) and self._depth[bIdx] >= (undeepest // 5)):
                return (192 * 3//4, 255 * 3 //4 , 128 * 3 // 4)
            elif(self._depth[aIdx] >= (undeepest * 2 // 5) and self._depth[bIdx] >= (undeepest * 2 // 5)):
                return (192//2, 255//2, 128//2)
            elif(self._depth[aIdx] >= (undeepest * 3 // 5) and self._depth[bIdx] >= (undeepest * 3 // 5)):
                return (255//2, 192//2, 128//2)
            elif(self._depth[aIdx] >= (undeepest * 4 // 5) and self._depth[bIdx] >= (undeepest * 4 // 5)):
                return (255 * 3//4, 192*3//4, 128*3//4)
            else:
                return (255, 255, 255)

        def _colorfn(aIdx, bIdx):
            color = _pcolor(aIdx, bIdx)
            return (color, color)
        
        return _colorfn
    
    def elev_color_fn(self):
        def _ecolor(idx):
            return((255,  0,255) if idx not in self._height else
                   (  0,  0,255) if self._height[idx] < -20 else
                   (  0,128,255) if self._height[idx] < -10 else
                   (  96,192,192) if self._height[idx] < 0 else
                   (192,192, 96) if self._height[idx] < 10 else
                   (128,256,128) if self._height[idx] < 100 else
                   ( 64,224, 64) if self._height[idx] < 200 else
                   (192,225,128) if self._height[idx] < 300 else
                   (192,192, 64) if self._height[idx] < 400 else
                   (224,128, 64) if self._height[idx] < 500 else
                   (128, 64, 32) if self._height[idx] < 600 else
                   (224, 64,  0) if self._height[idx] < 700 else
                   (224,192, 192) if self._height[idx] < 800 else
                   (256,256,256) if self._height[idx] < 984 else
                   (256,  0,  0))

        def _colorfn(aIdx, bIdx):
            return(_ecolor(aIdx), _ecolor(bIdx))

        return _colorfn

if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(
        description='Triangulation Map Generator');

    parser.add_argument('--mapwidth', type=float, metavar='meters', default=18000)
    parser.add_argument('--mapheight', type=float, metavar='meters', default=18000)
    parser.add_argument('--viewwidth', type=int, metavar='pixels', default=1024)
    parser.add_argument('--viewheight', type=int, metavar='pixels', default=1024)
    parser.add_argument('--elevwidth', type=int, metavar='pixels', default=1081)
    parser.add_argument('--elevheight', type=int, metavar='pixels', default=1081)
    parser.add_argument('--nodeseparation', type=float, metavar='meters', default=50)
    parser.add_argument('--riverseparation', type=float, metavar='meters', default=900)
    parser.add_argument('--rivermouths', type=int, metavar='mouths')
    parser.add_argument('--riverslopes', type=float, nargs='+', default=(0.01, 0.015,0.02, 0.03, 0.04))
    parser.add_argument('--riverdepth', type=float, metavar='meters', default=6)
    parser.add_argument('--landslopes', type=float, nargs='+', default=(0.05, 0.09, 0.15, 0.30, 0.75))
    parser.add_argument('--showshore', action='store_true')
    parser.add_argument('--showrivers', action='store_true')
    parser.add_argument('--showdownto', action='store_true')
    parser.add_argument('--showplan', action='store_true')
    parser.add_argument('--showheight', action='store_true')

    args = parser.parse_args()

    grid2viewscale = (args.viewwidth / args.mapwidth,
                      args.viewheight / args.mapheight)
    def grid2view(point):
        return(point[0] * grid2viewscale[0],
               point[1] * grid2viewscale[1])

    elev2gridscale = (args.mapwidth / args.elevwidth,
                      args.mapheight / args.elevheight)
    def view2elev(point):
        return(point[0] * elev2gridscale[0],
               point[1] * elev2gridscale[1])

    # dist is now args.nodeseparation
    # separate is math.ceil(args.riverseparation / args.nodeseparation) ???
    
    viewXY = (args.viewwidth, args.viewheight)

    ## Make some points
    points = list()
    distsq = args.nodeseparation * args.nodeseparation
    for idx in range(int(args.mapwidth * args.mapheight / distsq) + 3):
        newpoint = (random.uniform(0, args.mapwidth), random.uniform(0, args.mapheight));
        points.append(newpoint)

    print(f'POINTS (ANTE):  {len(points)}')
    mapper = IlMapper(points)
    points = None  ## Don't accidently use the original point set!
    print(f'POINTS (POST):  {len(mapper._grid.points)}')

    mapper.set_boundary_shore()
    print(f'SHORE POINTS:   {len(mapper._depth)}')

    if(args.showshore):
        view = PIL.Image.new('RGB', viewXY)
        mapper.draw_grid(view, grid2view,
                         mapper.plan_color_fn(),
                         vertex_color=(128,128,128))
        view.show()

    mapper.create_rivers(math.ceil(args.riverseparation / args.nodeseparation),
                         mouths=args.rivermouths if args.rivermouths else None)
    print("Rivers are now flowing")
    
    if(args.showrivers):
        view = PIL.Image.new('RGB', viewXY)
        mapper.draw_grid(view, grid2view,
                         mapper.plan_color_fn(),
                         vertex_color=(128,128,128))
        view.show()

    mapper.levelize_land()
    print("Land levelized")

    if(args.showdownto):
        view = PIL.Image.new('RGB', viewXY)
        mapper.draw_grid(view, grid2view,
                         mapper.downto_color_fn(),
                         vertex_color=(128,128,128))
        view.show()

    if(args.showplan):
        view = PIL.Image.new('RGB', viewXY)
        mapper.draw_grid(view, grid2view,
                         mapper.plan_color_fn(),
                         vertex_color=(128,128,128))
        view.show()

    mapper.init_shore_heights()
    print("Shore lowered")
    
    #mapper.gen_method_heights(mapper.branch_height_fn(mapper.slope_fn()))
    mapper.gen_sweep_heights(mapper.slope_fn(river_slopes = args.riverslopes,
                                             land_slopes = args.landslopes))
        
    print("Land heightmapped")
    
    mapper.force_one_tile_shore(view2elev)
    print("Shoreline encouraged")
    
    mapper.punch_rivers(d=args.riverdepth, src_d=args.riverdepth/2)
    
    print("Rivers punched")
    
    if(args.showheight):
        view = PIL.Image.new('RGB', viewXY)
        mapper.draw_grid(view, grid2view,
                         mapper.elev_color_fn())
        view.show()
        
    if(1):
        view = PIL.Image.new('I;16', (args.elevwidth, args.elevheight));
        mapper.draw_heightmap(view, view2elev)
        view.show()
