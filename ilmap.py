#!/usr/bin/env python3

import heapq
import math
import PIL.Image
import random
from scipy.spatial import Delaunay
import unittest

def rand_sum_to_n(summa, elements):
    for pick in range(elements-1, -1, -1):
        if pick==0:
            yield summa
        else:
            x = summa - random.uniform(0,summa**pick)**(1/pick)
            summa -= x
            yield x

def punctillate_rect(xmin, xmax, ymin, ymax, distsq):
    xdist = xmax-xmin
    ydist = ymax-ymin
    
    for idx in range(int(xdist * ydist / distsq) + 3):
        yield (random.uniform(xmin, xmax),
               random.uniform(ymin, ymax))
        
def in_radial_fn(pCenter, base, evenAmplSeq, oddAmplSeq):
    def _in_radial(point):
        xdel = point[0]-pCenter[0]
        ydel = point[1]-pCenter[1]
        r = math.sqrt(xdel*xdel + ydel*ydel)
        th = math.atan2(ydel, xdel)
        ampl = base
        for mul, ampls in enumerate(zip(evenAmplSeq, oddAmplSeq),
                                    start=1):
            freq = th*(mul)
            ampl += ampls[0] * math.cos(freq) + ampls[1] * math.sin(freq)

        return(r < ampl)
    
    return (_in_radial)

def make_bounded_rfn(pCenter, mindist, maxdist, *, overtones=17,
                     tiltangle=None, tiltspread=math.pi/30, istilttoward=True):
    radical = (maxdist + mindist) / 2
    baseamp = (maxdist - mindist) / 2
    thetamin = tiltangle if tiltangle is not None else -math.pi
    thetamax = tiltangle+tiltspread if tiltangle is not None else +math.pi
    
    evens = []
    odds = []
    for ampl in rand_sum_to_n(baseamp, overtones):
        if not istilttoward:
            ampl *= -1
        theta = random.uniform(thetamin, thetamax)
        evens.append(ampl * math.cos(theta))
        odds.append(ampl * math.sin(theta))

    return(in_radial_fn(pCenter, radical, evens, odds))

## NOTE:  Does not yield pTo, in accordance with python range behavior
def koch_path(pFrom, pTo, maxStep, bendFactor):
    pDel = (pTo[0]-pFrom[0], pTo[1]-pFrom[1])
    if pDel[0]*pDel[0] + pDel[1]*pDel[1] <= maxStep:
        yield pFrom
    else:
        localFactor = random.uniform(-bendFactor, bendFactor)
        pOrtho = (-pDel[1] * localFactor * math.sqrt(3) / 2, pDel[0] * localFactor * math.sqrt(3) / 2)
        pMid = ((pFrom[0] + pTo[0]) / 2 + pOrtho[0],
                (pFrom[1] + pTo[1]) / 2 + pOrtho[1])
        yield from koch_path(pFrom, pMid, maxStep, bendFactor)
        yield from koch_path(pMid, pTo, maxStep, bendFactor)

class IlMapper:
    def __init__(self, pointSeq):
        self._grid = Delaunay(pointSeq)
        self._nindptr, self._nindices = self._grid.vertex_neighbor_vertices
        self._forbidden_edges = set()
    
    def forbid_edge(self, aIdx, bIdx):
        if aIdx < bIdx:
            self._forbidden_edges.add((aIdx, bIdx))
        else:
            self._forbidden_edges.add((bIdx, aIdx))

    def forbid_long_edges(self, limit):
        limitsq = limit * limit
        for aIdx, aPoint in enumerate(self._grid.points):
            for bIdx in self._nindices[self._nindptr[aIdx]:self._nindptr[aIdx+1]]:
                if bIdx < aIdx:
                    continue
                bPoint = self._grid.points[bIdx]
                aDelt, bDelt = aPoint[0] - bPoint[0], aPoint[1] - bPoint[1]
                if aDelt*aDelt + bDelt*bDelt > limitsq:
                    self.forbid_edge(aIdx, bIdx)

        for aIdx, aPoint in enumerate(self._grid.points):
            if all(self.is_edge_forbidden(aIdx, bIdx) for bIdx in self._nindices[self._nindptr[aIdx]:self._nindptr[aIdx+1]]):
                raise(RuntimeError("Too many forbidden edges!"))
    
    def is_edge_forbidden(self, aIdx, bIdx):
        return((aIdx, bIdx) in self._forbidden_edges if aIdx < bIdx else
               (bIdx, aIdx) in self._forbidden_edges)
        
    ## Call one of the next three functions to set the shoreline based on grid features
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

    def set_radial_shore(self, *radial_select_fn_list):
        self._downto = dict()
        for pIdx, point in enumerate(self._grid.points):
            if not any(inside(point) for inside in radial_select_fn_list):
                self._downto[pIdx] = None

        self._depth = dict(self._downto)

    def set_koch_shore(self, pFrom, pTo, maxStep, bendFactor):
        self._downto = dict()
        for kPoint in koch_path(pFrom, pTo, maxStep, bendFactor):
            smplx = self._grid.find_simplex((kPoint,))
            if smplx[0] != -1:
                for pIdx in self._grid.simplices[smplx[0]]:
                    self._downto[pIdx] = None
        self._depth = dict(self._downto)
        
    ## region:     A container of point indices to extend
    ## dist:       An integer distance (in links) to extend
    ## exclusion:  A container of point indices which extension must not include
    ##---
    ## adjacents:  A dictionary of extension point indices to an index
    ##             tuple representing a path from region to adjacent
    def extend_region(self, region, dist=1, exclusion=()):
        adjacents = dict((r, (r,)) for r in region)
        adjacentweights = dict((r, 1) for r in adjacents.keys())
        exclusion = set(exclusion);
        
        while(dist > 0):
            exclusion.update(adjacents.keys());
            localadjacents = dict()
            localweights = dict()
            dist -= 1;

            for rIdx in adjacents.keys():
                for aIdx in self._nindices[self._nindptr[rIdx]:self._nindptr[rIdx+1]]:
                    if(self.is_edge_forbidden(aIdx, rIdx)):
                        continue
                    if(aIdx not in region and aIdx not in exclusion):
                        if aIdx in localadjacents:
                            pick = random.randrange(localweights[aIdx] + adjacentweights[rIdx])
                            if(pick < adjacentweights[rIdx]):
                                localadjacents[aIdx] = adjacents[rIdx] + (aIdx,)
                                pass
                            localweights[aIdx] += adjacentweights[rIdx]
                        else:
                            localadjacents[aIdx] = adjacents[rIdx] + (aIdx,)
                            localweights[aIdx] = adjacentweights[rIdx]
                            

            adjacents = localadjacents
            adjacentweights = localweights

        return(adjacents)

               
    ## Call this to create rivers separated by the given grid distance. If
    ## prune is given, prune river stubs which join into another rive in less
    ## than prune steps.
    def create_rivers(self, dist=1, mouths=None, prune=None):
        protshore = set(self._downto)
        if mouths is not None:
            for count in range(dist-1):
                protshore.update(self.extend_region(protshore))
                                                                        
        ## Create the riverine structures
        while(True):
            if mouths is None or mouths > 0:
                possibilities = self.extend_region(self._downto.keys(), dist);
            else:
                possibilities = self.extend_region((k for k, d in self._downto.items() if d is not None),
                                                   dist,
                                                   protshore)
            
            if mouths is not None and mouths > 0:
                mouths -= 1
            
            if(not possibilities):
                break

            pathIdx, path = random.choice(list(possibilities.items()))

            for pathIdx in range(1,len(path)):
                self._downto[path[pathIdx]] = path[pathIdx-1]
                self._depth[path[pathIdx]] = 0

        if prune:
            upto = dict()
            for idx, tgt in self._downto.items():
                if tgt is not None:
                    upto.setdefault(tgt, set()).add(idx)

            heads = list(idx for idx, tgt in self._downto.items() if tgt is not None and idx not in upto)
            random.shuffle(heads)
            for head in heads:
                parts = [head]
                for cnt in range(prune):
                    npart = self._downto[parts[-1]]
                    if npart is None or len(upto[npart]) > 1:
                        if npart is not None:
                            upto[npart].remove(parts[-1])
                        for idx in parts:
                            del self._depth[idx]
                            del self._downto[idx]
                        break
                    parts.append(npart)
    
        for rIdx in self._depth.keys():
            qIdx = rIdx
            while(qIdx is not None and
                  self._depth[qIdx] is not None):
                self._depth[qIdx] += 1
                qIdx = self._downto[qIdx]

    ## After creating rivers, levelize the remaining land
    def levelize_land(self, taper=0):
        exshore = 0 if taper == 0 else random.randrange(taper)
        shoredepths = set(idx for idx in self._depth.keys() if self._depth[idx] is None)
        riverdepths = sorted((idx for idx in self._depth.keys() if self._depth[idx] is not None),
                             key=lambda idx:self._depth[idx])
        level = 0
        while(True):
            exfrom = self._downto.keys()
            exclud = set()
            canstop = True
            if taper and riverdepths:
                riverdepths = riverdepths[:int(taper * len(riverdepths) / (taper + 1))]
                exfrom = set(exfrom)
                exfrom.difference_update(riverdepths)
                exclud.update(riverdepths)
                
                if taper > exshore:
                    exfrom.difference_update(shoredepths)
                    exclud.update(shoredepths)
                    
                taper-=1
                canstop = False
                    
            steppe = self.extend_region(exfrom, 1, exclud);
            if(canstop and not steppe):
                break
            for pIdx, path in steppe.items():
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
        futures = sorted((height, idx) for idx, height in self._height.items())
        heapq.heapify(futures)
        
        while(futures):
            pHeight, pIdx = heapq.heappop(futures)
            pPoint = self._grid.points[pIdx]
            for qIdx in self._nindices[self._nindptr[pIdx]:self._nindptr[pIdx+1]]:
                if self.is_edge_forbidden(pIdx, qIdx):
                    continue
                if qIdx in self._height:
                    continue
                qPoint = self._grid.points[qIdx]
                xdel = pPoint[0] - qPoint[0]
                ydel = pPoint[1] - qPoint[1]
                xydist = math.sqrt(xdel * xdel + ydel * ydel)                
                self._height[qIdx] = pHeight + slope_fn(self._depth[qIdx]) * xydist
                heapq.heappush(futures, (self._height[qIdx], qIdx))
            
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

            if lowest_setmax is None or highest_in_set is not None and highest_in_set < lowest_setmax:
                lowest_setmax = highest_in_set

        if(lowest_setmax is not None and lowest_setmax > 0):
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

    def force_one_side_shore(self, view2grid_fn):
        pointsets = list()
        omicron = 1081 * 2 / 9
        omega = 1081 * 7 / 9
        northset = []
        southset = []
        eastset = []
        westset = []

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
    
    def force_three_side_shore(self, view2grid_fn):
        pointsets = list()
        omicron = 1081 * 2 / 9
        omega = 1081 * 7 / 9
        northset = []
        southset = []
        eastset = []
        westset = []

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

        pointsets.append(northset + southset + eastset)
        pointsets.append(northset + southset + westset)
        pointsets.append(northset + eastset + westset)
        pointsets.append(southset + eastset + westset)
            
        self.force_one_underwater(pointsets)

    def height_compress(self, tgtmaxheight=1024-40, regionfn=None):
        maxheight = max(v for pIdx, v in self._height.items() if regionfn is None or regionfn(self._grid.points[pIdx]))
        if maxheight > tgtmaxheight:
            scale = tgtmaxheight / maxheight
            for pIdx, pHeight in self._height.items():
                if pHeight > 0:
                    self._height[pIdx] *= scale
    

    ## This takes a fully developed set of heights and reconstructs the set of
    ## depths and downtos based on flow.  The depths are normalized such that the same
    ## number of river nodes and shore nodes exist.
    def relevelize(self):
        ## This gives the number of nodes involved in a river
        numrivernodes = sum(1 for v in self._depth.values() if v is not None and v > 0)

        for pIdx, pHeight in self._height.items():
            if self._depth[pIdx] is None:
                continue
            steepest = None
            pPoint = self._grid.points[pIdx]
            self._depth[pIdx] = 0
            self._downto[pIdx] = None
            
            for qIdx in self._nindices[self._nindptr[pIdx]:self._nindptr[pIdx+1]]:
                if self.is_edge_forbidden(pIdx, qIdx):
                    continue
                qHeight = self._height[qIdx]
                qRise = pHeight - qHeight
                if(qRise < 0):
                    continue
                qPoint = self._grid.points[qIdx]
                qDel = (qPoint[0]-pPoint[0], qPoint[1]-pPoint[1])
                qRun = math.sqrt(qDel[0] * qDel[0] + qDel[1] * qDel[1])
                qSlope = qRise / qRun
                if(steepest is None or steepest < qSlope):
                    steepest = qSlope
                    self._downto[pIdx] = qIdx

            if self._downto[pIdx] is None:
                raise(RuntimeError("Created an unintended sink node!  How?"))

        for pIdx, pDepth in self._depth.items():
            if pDepth is None:
                continue
            qIdx = pIdx
            while(qIdx is not None and self._depth[qIdx] is not None):
                self._depth[qIdx] += 1
                qIdx = self._downto[qIdx]

        # Renormalize so the number of river nodes is the same (or nearly so)
        alldepths = sorted(v for v in self._depth.values() if v is not None)
        minriverflow = alldepths[-numrivernodes] - 1  # Rivers have flow of 1
        for pIdx in self._depth:
            if self._depth[pIdx] is not None:
                self._depth[pIdx] -= minriverflow

    def punch_rivers(self, d=6, src_d=0, minsegs=0, maxsegs=1):
        maxflow = max(v for v in self._depth.values() if v is not None)
        riverNodes = set()
        presumedDepth = dict()

        tapersegs = maxsegs - minsegs
        for segment in range(tapersegs):
            minsegflow = maxflow * (1 -  segment / tapersegs)
            additions = self.extend_region(riverNodes)
            for nIdx, path in additions.items():
                presumedDepth[nIdx] = presumedDepth[path[0]]
            
            riverNodes.update(additions)
            for rIdx, dp in self._depth.items():
                if dp is not None and dp > minsegflow:
                    riverNodes.add(rIdx)
                    presumedDepth[rIdx] = dp

        if(minsegs > 0):
            for rIdx, dp in self._depth.items():
                if dp is not None and dp > 0:
                    riverNodes.add(rIdx)
                    presumedDepth[rIdx] = dp

            additions = self.extend_region(riverNodes, dist=minsegs)
            for nIdx, path in additions.items():
                presumedDepth[nIdx] = presumedDepth[path[0]]
                    
            riverNodes.update(additions)
            
        for rIdx in riverNodes:
            rHeight = self._height[rIdx]
            rDepth = presumedDepth[rIdx]
            assert(rDepth > 0)
            
            fd = (d - src_d) * (1 - (1/2)**(rDepth - 1)) + src_d
            
            self._height[rIdx] -= (fd if rHeight >=0 else
                                   (40 + rHeight) / 40 * fd if rHeight >= -40 else
                                   0)

    def draw_heightmap(self, view, view2grid_fn):
        for x in range(view.width):
            for y in range(view.height):
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
                    if(0 < pXY[0] < view.width and
                       1 < pXY[1] < view.height):
                        pColor = tuple(int(weight * a + unweight * b) for a,b in zip(aColor, bColor))
                        view.putpixel(pXY, pColor)

            if(vertex_color is not None and
               0 < aXY[0] < view.width and
               0 < aXY[1] < view.height):                
                view.putpixel((int(aXY[0]), int(aXY[1])), vertex_color)

    def downto_color_fn(self, pathcolor=(255,255,128), offcolor=(128, 128, 255)):
        def _color(aIdx, bIdx):
            if((aIdx in self._downto and self._downto[aIdx] == bIdx) or
               (bIdx in self._downto and self._downto[bIdx] == aIdx)):
                return(pathcolor, pathcolor)
            else:
                return(offcolor, offcolor)

        return _color        
    
    def plan_color_fn(self, undeepest=None, levelCheck=True):
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
            elif(levelCheck and
                 self._depth[aIdx] <=0 and self._depth[bIdx] <=0 and
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

    pointsel = parser.add_mutually_exclusive_group()

    pointsel.add_argument('--shore', dest='pointsel', const='shore', action='store_const')
    pointsel.add_argument('--peninsula', dest='pointsel', const='peninsula', action='store_const')
    pointsel.add_argument('--island', dest='pointsel', const='island', action='store_const')
    pointsel.add_argument('--bay', dest='pointsel', const='bay', action='store_const')
    pointsel.add_argument('--strait', dest='pointsel', const='strait', action='store_const')
    
    parser.add_argument('--mapwidth', type=float, metavar='meters', default=18000)
    parser.add_argument('--mapheight', type=float, metavar='meters', default=18000)
    parser.add_argument('--viewwidth', type=int, metavar='pixels', default=1024)
    parser.add_argument('--viewheight', type=int, metavar='pixels', default=1024)
    parser.add_argument('--elevwidth', type=int, metavar='pixels', default=1081)
    parser.add_argument('--elevheight', type=int, metavar='pixels', default=1081)
    parser.add_argument('--nodeseparation', type=float, metavar='meters', default=50)
    parser.add_argument('--forbidedgefactor', type=float, default=5)
    parser.add_argument('--riverseparation', type=float, metavar='meters', default=900)
    parser.add_argument('--riverprune', type=float, metavar='meters')
    parser.add_argument('--rivermouths', type=int, metavar='mouths')
    parser.add_argument('--riverslopes', type=float, nargs='+', default=(0.01, 0.015,0.02, 0.03, 0.04))
    parser.add_argument('--riverdepth', type=float, metavar='meters', default=6)
    parser.add_argument('--riverwidth', type=float, metavar='meters', default=10)
    parser.add_argument('--mouthwidth', type=float, metavar='meters', default=150)
    parser.add_argument('--landslopes', type=float, nargs='+', default=(0.05, 0.09, 0.15, 0.30, 0.75))
    parser.add_argument('--landtaper', type=float, metavar='meters', default=0)
    parser.add_argument('--forceshore', action='store_true')
    parser.add_argument('--showshore', action='store_true')
    parser.add_argument('--showrivers', action='store_true')
    parser.add_argument('--showdownto', action='store_true')
    parser.add_argument('--showplan', action='store_true')
    parser.add_argument('--showrelevelize', action='store_true')
    parser.add_argument('--showheight', action='store_true')
    parser.add_argument('--output', metavar='filename', type=str)
    
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
    overflow = 1/9
    pminx = args.mapwidth * -overflow
    pmaxx = args.mapwidth * (1 + overflow)
    pminy = args.mapheight * -overflow
    pmaxy = args.mapheight * (1 + overflow)

    points.extend(punctillate_rect(pminx, pmaxx, pminy, pmaxy, distsq))

    print(f'POINTS (ANTE):  {len(points)}')
    mapper = IlMapper(points)
    points = None  ## Don't accidently use the original point set!
    print(f'POINTS (POST):  {len(mapper._grid.points)}')

    mapper.forbid_long_edges(args.nodeseparation * args.forbidedgefactor)

    print(f'{len(mapper._forbidden_edges)} edges forbidden!')
    
    
    radials = []
    mapninth = args.mapwidth / 9
    mapcenter = ((pmaxx+pminx)/2,
                 (pmaxy+pminy)/2)
    if(args.pointsel == 'island'):
        # Centers
        numIslands = random.randrange(1, 10)
        totArea = math.pi * mapninth * mapninth * 2.25 * 2.25
        for nIdx, isleArea in enumerate(rand_sum_to_n(totArea, numIslands)):
            isleMidRad = math.sqrt(isleArea / math.pi)
            isleMinRad = 0
            if isleMinRad >= mapninth * 1.25:
                isleMinRad = 0.5 * isleMidRad
            isleMaxRad = 2 * isleMidRad - isleMinRad
            isleMaxPos = 2.5 * mapninth - isleMaxRad
            isleMinPos = isleMaxPos * 0.8

            islePos = random.uniform(isleMinPos, isleMaxPos)
            isleTheta = math.tau * nIdx / numIslands
            radials.append(make_bounded_rfn((4.5 * mapninth + islePos * math.cos(isleTheta),
                                             4.5 * mapninth + islePos * math.sin(isleTheta)),
                                            isleMinRad, isleMaxRad))

        # Edges
        for centerPoint in ((4.5 * mapninth, 0),
                            (0, 4.5 * mapninth),
                            (4.5 * mapninth, args.mapwidth),
                            (args.mapwidth, 4.5 * mapninth)):
            radials.append(make_bounded_rfn(centerPoint,
                                            1.5 * mapninth,
                                            2.0 * mapninth))
       
    elif(args.pointsel == 'peninsula'):
        radials.append(make_bounded_rfn(mapcenter,
                                        args.mapwidth * 1 / 9,
                                        args.mapwidth * 7 / 9,
                                        tiltangle=random.uniform(-math.pi, +math.pi),
                                        tiltspread=random.uniform(math.pi / 360, math.pi / 60)))
    elif(args.pointsel == 'shore'):
        radials.append(make_bounded_rfn(mapcenter,
                                        args.mapwidth * 1 / 9,
                                        args.mapwidth * 11 / 9,
                                        istilttoward=False,
                                        tiltangle=random.uniform(-math.pi, +math.pi),
                                        tiltspread=random.uniform(math.pi / 360, math.pi / 60)))
    elif(args.pointsel == 'bay'):
        radials = [make_bounded_rfn(((pmaxx+pminx) * 1.5 / 9,
                                     (pmaxy+pminy) * 3.5 / 9),
                                    args.mapwidth * 3 / 9,
                                    args.mapwidth * 4 / 9),
                   make_bounded_rfn(((pmaxx+pminx) * 7.5 / 9,
                                     (pmaxy+pminy) * 3.5 / 9),
                                    args.mapwidth * 3 / 9,
                                    args.mapwidth * 4 / 9),
                   make_bounded_rfn(((pmaxy+pminy) / 2,
                                     0),
                                    args.mapwidth * 3.5 / 9,
                                    args.mapwidth * 5 / 9),
                   make_bounded_rfn((0,0),
                                    args.mapwidth * 4 / 9,
                                    args.mapwidth * 5 / 9),
                   make_bounded_rfn((args.mapwidth,0),
                                    args.mapwidth * 4 / 9,
                                    args.mapwidth * 5 / 9)]
    elif(args.pointsel == 'strait'):
        radius = max(abs(pmaxx - mapcenter[0]),
                     abs(pminx - mapcenter[0]),
                     abs(pmaxy - mapcenter[1]),
                     abs(pminy - mapcenter[1])) * math.sqrt(2)
        theta = random.uniform(-math.pi, math.pi)

        pDel = (radius * math.cos(theta), radius * math.sin(theta))
        pFrom = (mapcenter[0] - pDel[0], mapcenter[1] - pDel[1])
        pTo = (mapcenter[0] + pDel[0], mapcenter[1] + pDel[1])
        mapper.set_koch_shore(pFrom, pTo, 2 ** random.uniform(-1, 1) * args.nodeseparation,
                              random.uniform(0.2, 0.8))
        
    else:
        radials.append(make_bounded_rfn(mapcenter,
                                        args.mapwidth * 1.5 / 9,
                                        args.mapwidth * 10.5 / 9))

    if radials:
        mapper.set_radial_shore(*radials)
                          
    #mapper.set_boundary_shore()
    print(f'SHORE POINTS:   {len(mapper._depth)}')

    if(args.showshore):
        view = PIL.Image.new('RGB', viewXY)
        mapper.draw_grid(view, grid2view,
                         mapper.plan_color_fn(),
                         vertex_color=(128,128,128))
        view.show()
        
    mapper.create_rivers(math.ceil(args.riverseparation / args.nodeseparation),
                         prune=math.ceil(args.riverprune / args.nodeseparation) if args.riverprune else None,
                         mouths=args.rivermouths if args.rivermouths else None)
    print("Rivers are now flowing")
    
    if(args.showrivers):
        view = PIL.Image.new('RGB', viewXY)
        mapper.draw_grid(view, grid2view,
                         mapper.plan_color_fn(),
                         vertex_color=(128,128,128))
        view.show()

    mapper.levelize_land(taper=int(args.landtaper / args.nodeseparation))
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

    if(args.forceshore):
        mapper.force_one_tile_shore(view2elev)
        print("Shoreline encouraged")

    mapper.height_compress(regionfn=lambda p: (0 < p[0] < args.mapwidth) and (0 < p[1] < args.mapheight))
    mapper.relevelize()
    print("New levelization for new river flows")
    if(args.showrelevelize):
        view = PIL.Image.new('RGB', viewXY)
        mapper.draw_grid(view, grid2view,
                         mapper.plan_color_fn(levelCheck=False))
        view.show()

    
    mapper.punch_rivers(d=args.riverdepth, src_d=args.riverdepth/2,
                        minsegs=int(math.ceil(args.riverwidth/args.nodeseparation/2)),
                        maxsegs=int(math.ceil(args.mouthwidth/args.nodeseparation/2)))
    print("Rivers punched")
    
    if(args.showheight):
        view = PIL.Image.new('RGB', viewXY)
        mapper.draw_grid(view, grid2view,
                         mapper.elev_color_fn())
        view.show()
        
    if(1):
        view = PIL.Image.new('I;16', (args.elevwidth, args.elevheight));
        mapper.draw_heightmap(view, view2elev)
        if(args.output is None):
            view.show()
        else:
            view.save(args.output)
