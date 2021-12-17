#!/usr/bin/env python3

import unittest

## For our purposes:
##   pnt : point, simply a complex value
##   seg : line segment, simply a pair of points
##   tri : triangle, simpy a trio of points

def is_tri_rh(a_tri):
    ## Determine if the points of a_tri are oriented right-handedly
    return ((a_tri[1].imag - a_tri[0].imag) * (a_tri[2].real - a_tri[0].real) <
            (a_tri[2].imag - a_tri[0].imag) * (a_tri[1].real - a_tri[0].real))

def is_pnt_inside_tri(a_pnt, b_tri):
    orient = is_tri_rh((a_pnt, b_tri[0], b_tri[1]))
    return (is_tri_rh((a_pnt, b_tri[1], b_tri[2])) == orient and
            is_tri_rh((a_pnt, b_tri[2], b_tri[0])) == orient)

def does_intersect(a_seg, b_seg):
    ## The line A defined by two points intersects the segment B defined by two
    ## points if the triangle defined by A and b[0] and the triable defined by
    ## A and b[1] have opposite orientations.
    return(is_tri_rh((*a_seg, b_seg[0])) != is_tri_rh((*a_seg, b_seg[1])) and
           is_tri_rh((*b_seg, a_seg[0])) != is_tri_rh((*b_seg, a_seg[1])))

def zorder(seq):
    return(tuple(sorted(seq, key=lambda z:(z.real, z.imag))))

class Triangulation:
    def __init__(self):
        # Points:  v => { 'seg': set() }  ## Other side of linked segments        
        self._points = dict();

        # Segments:  (v, v) => { 'tri': set() } ## Third point of attached
        #                                       ## triangles.  Should be at most 2
        # (v,v) must be ordered
        self._segments = dict();

        # Triangles: (v, v, v) => dict()
        # Values must be ordered
        self._triangles = dict();

    def points(self):
        yield from self._points.keys()

    def segments(self):
        yield from self._segments.keys()

    def triangles(self):
        yield from self._triangles.keys()

    def hull_segs(self):
        for key, data in self._segments.item():
            if len(data['tri']) != 2:
                yield key
        
    def _add_point(self, pnt):
        self._points[pnt] = {'seg': set()}
        return pnt

    def _add_segment(self, a_pnt, b_pnt):
        seg = zorder((a_pnt, b_pnt))
        self._points[a_pnt]['seg'].add(b_pnt)
        self._points[b_pnt]['seg'].add(a_pnt)
        self._segments[seg] = { 'tri': set() }
        return seg
        

    def _add_triangle(self, a_pnt, b_pnt, c_pnt):
        tri = zorder((a_pnt, b_pnt, c_pnt))
        self._segments[(tri[0], tri[1])]['tri'].add(tri[2])
        self._segments[(tri[0], tri[2])]['tri'].add(tri[1])
        self._segments[(tri[1], tri[2])]['tri'].add(tri[0])
        self._triangles[tri] = dict()
        return tri

    def _del_triangle(self, tri):
        del self._triangles[tri]
        self._segments[(tri[0], tri[1])]['tri'].remove(tri[2])
        self._segments[(tri[0], tri[2])]['tri'].remove(tri[1])
        self._segments[(tri[1], tri[2])]['tri'].remove(tri[0])

    def _del_segment(self, seg):
        del self._segments[seg]
        self._points[seg[0]].remove(seg[1])
        self._points[seg[1]].remove(seg[0])
    
    def add_point(self, pnt):
        self._add_point(pnt);
        
        checksegs = list()
        for tri in self._triangles.keys():
            if(is_pnt_inside_tri(pnt, tri)):
                self._del_triangle(tri)
                checksegs.append(self._add_segment(pnt, tri[0]))
                checksegs.append(self._add_segment(pnt, tri[1]))
                checksegs.append(self._add_segment(pnt, tri[2]))

                self._add_triangle(pnt, tri[0], tri[1])
                self._add_triangle(pnt, tri[1], tri[2])
                self._add_triangle(pnt, tri[2], tri[0])

                break
            else:
                pass  ## safety
        else:
            ## Not inside any triangle.
            ## BOZO
            
                
                
            

class ut_ilmap(unittest.TestCase):
    def test_orient(self):
        self.assertTrue(is_tri_rh((0, 1, 2+1j)))
        self.assertTrue(is_tri_rh((0, 1, 1+1j)))
        self.assertTrue(is_tri_rh((0, 1, 0.5 + 1j)))
        self.assertTrue(is_tri_rh((0, 1, 0 + 1j)))
        self.assertTrue(is_tri_rh((0, 1, -1 + 1j)))
        self.assertTrue(is_tri_rh((1, 2+1j, -3+2j)))
        self.assertFalse(is_tri_rh((0, 1, 1-1j)))

    def test_inside(self):
        self.assertTrue(is_pnt_inside_tri(1.5+1j, (0, 2, 2+2j)))
        self.assertFalse(is_pnt_inside_tri(3+1j, (0, 2, 2+2j)))
        self.assertFalse(is_pnt_inside_tri(3+3j, (0, 2, 2+2j)))

    def test_intersect(self):
        self.assertTrue(does_intersect((0, 2), (1-1j, 1+1j)))
        self.assertFalse(does_intersect((0, 0+2j), (1-1j, 1+1j)))
        
