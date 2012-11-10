#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
import numpy as n, pylab as p, time

# simple routines for calculating convex hull given a set of points in 2D

# calculate angle in 2-D between points and x axis
def _angle_to_point(point, centre):
    delta = point - centre
    res = n.arctan(delta[1] / delta[0])
    if delta[0] < 0:
        res += n.pi
    return res

def _draw_triangle(p1, p2, p3, **kwargs):
    tmp = n.vstack((p1,p2,p3))
    x,y = [x[0] for x in zip(tmp.transpose())]
    p.fill(x,y, **kwargs)

# calculate area of any triangle given co-ordinates of the corners in p1, p2 and p3
def area_of_triangle(p1, p2, p3):
    return n.linalg.norm(n.cross((p2 - p1), (p3 - p1)))/2.


# Calculate subset of points that make a convex hull around points
# Recursively eliminates points that lie inside two neighbouring points
# until only convex hull is remaining.
# :Parameters:
#    points : ndarray (2 x m)
#        array of points for which to find hull
#
#:Returns:
#    hull_points : ndarray (2 x n)
#        convex hull surrounding points
def convex_hull(points):
    n_pts = points.shape[1]
    assert(n_pts > 5)
    centre = points.mean(1)
    angles = n.apply_along_axis(_angle_to_point, 0, points, centre)
    pts_ord = points[:,angles.argsort()]
    pts = [x[0] for x in zip(pts_ord.transpose())]
    prev_pts = len(pts) + 1
    k = 0
    while prev_pts > n_pts:
        prev_pts = n_pts
        n_pts = len(pts)
        i = -2
        while i < (n_pts - 2):
            Aij = area_of_triangle(centre, pts[i],     pts[(i + 1) % n_pts])
            Ajk = area_of_triangle(centre, pts[(i + 1) % n_pts], \
                                   pts[(i + 2) % n_pts])
            Aik = area_of_triangle(centre, pts[i],     pts[(i + 2) % n_pts])
            if Aij + Ajk < Aik:
                del pts[i+1]
            i += 1
            n_pts = len(pts)
        k += 1
    return n.asarray(pts)