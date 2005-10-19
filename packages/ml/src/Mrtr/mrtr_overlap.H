/*
#@HEADER
# ************************************************************************
#
#                 Copyright (2002) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2, or (at your option)
# any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.
#
# Questions? Contact Jonathan Hu (jhu@sandia.gov) or Ray Tuminaro 
# (rstumin@sandia.gov).
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#ifdef TRILINOS_PACKAGE

#ifndef MRTR_OVERLAP_H
#define MRTR_OVERLAP_H

#include <ctime>
#include <iostream>
#include <iomanip>
#include <map>

#include "mrtr_point.H"
#include "mrtr_edge.H"


using namespace std;

// ----------   User Defined Includes   ----------

namespace MRTR
{

// forward declarations
class MRTR::Interface;
class MRTR::Segment;
class MRTR::Node;

class  Overlap 
{
public:

  enum Polygon_Type
  {
    poly_none,
    poly_1_2_same,
    poly_1_2_dif,
    poly_1_4_same,
    poly_1_4_dif,
    poly_2_2_same,
    poly_2_2_dif,
    poly_3_0,
    poly_0_0,
    poly_0_4_same,
    poly_0_4_dif,
    poly_0_2
  };
  
  // ctor
  explicit Overlap(MRTR::Segment& sseg, MRTR::Segment& mseg, MRTR::Interface& inter);
  
  // dtor
  virtual ~Overlap();
  
  // compute overlap between 2 segments and it's triangulation
  bool ComputeOverlap();
  
private:  
  // don't want = operator
  Overlap operator = (const Overlap& old);
  // don't want copy-ctor
  Overlap(MRTR::Overlap& old);

private:

  // build line information from triangles
  bool build_lines();
  // build projection of master nodes onto slave segment
  bool build_mxi();
  // build projection of slave nodes onto master segment
  bool build_sxi();
  // build the outward normal to the sseg
  bool build_normal();
  
  // perform clipping algorithm
  bool Clipelements();
  // test whether a point is inside or outside of a clip edge
  bool Clip_TestPoint(const double* N, const double* PE, const double* P);
  // find intersection of clipping edge with line
  bool Clip_Intersect(double* N,double* PE,double* P0,double* P1,double* xi);

  // add a point to the polygon
  bool Clip_AddPointtoPolygon(int id, double* P);
  // remove a point from the polygon
  bool Clip_RemovePointfromPolygon(const int id,const double* P);
  // get the size of the nodal polygon
  int Clip_SizePolygon() { return p_.size(); }
  // get view of point in polygon, calling routine is responsible for freeing
  MRTR::Point** Clip_PointView();

  // add edge to polygon
  bool Clip_AddEdgetoPolygon(int id, int node0, int node1);
  // get view from an edge
  MRTR::Edge* GetEdgeViewfromPolygon(int id);
  
private:

MRTR::Interface&    inter_;
MRTR::Segment&      sseg_;
MRTR::Segment&      mseg_;

bool                overlap_;   // true if segments overlap
int                 nnode_;     // # of master nodes

double              mxi_[3][2]; // local coords of mnodes projected onto sseg
bool                havemxi_;   // flag indicating whether mxi_ and min_ have been build

double              sxi_[3][2]; // local coords of snodes 
double              sn_[3][2];   // outward normal to sseg's edges in local coords
bool                havesxi_;   // flag indicating whether sxi_ and sin_ have been build

double              sline_[3][4]; // 3 lines of sseg_ in sseg_'s local coords
double              mline_[3][4]; // 3 lines of mseg_in sseg_'s local coords
bool                haveline_;    // flag indicating whether line info has been build

map<int,MRTR::Point*> p_;
map<int,MRTR::Edge*>  e_;

}; // class Overlap

} // namespace MRTR
#endif // MRTR_OVERLAP_H
#endif // TRILINOS_PACKAGE
