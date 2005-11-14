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

#include <Epetra_SerialComm.h>
#include "Teuchos_RefCountPtr.hpp"


#include "mrtr_point.H"


using namespace std;
using namespace Teuchos;

// ----------   User Defined Includes   ----------

namespace MRTR
{

// forward declarations
class MRTR::Interface;
class MRTR::Segment;
class MRTR::Node;
class MRTR::Integrator;

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
    poly_0_6,
    poly_0_0,
    poly_0_4_same,
    poly_0_4_dif,
    poly_0_2
  };

  // the interface class is a friend to this class
  friend class Interface;
  
  // the integrator class is a friend to this class
  friend class Integrator;
  
  // ctor
  explicit Overlap(MRTR::Segment& sseg, MRTR::Segment& mseg, MRTR::Interface& inter);
  
  // dtor
  virtual ~Overlap();
  
  // compute overlap between 2 segments and it's polygon
  bool ComputeOverlap();
    
private:  
  // don't want = operator
  Overlap operator = (const Overlap& old);
  // don't want copy-ctor
  Overlap(MRTR::Overlap& old);

private:

  // build line information from triangles in slave coord system
  bool build_lines_s();
  // build line information from triangles in master coord system
  bool build_lines_m();
  // build projection of master nodes onto slave segment
  bool build_mxi();
  // build projection of slave nodes onto master segment
  bool build_sxim();
  // build projection of slave nodes onto master segment
  bool build_sxi();
  // build the outward normal to the sseg
  bool build_normal();
  
  //------------------------------------------ perform quick coarse search
  bool QuickOverlapTest();

  //-------------------------------------------- perform clipping algorithm
  bool Clipelements();
  // test whether a point is inside or outside of a clip edge
  bool Clip_TestPoint(const double* N, const double* PE, const double* P, double eps);
  // find intersection of clipping edge with line
  bool Clip_Intersect(const double* N,const double* PE,const double* P0,const double* P1,double* xi);
  // find parameterization alpha for point on line
  double Clip_ParameterPointOnLine(const double* P0,const double* P1,const double* P);

  // -------------------------------------finding a convex hull for the polygon
  // create a convex hull of a set of given points in 2D
  bool ConvexHull(map<int,RefCountPtr<MRTR::Point> >& p);
  bool MakeRightTurnUpper(int i,map<int,RefCountPtr<MRTR::Point> >& hull);
  bool MakeRightTurnLower(int i,map<int,RefCountPtr<MRTR::Point> >& hull);
  void RemovePointBefore(int i,map<int,RefCountPtr<MRTR::Point> >& hull);
  void RemovePointAfter(int i,map<int,RefCountPtr<MRTR::Point> >& hull);

  //-------------------------------------------make triangulization of polygon
  bool Triangulization();
  
  //-----------------------------------------------------------collapse points
  bool CollapsePoints(map<int,RefCountPtr<MRTR::Point> >& p, const double eps);
  
protected:
  //----------------------------------- tools to handle the triangulation map
  // add a segment to the triangulization
  bool AddSegment(int id, MRTR::Segment* seg);
  // return # segments the overlap was discretized with
  int Nseg() { return s_.size(); }
  // get view of segments in triangulation map
  void SegmentView(vector<RefCountPtr<MRTR::Segment> >& segs);

  //------------------------------------------ tools to handle the polygon map
  // add a point to the polygon
  bool AddPointtoPolygon(const int id, const double* P);
  bool AddPointtoPolygon(map<int,RefCountPtr<MRTR::Point> >& p, const int id, const double* P);
  // remove a point from the polygon
  bool RemovePointfromPolygon(const int id,const double* P);
  // get the size of the nodal polygon
  int SizePointPolygon() { return p_.size(); }
  // get view of point in polygon, calling routine is responsible for freeing
  void PointView(vector<RefCountPtr<MRTR::Point> >& points);
  void PointView(map<int,RefCountPtr<MRTR::Point> >& p, vector<RefCountPtr<MRTR::Point> >& points);
  void PointView(vector<MRTR::Point*>& p,const int* nodeids,const int np);
  // copy a point polygon to another polygon
  bool CopyPointPolygon(map<int,RefCountPtr<MRTR::Point> >& from, map<int,RefCountPtr<MRTR::Point> >& to);

private:

MRTR::Interface&                     inter_;
MRTR::Segment&                       sseg_;
MRTR::Segment&                       mseg_;

bool                                 overlap_;   // true if segments overlap
int                                  nnode_;     // # of master nodes

//------------------ data in sseg_'s local coord system
double                               mxi_[3][2]; // local coords of mnodes projected onto sseg
bool                                 havemxi_;   // flag indicating whether mxi_ and min_ have been build

double                               sxi_[3][2]; // local coords of snodes in sseg's coord system
double                               n_[3][2];   // outward normal to segment's edges in local coords
bool                                 havesxi_;   // flag indicating whether sxi_ and sin_ have been build

double                               sline_[3][4]; // 3 lines of sseg_ in sseg_'s local coords
double                               mline_[3][4]; // 3 lines of mseg_in sseg_'s local coords
bool                                 havelines_;   // flag indicating whether line info has been build

//------------------ data in mseg_'s local coord system
double                               sxim_[3][2];  // coords of slave nodes on master segment
bool                                 havesxim_;    // flag indicating whether sxim_ was built

double                               slinem_[3][4]; // 3 lines of sseg_ in mseg_'s local coords
double                               mlinem_[3][4]; // 3 lines of mseg_in mseg_'s local coords
bool                                 havelinem_;   // flag indicating built of line information in master coords


map<int,RefCountPtr<MRTR::Point> >   p_;       // map holding polygon points and later points of triangulation
map<int,RefCountPtr<MRTR::Segment> > s_;       // map holding segments of polygon triangulation


}; // class Overlap

} // namespace MRTR
#endif // MRTR_OVERLAP_H
#endif // TRILINOS_PACKAGE
