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


using namespace std;

// ----------   User Defined Includes   ----------

namespace MRTR
{

// forward declarations
class MRTR::Interface;
class MRTR::Segment;
class MRTR::Node;
class MRTR::Point;
class MRTR::Edge;

class  Overlap 
{
public:
  
  // ctor
  explicit Overlap(MRTR::Segment& sseg, MRTR::Segment& mseg, MRTR::Interface& inter);
  
  // dtor
  virtual ~Overlap();
  
  // determine whether overlap between 2 segments exists
  bool HaveOverlap();
  
  // compute overlap between 2 segments
  bool ComputeOverlap();
  
private:  
  // don't want = operator
  Overlap operator = (const Overlap& old);
  // don't want copy-ctor
  Overlap(MRTR::Overlap& old);

  // compute intersection between 2 lines in local coordinates
  bool intersect_line2D(int s, int m, double* alpha, double* beta, double* xi);

  // build line information from triangles
  bool build_lines();
  // build projection of master nodes onto slave segment
  bool build_mxi();
  // build projection of slave nodes onto master segment
  bool build_sxi();

  // add point to the polygon
  bool AddPoint(MRTR::Point* p);
  // add edge to the polygon
  bool AddEdge(MRTR::Edge* e);
  // get view point from polygon
  MRTR::Point* PointView(int id);
  // get view from edge from the polygon
  MRTR::Edge* EdgeView(int id);
  
private:

MRTR::Interface&    inter_;
MRTR::Segment&      sseg_;
MRTR::Segment&      mseg_;

bool                overlap_;   // true if segments overlap
int                 nnode_;     // # of master nodes

double              mxi_[3][2]; // local coords of mnodes projected onto sseg
bool                min_[3];    // true if projection of mnodes is inside sseg
bool                havemxi_;   // flag indicating whether mxi_ and min_ have been build

double              sxi_[3][2]; // local coords of snodes projected onto mseg
bool                sin_[3];    // true if projection of snodes is inside mseg
bool                havesxi_;   // flag indicating whether sxi_ and sin_ have been build

double              sline_[3][4]; // 3 lines of sseg_ in sseg_'s local coords
double              mline_[3][4]; // 3 lines of mseg_in sseg_'s local coords
bool                haveline_;    // flag indicating whether line info has been build

map<int,MRTR::Point*> p_;         // map of point on the polygon
map<int,MRTR::Edge*>  e_;         // map of edges on the polygon

}; // class Overlap

} // namespace MRTR
#endif // MRTR_OVERLAP_H
#endif // TRILINOS_PACKAGE
