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

#ifndef MRTR_PROJECTOR_H
#define MRTR_PROJECTOR_H

#include <ctime>
#include <iostream>
#include <iomanip>


using namespace std;

// ----------   User Defined Includes   ----------

namespace MRTR
{

// forward declarations
class MRTR::Interface;
class MRTR::Segment;
class MRTR::Node;

class  Projector 
{
public:
  
  // ctor
  explicit Projector(bool twoD);
  
  // dtor
  virtual ~Projector();
  
  // return true if this is a 2D projector, otherwise false
  bool IsTwoDimensional() { return twoD_; }

  // project a node onto a segment along the node's normal
  bool ProjectNodetoSegment_NodalNormal(MRTR::Node& node, MRTR::Segment& seg, double xi[]);
  
  // project a node onto a segment along the interpolated segment's normal field
  bool ProjectNodetoSegment_SegmentNormal(MRTR::Node& node, MRTR::Segment& seg, double xi[]);

private:  
  // don't want = operator
  Projector operator = (const Projector& old);
  // don't want copy-ctor
  Projector(MRTR::Projector& old);

  // evaluate F and gradF functions for ProjectNodetoSegment_NodalNormal
  double evaluate_F_2D_NodalNormal(MRTR::Node& node, MRTR::Segment& seg, double eta);  
  double evaluate_gradF_2D_NodalNormal(MRTR::Node& node, MRTR::Segment& seg, double eta);

  // evaluate F and gradF functions for ProjectNodetoSegment_SegmentNormal
  double evaluate_F_2D_SegmentNormal(MRTR::Node& node, MRTR::Segment& seg, double eta);  
  double evaluate_gradF_2D_SegmentNormal(MRTR::Node& node, MRTR::Segment& seg, double eta);


private:

  bool twoD_;   // dimension of the projection, true if 2-dimensional

};

} // namespace MRTR

#endif // MRTR_PROJECTOR_H
#endif // TRILINOS_PACKAGE
