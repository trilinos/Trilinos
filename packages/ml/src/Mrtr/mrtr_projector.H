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

  //====2D/3D projection methods
  
  // project a node onto a segment along the node's normal
  bool ProjectNodetoSegment_NodalNormal(const MRTR::Node& node, MRTR::Segment& seg, double xi[]);
  // project a node onto a segment along the interpolated segment's normal field
  bool ProjectNodetoSegment_SegmentNormal(const MRTR::Node& node, MRTR::Segment& seg, double xi[]);

  //====2D projection methods
  
  // project a node onto a segment orthogonal to the segment
  bool ProjectNodetoSegment_SegmentOrthogonal(const MRTR::Node& node, MRTR::Segment& seg, double xi[]);
  bool ProjectNodetoSegment_Orthogonal_to_Slave(const MRTR::Node& snode, MRTR::Segment& seg, double xi[], MRTR::Segment& sseg);

private:  
  // don't want = operator
  Projector operator = (const Projector& old);
  // don't want copy-ctor
  Projector(MRTR::Projector& old);

  //====2D projection methods
  
  // evaluate F and gradF functions for ProjectNodetoSegment_NodalNormal in 2D
  double evaluate_F_2D_NodalNormal(const MRTR::Node& node, MRTR::Segment& seg, double eta);  
  double evaluate_gradF_2D_NodalNormal(const MRTR::Node& node, MRTR::Segment& seg, double eta);

  // evaluate F and gradF functions for ProjectNodetoSegment_SegmentNormal in 2D
  double evaluate_F_2D_SegmentNormal(const MRTR::Node& node, MRTR::Segment& seg, double eta);  
  double evaluate_gradF_2D_SegmentNormal(const MRTR::Node& node, MRTR::Segment& seg, double eta);

  // evaluate F and gradF functions for ProjectNodetoSegment_SegmentOrthogonal in 2D
  double evaluate_F_2D_SegmentOrthogonal(const MRTR::Node& node, MRTR::Segment& seg, double eta);  
  double evaluate_gradF_2D_SegmentOrthogonal(const MRTR::Node& node, MRTR::Segment& seg, double eta);
  
  // evalauate F and gradF functions for ProjectNodetoSegment_Orthogonal_to_Slave in 2D
  double evaluate_F_2D_SegmentOrthogonal_to_g(const MRTR::Node& node, MRTR::Segment& seg, double eta, double* g);  
  double evaluate_gradF_2D_SegmentOrthogonal_to_g(const MRTR::Node& node, MRTR::Segment& seg, double eta, double* g);
  
  //====3D projection methods

  // evaluate F and gradF functions for ProjectNodetoSegment_NodalNormal in 3D
  bool evaluate_FgradF_3D_NodalNormal(double* F,double dF[][3],const MRTR::Node& node, 
                                      MRTR::Segment& seg,double* eta,double alpha);
  // evaluate F and gradF functions for ProjectNodetoSegment_SegmentNormal in 3D
  bool evaluate_FgradF_3D_SegmentNormal(double* F,double dF[][3],const MRTR::Node& node, 
                                        MRTR::Segment& seg,double* eta,double alpha);

private:

  bool twoD_;   // dimension of the projection, true if 2-dimensional

};

} // namespace MRTR

#endif // MRTR_PROJECTOR_H
#endif // TRILINOS_PACKAGE
