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

#ifndef MRTR_SEGMENT_LINEAR1D_H
#define MRTR_SEGMENT_LINEAR1D_H

#include <ctime>
#include <iostream>
#include <map>

#include "mrtr_segment.H"
#include "mrtr_functions.H"
#include "mrtr_utils.H"

using namespace std;

namespace MRTR
{

class  Segment_Linear1D : public MRTR::Segment
{
public:
  
  // ctor
  explicit Segment_Linear1D(int id, int nnode, int* nodeId);
  
  // ctor (for internal use only) to be used together with Pack()/UnPack()
  explicit Segment_Linear1D();
  
  // copy-ctor
  Segment_Linear1D(MRTR::Segment_Linear1D& old);
  
  // dtor
  virtual ~Segment_Linear1D();

  // clone
  MRTR::Segment* Clone();

  // (derived) pack function for communication
  int* Pack(int* size);
  
  // (derived) unpack function for communication
  bool UnPack(int* pack);
  
  // (derived) build an outward normal at a node adjacent to this
  double* BuildNormal(double* xi);
  
  // (derived) integrate the are of this segment
  double Area();
  
  // (derived) build the basis vector and metrix tensor at a given point in element
  double Metric(double* xi, double g[], double G[][3]);

  // (derived) get local coordinates at a node in LOCAL numbering
  bool LocalCoordinatesOfNode(int lid, double* xi);
  

protected:  

  // don't want = operator
  Segment_Linear1D operator = (const Segment& old);

};

} // namespace MRTR

// << operator
ostream& operator << (ostream& os, const MRTR::Segment_Linear1D& seg);


#endif // MRTR_SEGMENT_LINEAR1D_H
#endif // TRILINOS_PACKAGE
