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

#ifndef MOERTEL_SEGMENT_BILINEARTRI_H
#define MOERTEL_SEGMENT_BILINEARTRI_H

#include <ctime>
#include <iostream>
#include <map>

#include "mrtr_segment.H"
#include "mrtr_functions.H"
#include "mrtr_utils.H"

using namespace std;

namespace MOERTEL
{

class  Segment_BiLinearTri : public MOERTEL::Segment
{
public:
  
  // ctor
  explicit Segment_BiLinearTri(int id, int nnode, int* nodeId, int out);
  
  // ctor (for internal use only) to be used together with Pack()/UnPack()
  explicit Segment_BiLinearTri(int out);
  
  // copy-ctor
  Segment_BiLinearTri(MOERTEL::Segment_BiLinearTri& old);
  
  // dtor
  virtual ~Segment_BiLinearTri();

  // clone
  MOERTEL::Segment* Clone();

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
  Segment_BiLinearTri operator = (const Segment& old);

};

} // namespace MOERTEL

// << operator
ostream& operator << (ostream& os, const MOERTEL::Segment_BiLinearTri& seg);


#endif // MOERTEL_SEGMENT_BILINEARTRI_H
#endif // TRILINOS_PACKAGE
