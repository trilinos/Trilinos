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

#ifndef MRTR_PROJECTNODE_H
#define MRTR_PROJECTNODE_H

#include <ctime>
#include <iostream>
#include <iomanip>


#include "mrtr_node.H"

using namespace std;

// ----------   User Defined Includes   ----------

namespace MRTR
{

// forward declarations
class MRTR::Interface;
class MRTR::Segment;

class  ProjectedNode : public MRTR::Node
{
public:
  
  // ctor 
  explicit ProjectedNode(MRTR::Node& basenode, const double* xi, MRTR::Segment* pseg);
  
  // ctor for the case of orthogonal projection 
  explicit ProjectedNode(MRTR::Node& basenode, const double* xi, MRTR::Segment* pseg, int orthseg);
  
  // copy-ctor
  ProjectedNode(MRTR::ProjectedNode& old);
  
  // dtor
  virtual ~ProjectedNode();

  // print node
  bool Print() const;
  
  // get view of local coordinates xi
  const double* Xi() { return xi_; }
  
  // get view of segment I am on
  MRTR::Segment* Segment() { return pseg_; }
  
  // get id of segment this projection is orthogonal to
  int OrthoSegment() { return orthseg_; }
  
protected:  
  // don't want = operator
  ProjectedNode operator = (const ProjectedNode& old);

protected:

  double xi_[2];            // local coordinates of this projected node on pseg_;
  MRTR::Segment* pseg_;     // segment this projected node is on
  int orthseg_;             // id of segment this projection is orthogonal to
                            // (used only in orth. projection, otherwise -1)
};

} // namespace MRTR

// << operator
std::ostream& operator << (ostream& os, const MRTR::ProjectedNode& pnode);

#endif // MRTR_PROJECTNODE_H
#endif // TRILINOS_PACKAGE
