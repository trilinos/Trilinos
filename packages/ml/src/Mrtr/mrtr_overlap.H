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
  
  // ctor
  explicit Overlap(MRTR::Segment& sseg, MRTR::Segment& mseg, MRTR::Interface& inter);
  
  // determine overlap between 2 segments
  bool HaveOverlap();
  
  // dtor
  virtual ~Overlap();
  
private:  
  // don't want = operator
  Overlap operator = (const Overlap& old);
  // don't want copy-ctor
  Overlap(MRTR::Overlap& old);

private:

MRTR::Interface&    inter_;
MRTR::Segment&      sseg_;
MRTR::Segment&      mseg_;

bool                overlap_;  // true if segments overlap
int                 nnode_;    // # of master nodes
double              xi_[4][2]; // local coords of mnodes projected onto sseg
bool                in_[4];    // true if projection of mnodes is inside sseg

};

} // namespace MRTR

#endif // MRTR_OVERLAP_H
#endif // TRILINOS_PACKAGE
