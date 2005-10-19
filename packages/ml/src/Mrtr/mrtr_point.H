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

#ifndef MRTR_POINT_H
#define MRTR_POINT_H

#include <ctime>
#include <iostream>
#include <iomanip>


using namespace std;

// ----------   User Defined Includes   ----------

namespace MRTR
{

class Point
{
public:

  // ctor
  Point(int id, double* xi);
  
  // dtor
  virtual ~Point();
  
  // get point id
  int Id() { return id_; }
  
  // return view of coords
  const double* Xi() { return &xi_[0]; }

  // set Xi coords
  bool SetXi(double* xi) { xi_[0] = xi[0]; xi_[1] = xi[1]; return true; }

private:

  // don't want = operator
  Point operator = (const Point& old);
  // don't want copy-ctor
  Point(MRTR::Point& old);

private:

  int     id_;
  double  xi_[2];


}; // class Point


} // namespace MRTR
#endif // MRTR_POINT_H
#endif // TRILINOS_PACKAGE
