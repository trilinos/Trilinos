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

#ifndef MOERTEL_POINT_H
#define MOERTEL_POINT_H

#include <ctime>
#include <iostream>
#include <iomanip>

#include "mrtr_node.H"

using namespace std;

// ----------   User Defined Includes   ----------

namespace MOERTEL
{


class Point
{
public:

  // ctor
  Point(const int id, const double* xi, int out);
  
  // dtor
  virtual ~Point();
  
  // return level of output
  int OutLevel() { return outputlevel_; }
  
  // print
  void Print() const;

  // get point id
  inline int Id() { return id_; }
  
  // return view of local coords
  inline const double* Xi() { return &xi_[0]; }

  // return view of global coords
  inline const double* X() { if (node_ != null) return node_->X(); else return NULL; }

  // return view of node
  inline RefCountPtr<MOERTEL::Node> Node() { return node_; }
  
  // set Xi coords
  bool SetXi(const double* xi) { xi_[0] = xi[0]; xi_[1] = xi[1]; return true; }

  // set x coords
  bool SetNode(MOERTEL::Node* node) { node_ = rcp(node); return true; }

  // store function values in the vals_ vector
  void StoreFunctionValues(int place, double* val, int valdim);

  // return ptr to vals_ vector
  vector<double>* FunctionValues() { return vals_; }

private:

  // don't want = operator
  Point operator = (const Point& old);
  // don't want copy-ctor
  Point(MOERTEL::Point& old);

private:

  int                      id_;      // id of this point
  int                      outputlevel_; 
  double                   xi_[2];   // local coords in some slave elements coord system
  RefCountPtr<MOERTEL::Node>  node_;    // a node at this point (contains real world coords)
  vector<double>           vals_[3]; // [0] values of shape function 0 from sseg
                                     // [1] values of shape function 1 from sseg
                                     // [2] values of shape function 0 from mseg
}; // class Point


} // namespace MOERTEL

// << operator
std::ostream& operator << (ostream& os, const MOERTEL::Point& point);


#endif // MOERTEL_POINT_H
#endif // TRILINOS_PACKAGE
