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

#ifndef MRTR_EDGE_H
#define MRTR_EDGE_H

#include <ctime>
#include <iostream>
#include <iomanip>


using namespace std;

// ----------   User Defined Includes   ----------

namespace MRTR
{

class Edge
{
public:

  // ctor
  Edge(int id, int node0, int node1);
  
  // dtor
  virtual ~Edge();

  // get edge id
  int Id() { return id_; }

  // set node ids
  bool SetNodes(int node0, int node1) 
  { nodeid_[0] = node0;
    nodeid_[1] = node1; return true;}
    
  // get view from nodal ids
  const int* GetNodesView() { return &nodeid_[0]; }

private:

  // don't want = operator
  Edge operator = (const Edge& old);
  // don't want copy-ctor
  Edge(MRTR::Edge& old);

private:

  int   id_;
  int   nodeid_[2];

}; // class Edge


} // namespace MRTR
#endif // MRTR_EDGE_H
#endif // TRILINOS_PACKAGE
