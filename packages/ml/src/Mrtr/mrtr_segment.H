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

#ifndef MRTR_SEGMENT_H
#define MRTR_SEGMENT_H

#include <ctime>
#include <iostream>
#include <map>
#include <vector>

#include "Teuchos_RefCountPtr.hpp"
#include "mrtr_function.H"

using namespace std;
using namespace Teuchos;

namespace MRTR
{

class MRTR::Interface;
class MRTR::Node;

class Segment 
{
public:
  
  // types of segments
  enum SegmentType 
  {
    seg_none,
    seg_Linear1D,
    seg_Quadratic1D,
    seg_BiLinearQuad,
    seg_BiLinearTri
  };

  //------------- methods used from base class

  // ctor
  Segment(int id, int nnode, int* nodeId);
  
  // ctor
  Segment();
  
  // copy-ctor
  Segment(MRTR::Segment& old);
  
  // dtor
  virtual ~Segment();
  
  // get id of this segment
  int Id() const { return Id_; }
  
  // get number of nodes on this segment
  int Nnode() const { return nodeId_.size(); }
  
  // get type of segment
  MRTR::Segment::SegmentType Type() const { return stype_; }
  
  // get view of node Ids on this segment
  const int* NodeIds() const { return &(nodeId_[0]); }
  
  // return ptr vector to my nodes (might be NULL before call to Interface::Complete())
  MRTR::Node** Nodes() { return &(nodeptr_[0]); }
  
  // return number of functions set
  int Nfunctions() { return functions_.size(); }
  
  // return function type for a certain function id
  MRTR::Function::FunctionType FunctionType(int id);
  
  // attach a function to this segment
  bool SetFunction(int id, MRTR::Function* func);
  bool SetFunction(int id, RefCountPtr<MRTR::Function> func);

#if 0  
  // get ptr to function with id id
  MRTR::Function* GetFunction(int id);
#endif
  
  // evaluate shape function with certain id;
  bool EvaluateFunction(int id, const double* xi, double* val, int valdim, double* deriv);
  
  // build an outward normal at a node adjacent to this
  double* BuildNormalAtNode(int nid);
  
  // construct ptrs to redundant nodes from my node id list
  bool GetPtrstoNodes(MRTR::Interface& interface);
  bool GetPtrstoNodes(vector<MRTR::Node*>& nodes);

  // printout Segment
  virtual bool Print() const;
  
  // get local numbering id for global node id adjacent to segment
  int GetLocalNodeId(int nid);


  //------------- pure virtual methods to be implemented in derived class

  // clone
  virtual MRTR::Segment* Clone() = 0;

  // pack function for communication
  virtual int* Pack(int* size) = 0;
  
  // unpack function for communication
  virtual bool UnPack(int* pack) = 0;
  
  // build an outward normal
  virtual double* BuildNormal(double* xi) = 0;
  
  // integrate the are of this segment
  virtual double Area() = 0;

  // build the basis vector and metrix tensor at a given point in element
  virtual double Metric(double* xi, double g[], double G[][3]) = 0;

  // get local coordinates at a node in LOCAL numbering
  virtual bool LocalCoordinatesOfNode(int lid, double* xi) = 0;

protected:

  int                      Id_;         // this segments unique id
  vector<int>              nodeId_;     // vector of unique node ids 
  vector<MRTR::Node*>      nodeptr_;    // vector with ptrs to nodes adj to me

  SegmentType              stype_;      // segment type
  map< int,RefCountPtr<MRTR::Function> > functions_;  // functions that live on this segment

};

} // namespace MRTR

// << operator
ostream& operator << (ostream& os, const MRTR::Segment& seg);

#endif // MRTR_SEGMENT_H
#endif // TRILINOS_PACKAGE
