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

#ifndef MRTR_NODE_H
#define MRTR_NODE_H

#include <ctime>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>


using namespace std;

// ----------   User Defined Includes   ----------

namespace MRTR
{

// forward declarations
class MRTR::Interface;
class MRTR::Segment;
class MRTR::ProjectedNode;

class  Node 
{
public:
  
  // ctor
  explicit Node(int Id, const double* x, int ndof, const int* dof);
  
  // ctor (for internal use only) to be used together with Pack()/UnPack()
  explicit Node();
  
  // copy-ctor
  Node(MRTR::Node& old);
  
  // dtor
  virtual ~Node();

  // get id of this node
  int Id() const { return Id_; }
  
  // print node
  bool Print() const;
  
  // get coords of this node
  const double* X() const { return x_; }
  
  // get normal of this node
  const double* N() const { return n_; }
  
  // set normal to this node
  bool SetN(double* n) { for (int i=0; i<3; ++i) n_[i] = n[i]; return true; }
  
  // get number of degrees of freedom on this node
  int Ndof() const { return dof_.size(); }

  // get number of lagrange multipliers on this node
  int Nlmdof() const { return LMdof_.size(); }
  
  // get vector of dofs on this node
  const int* Dof() const { return &(dof_[0]); }
  
  // get vector of lagrange multipliers on this node
  const int* LMDof() const { return &(LMdof_[0]); }
  
  // get number of segments adjacent to me
  int Nseg() const { return seg_.size(); }
  
  // get view of segment ids adjacent to this node
  int* SegmentIds() { return &(seg_[0]); }
  
  // get view of segments adjacent to this node
  MRTR::Segment** Segments() { return &(segptr_[0]); }
    
  // add a segment id to my adjacency list (checking whether already present)
  bool AddSegment(int sid);
  
  // pack function for communication
  double* Pack(int* size);
  
  // unpack function for communication
  bool UnPack(double* pack);
  
  // build an averaged normal from neighboring segments from this node
  bool BuildAveragedNormal();
  
  // construct ptrs to redundant segments from my segment id list
  bool GetPtrstoSegments(MRTR::Interface& interface);
  
  // store a ptr to a projeted node
  bool SetProjectedNode(MRTR::ProjectedNode* pnode);
  
  // get view of projected nodes (if there is any, otherwise NULL)
  MRTR::ProjectedNode** GetProjectedNodes(int& length);
  // get a view of the first projected node in the vector of projected nodes
  MRTR::ProjectedNode*  GetProjectedNode();
  
  // set dof number for lagrange multiplier
  bool SetLagrangeMultiplierId(int LMId);
  
protected:  
  // don't want = operator
  Node operator = (const Node& old);

protected:

  int                    Id_;      // Id of this node
  double                 x_[3];    // cordinates of this node 
  double                 n_[3];    // a nodal outward normal

  vector<int>            dof_;     // dof ids on this node 
  vector<int>            LMdof_;   // dofs ids of LM if this is slave side

  vector<int>            seg_;     // vector of segment ids adjacent to me
  vector<MRTR::Segment*> segptr_;  // ptrs to segments length nseg_

  //MRTR::ProjectedNode*   pnode_;   // the projected node on some segment
  vector<MRTR::ProjectedNode*>   pnode_;   // the projected node on some segment

};

} // namespace MRTR

// << operator
std::ostream& operator << (ostream& os, const MRTR::Node& node);

#endif // MRTR_NODE_H
#endif // TRILINOS_PACKAGE
