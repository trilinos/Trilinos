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

#ifndef MOERTEL_NODE_H
#define MOERTEL_NODE_H

#include <ctime>
#include <iostream>
#include <iomanip>
#include <vector>
#include <map>

#include "Teuchos_RefCountPtr.hpp"

using namespace std;
using namespace Teuchos;

// ----------   User Defined Includes   ----------

namespace MOERTEL
{

// forward declarations
class MOERTEL::Interface;
class MOERTEL::Segment;
class MOERTEL::ProjectedNode;

class  Node 
{
public:
  
  // ctor
  explicit Node(int Id, const double* x, int ndof, const int* dof, bool isonboundary, int out);
  
  // ctor (for internal use only) to be used together with Pack()/UnPack()
  explicit Node(int out);
  
  // copy-ctor
  Node(const MOERTEL::Node& old);
  
  // dtor
  virtual ~Node();

  // return level of output
  inline int OutLevel() { return outputlevel_; }
  
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
  inline int Ndof() const { return dof_.size(); }

  // get number of lagrange multipliers on this node
  inline int Nlmdof() const { return LMdof_.size(); }
  
  // get vector of dofs on this node
  inline const int* Dof() const { return &(dof_[0]); }
  
  // get vector of lagrange multipliers on this node
  inline const int* LMDof() const { return &(LMdof_[0]); }
  
  // get number of segments adjacent to me
  inline int Nseg() const { return seg_.size(); }
  
  // get view of segment ids adjacent to this node
  inline int* SegmentIds() { return &(seg_[0]); }
  
  // get view of segments adjacent to this node
  inline MOERTEL::Segment** Segments() { return &(segptr_[0]); }
    
  // add a segment id to my adjacency list (checking whether already present)
  bool AddSegment(int sid);
  
  // pack function for communication
  double* Pack(int* size);
  
  // unpack function for communication
  bool UnPack(double* pack);
  
  // build an averaged normal from neighboring segments from this node
  bool BuildAveragedNormal();
  
  // construct ptrs to redundant segments from my segment id list
  bool GetPtrstoSegments(MOERTEL::Interface& interface);
  
  // store a ptr to a projeted node
  bool SetProjectedNode(MOERTEL::ProjectedNode* pnode);
  
  // get view of projected nodes (if there is any, otherwise NULL)
  RefCountPtr<MOERTEL::ProjectedNode>* GetProjectedNode(int& length);
  // get a view of the first projected node in the vector of projected nodes
  RefCountPtr<MOERTEL::ProjectedNode>  GetProjectedNode();
  
  // set dof number for lagrange multiplier
  bool SetLagrangeMultiplierId(int LMId);
  
  // add a value to the Drow_ map
  void AddDValue(double val, int col);
  // add a value to the Mrow_ map
  void AddMValue(double val, int col);
  // get M or D map
  inline RefCountPtr<map<int,double> > GetD() { return Drow_; }
  inline RefCountPtr<map<int,double> > GetM() { return Mrow_; }

  // set/get the iscorner_ flag
  inline void SetCorner() { iscorner_ = true; return; }
  inline bool IsCorner() const { return iscorner_; }

  // get the isonboundary_ flag
  inline bool IsOnBoundary() const { return isonboundary_; }
  
  // get the number of internal nodes supporting this boundary node
  inline int NSupportSet() const { return supportedby_.size(); }
  
  // set an internal node supporting this boundary node
  void AddSupportedByNode(MOERTEL::Node* suppnode) 
  { supportedby_[suppnode->Id()] = suppnode; return; }
  
protected:  
  // don't want = operator
  Node operator = (const Node& old);

protected:

  int                          Id_;       // Id of this node
  double                       x_[3];     // cordinates of this node 
  double                       n_[3];     // a nodal outward normal

  int                          outputlevel_;

  bool                         iscorner_; // flag indicating whether this is a node
                                          // shared by several interfaces in 2D 

  bool                         isonboundary_; // flag indicating whether this node is on the boundary
                                              // of a 3D interface
                                              // Info only matters for the slave side
  map<int,MOERTEL::Node*>      supportedby_;
  map<int,MOERTEL::Node*>      supporting_;

  vector<int>                  dof_;      // dof ids on this node 
  vector<int>                  LMdof_;    // dofs ids of LM if this is slave side

  vector<int>                  seg_;      // vector of segment ids adjacent to me
  vector<MOERTEL::Segment*>    segptr_;   // ptrs to segments length nseg_ (not in charge of destorying the ptrs))

  RefCountPtr<map<int,double> > Drow_;    // a nodal block row of D
  RefCountPtr<map<int,double> > Mrow_;    // a nodal block row of D

  vector< RefCountPtr<MOERTEL::ProjectedNode> > pnode_;    // the projected node on some segment

};

} // namespace MOERTEL

// << operator
std::ostream& operator << (ostream& os, const MOERTEL::Node& node);

#endif // MOERTEL_NODE_H
#endif // TRILINOS_PACKAGE
