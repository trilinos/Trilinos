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

#include "mrtr_segment.H"
#include "mrtr_interface.H"
#include "mrtr_utils.H"

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 06/05|
 |  id               (in)  a unique segment id                          |
 |  nnode            (in)  number of nodes on this segment              |
 |  nodeId           (in)  unique node ids of nodes on this segment     |
 |                         nodeIds have to be sorted on input such that |
 |                         1D case: end nodes of segments first going   |
 |                                  through segment in mathematical     |
 |                                  positive sense. This is             |
 |                                  important to compute the direction  |
 |                                  of the outward normal of the segment|
 |                         2D case: corner nodes of segment first in    |
 |                                  counterclockwise order              |
 |                                  This is                             |
 |                                  important to compute the direction  |
 |                                  of the outward normal of the segment|
 *----------------------------------------------------------------------*/
MRTR::Segment::Segment(int id, int nnode, int* nodeId)
{
  Id_ = id;
  stype_ = MRTR::Segment::seg_none;
  
  nodeId_.resize(nnode);
  for (int i=0; i<nnode; ++i) nodeId_[i] = nodeId[i];
  nodeptr_.resize(0);
}

/*----------------------------------------------------------------------*
 | base class constructor                                    mwgee 07/05|
 *----------------------------------------------------------------------*/
MRTR::Segment::Segment() :
Id_(-1),
stype_(MRTR::Segment::seg_none)
{
  nodeId_.resize(0);
  nodeptr_.resize(0); 
}

/*----------------------------------------------------------------------*
 | base class copy ctor                                      mwgee 07/05|
 *----------------------------------------------------------------------*/
MRTR::Segment::Segment(MRTR::Segment& old)
{ 
  Id_      = old.Id_;
  stype_   = old.stype_; 
  nodeId_  = old.nodeId_;
  nodeptr_ = old.nodeptr_; 
  
  // copy the functions
  // this is not a deep copy but we simply copy the refcountptr
  map< int,RefCountPtr<MRTR::Function> >::iterator curr;
  for (curr = old.functions_.begin(); curr != old.functions_.end(); ++curr)
  {
    if (curr->second == null)
    {
      cout << "***ERR*** MRTR::Segment::BaseClone(MRTR::Segment& old):\n"
           << "***ERR*** function id " << curr->first << " is null\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);     
    }
    RefCountPtr<MRTR::Function> newfunc = curr->second;
    functions_.insert(pair< int,RefCountPtr<MRTR::Function> >(curr->first,newfunc));
  }
  
}

/*----------------------------------------------------------------------*
 | base class destructor                                     mwgee 07/05|
 *----------------------------------------------------------------------*/
MRTR::Segment::~Segment()
{
  nodeId_.clear();
  nodeptr_.clear();
  functions_.clear();
}

/*----------------------------------------------------------------------*
 |  print segment                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
bool MRTR::Segment::Print() const
{ 
  cout << "Segment " << setw(6) << Id_; 
  if (stype_ == MRTR::Segment::seg_Linear1D)
    cout << " Typ Linear1D   ";
  if (stype_ == MRTR::Segment::seg_Quadratic1D)
    cout << " Typ Quadratic1D";
  if (stype_ == MRTR::Segment::seg_BiLinearQuad)
    cout << " Typ BiLinearQuad ";
  if (stype_ == MRTR::Segment::seg_BiLinearTri)
    cout << " Typ BiLinearTri";
  if (stype_ == MRTR::Segment::seg_none)
    cout << " Typ NONE       ";
  cout << " #Nodes " << nodeId_.size() << " Nodes: ";
  for (int i=0; i<nodeId_.size(); ++i)
    cout << setw(6) << nodeId_[i] << "  ";
  cout << "  #Functions " << functions_.size() << "  Types: ";
  map<int,RefCountPtr<MRTR::Function> >::const_iterator curr;
  for (curr=functions_.begin(); curr != functions_.end(); ++curr)
    cout << curr->second->Type() << "  ";
  cout << endl;
  return true;
}

/*----------------------------------------------------------------------*
 | attach a certain shape function to this segment           mwgee 06/05|
 | the user is not supposed to destroy func!                            |
 | the user can set func to several segments!                           |
 *----------------------------------------------------------------------*/
bool MRTR::Segment::SetFunction(int id, MRTR::Function* func)
{ 
  if (id<0)
  {
    cout << "***ERR*** MRTR::Segment::SetFunction:\n"
         << "***ERR*** id = " << id << " < 0 (out of range)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!func)
  {
    cout << "***ERR*** MRTR::Segment::SetFunction:\n"
         << "***ERR*** func = NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  // check for existing function with this id and evtentually overwrite
  map<int,RefCountPtr<MRTR::Function> >::iterator curr = functions_.find(id);
  if (curr != functions_.end())
  {
    curr->second = rcp(func);
    return true;
  }
  RefCountPtr<MRTR::Function> newfunc = rcp(func);
  functions_.insert(pair<int,RefCountPtr<MRTR::Function> >(id,newfunc));
  return true;
}

/*----------------------------------------------------------------------*
 | attach a certain shape function to this segment           mwgee 11/05|
 | the user is not supposed to destroy func!                            |
 | the user can set func to several segments!                           |
 *----------------------------------------------------------------------*/
bool MRTR::Segment::SetFunction(int id, RefCountPtr<MRTR::Function> func)
{ 
  if (id<0)
  {
    cout << "***ERR*** MRTR::Segment::SetFunction:\n"
         << "***ERR*** id = " << id << " < 0 (out of range)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (func==null)
  {
    cout << "***ERR*** MRTR::Segment::SetFunction:\n"
         << "***ERR*** func = NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  // check for existing function with this id and evtentually overwrite
  map<int,RefCountPtr<MRTR::Function> >::iterator curr = functions_.find(id);
  if (curr != functions_.end())
  {
    curr->second = func;
    return true;
  }
  RefCountPtr<MRTR::Function> newfunc = func;
  functions_.insert(pair<int,RefCountPtr<MRTR::Function> >(id,newfunc));
  return true;
}

/*----------------------------------------------------------------------*
 | evaluate the shape function id at the point xi            mwgee 06/05|
 | id     (in)   id of the function to evaluate                         |
 | xi     (in)   natural coordinates -1<xi<1 where to eval the function |
 | val    (out)  function values, if NULL on input, no evaluation       |
 | valdim (in)   dimension of val                                       |
 | deriv  (out)  derivatives of functions at xi, if NULL on input,      |
 |               no evaluation                                          | 
 *----------------------------------------------------------------------*/
bool MRTR::Segment::EvaluateFunction(int id, const double* xi, double* val, 
                                     int valdim, double* deriv)
{ 
  map<int,RefCountPtr<MRTR::Function> >::iterator curr = functions_.find(id);
  if (curr == functions_.end())
  {
    cout << "***ERR*** MRTR::Segment::EvaluateFunction:\n"
         << "***ERR*** function id " << id << " does not exist on this segment\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);     
  }
  curr->second->EvaluateFunction(xi,val,valdim,deriv);
  
  return true;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 06/05|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const MRTR::Segment& seg)
{
  seg.Print(); 
  return os;
}

/*----------------------------------------------------------------------*
 | get local numbering id for global node id on this segment mwgee 07/05|
 | return -1 of nid is not adjacent to this segment                     |
 *----------------------------------------------------------------------*/
int MRTR::Segment::GetLocalNodeId(int nid)
{ 
  int lid=-1;
  for (int i=0; i<Nnode(); ++i)
    if (nodeId_[i]==nid)
    {
      lid = i;
      break;
    }
  if (lid<0)
  {
    cout << "***ERR*** MRTR::Segment::GetLocalNodeId:\n"
         << "***ERR*** cannot find node " << nid << " in segment " << this->Id() << " list of nodes\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);     
  }
  return lid;
}

/*----------------------------------------------------------------------*
 | build an outward normal at a node adjacent to this        mwgee 07/05|
 *----------------------------------------------------------------------*/
double* MRTR::Segment::BuildNormalAtNode(int nid)
{ 
  // find this node in my list of nodes and get local numbering for it
  int lid = GetLocalNodeId(nid);
  
  // depending on what type of segment I am get local coordinates
  double xi[2];
  LocalCoordinatesOfNode(lid,xi);
  
  // build an outward unit normal at xi
  return BuildNormal(xi);
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 | construct ptrs to redundant nodes from my node id list               |
 *----------------------------------------------------------------------*/
bool MRTR::Segment::GetPtrstoNodes(MRTR::Interface& interface)
{ 
  if (!interface.IsComplete()) return false;
  if (!interface.lComm()) return true;
  if (!nodeId_.size()) return false;
  
  // vector nodeptr_ might already exist, recreate it
  nodeptr_.clear();
  nodeptr_.resize(nodeId_.size());
  
  for (int i=0; i<nodeId_.size(); ++i)
  {
    nodeptr_[i] = interface.GetNodeView(nodeId_[i]).get();
    if (!nodeptr_[i])
    {
      cout << "***ERR*** MRTR::Segment::GetPtrstoNodes:\n"
           << "***ERR*** interface " << interface.Id() << " GetNodeView failed\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);     
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 10/05|
 | construct ptrs to nodes from vector                                  |
 *----------------------------------------------------------------------*/
bool MRTR::Segment::GetPtrstoNodes(vector<MRTR::Node*>& nodes)
{ 
  if (!nodeId_.size()) return false;
  
  // vector nodeptr_ might already exist, recreate it
  nodeptr_.clear();
  nodeptr_.resize(nodeId_.size());
  
  for (int i=0; i<nodeId_.size(); ++i)
  {
    bool foundit = true;
    for (int j=0; j<nodes.size(); ++j)
      if (nodes[j]->Id() == nodeId_[i])
      {
        foundit = true;
        nodeptr_[i] = nodes[j];
        break;
      }
    if (!foundit)
    {
      cout << "***ERR*** MRTR::Segment::GetPtrstoNodes:\n"
           << "***ERR*** cannot find node " << nodeId_[i] << " in vector\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);     
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 | return type of function with certain id                              |
 *----------------------------------------------------------------------*/
MRTR::Function::FunctionType MRTR::Segment::FunctionType(int id)
{ 
  // find the function with id id
  map<int,RefCountPtr<MRTR::Function> >::iterator curr = functions_.find(id);
  if (curr==functions_.end()) 
    return MRTR::Function::func_none;
  else
    return curr->second->Type();
}

#if 0
/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 | get ptr to function with id id                                       |
 *----------------------------------------------------------------------*/
MRTR::Function* MRTR::Segment::GetFunction(int id)
{ 
  // find the function with id id
  map<int,RefCountPtr<MRTR::Function> >::iterator curr = functions_.find(id);
  if (curr==functions_.end()) 
    return NULL;
  else
    return curr->second;
}
#endif

#endif // TRILINOS_PACKAGE
