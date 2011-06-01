/*
#@HEADER
# ************************************************************************
#
#                          Moertel FE Package
#                 Copyright (2006) Sandia Corporation
#
# Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
# license for use of this work by or on behalf of the U.S. Government.
#
# This library is free software; you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 2.1 of the
# License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
# USA
# Questions? Contact Glen Hansen (Glen.Hansen@inl.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
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
MOERTEL::Segment::Segment(int id, int nnode, int* nodeId, int outlevel) :
Id_(id),
outputlevel_(outlevel),
stype_(MOERTEL::Segment::seg_none)
{
  nodeId_.resize(nnode);
  for (int i=0; i<nnode; ++i) nodeId_[i] = nodeId[i];
  nodeptr_.resize(0);
}

/*----------------------------------------------------------------------*
 | base class constructor                                    mwgee 07/05|
 *----------------------------------------------------------------------*/
MOERTEL::Segment::Segment(int outlevel) :
Id_(-1),
outputlevel_(outlevel),
stype_(MOERTEL::Segment::seg_none)
{
  nodeId_.resize(0);
  nodeptr_.resize(0); 
}

/*----------------------------------------------------------------------*
 | base class copy ctor                                      mwgee 07/05|
 *----------------------------------------------------------------------*/
MOERTEL::Segment::Segment(MOERTEL::Segment& old)
{ 
  Id_          = old.Id_;
  outputlevel_ = old.outputlevel_;
  stype_       = old.stype_; 
  nodeId_      = old.nodeId_;
  nodeptr_     = old.nodeptr_; 
  
  // copy the functions
  // this is not a deep copy but we simply copy the refcountptr
  std::map< int,Teuchos::RCP<MOERTEL::Function> >::iterator curr;
  for (curr = old.functions_.begin(); curr != old.functions_.end(); ++curr)
  {
    if (curr->second == Teuchos::null)
    {
		std::stringstream oss;
		oss << "***ERR*** MOERTEL::Segment::BaseClone(MOERTEL::Segment& old):\n"
           << "***ERR*** function id " << curr->first << " is null\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
    }
	Teuchos::RCP<MOERTEL::Function> newfunc = curr->second;
    functions_.insert(std::pair< int,Teuchos::RCP<MOERTEL::Function> >(curr->first,newfunc));
  }
  
}

/*----------------------------------------------------------------------*
 | base class destructor                                     mwgee 07/05|
 *----------------------------------------------------------------------*/
MOERTEL::Segment::~Segment()
{
  nodeId_.clear();
  nodeptr_.clear();
  functions_.clear();
}

/*----------------------------------------------------------------------*
 |  print segment                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Segment::Print() const
{ 
  std::cout << "Segment " << std::setw(6) << Id_; 
  if (stype_ == MOERTEL::Segment::seg_Linear1D)
    std::cout << " Typ Linear1D   ";
  if (stype_ == MOERTEL::Segment::seg_BiLinearQuad)
    std::cout << " Typ BiLinearQuad ";
  if (stype_ == MOERTEL::Segment::seg_BiLinearTri)
    std::cout << " Typ BiLinearTri";
  if (stype_ == MOERTEL::Segment::seg_none)
    std::cout << " Typ NONE       ";
  std::cout << " #Nodes " << nodeId_.size() << " Nodes: ";
  for (int i=0; i<(int)nodeId_.size(); ++i)
    std::cout << std::setw(6) << nodeId_[i] << "  ";
  std::cout << "  #Functions " << functions_.size() << "  Types: ";
  std::map<int,Teuchos::RCP<MOERTEL::Function> >::const_iterator curr;
  for (curr=functions_.begin(); curr != functions_.end(); ++curr)
    std::cout << curr->second->Type() << "  ";
  std::cout << endl;
  return true;
}

/*----------------------------------------------------------------------*
 | attach a certain shape function to this segment           mwgee 06/05|
 | the user is not supposed to destroy func!                            |
 | the user can set func to several segments!                           |
 *----------------------------------------------------------------------*/
bool MOERTEL::Segment::SetFunction(int id, MOERTEL::Function* func)
{ 
  if (id<0)
  {
	std::cout << "***ERR*** MOERTEL::Segment::SetFunction:\n"
         << "***ERR*** id = " << id << " < 0 (out of range)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!func)
  {
	std::cout << "***ERR*** MOERTEL::Segment::SetFunction:\n"
         << "***ERR*** func = NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  // check for existing function with this id and evtentually overwrite
  std::map<int,Teuchos::RCP<MOERTEL::Function> >::iterator curr = functions_.find(id);
  if (curr != functions_.end())
  {
    curr->second = Teuchos::null;
    curr->second = Teuchos::rcp(func->Clone());
    return true;
  }
  Teuchos::RCP<MOERTEL::Function> newfunc = Teuchos::rcp(func->Clone());
  functions_.insert(std::pair<int,Teuchos::RCP<MOERTEL::Function> >(id,newfunc));
  return true;
}

#if 0
/*----------------------------------------------------------------------*
 | attach a certain shape function to this segment           mwgee 11/05|
 | the user is not supposed to destroy func!                            |
 | the user can set func to several segments!                           |
 *----------------------------------------------------------------------*/
bool MOERTEL::Segment::SetFunction(int id, Teuchos::RCP<MOERTEL::Function> func)
{ 
  if (id<0)
  {
	std::cout << "***ERR*** MOERTEL::Segment::SetFunction:\n"
         << "***ERR*** id = " << id << " < 0 (out of range)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (func==Teuchos::null)
  {
	std::cout << "***ERR*** MOERTEL::Segment::SetFunction:\n"
         << "***ERR*** func = NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  // check for existing function with this id and evtentually overwrite
  std::map<int,Teuchos::RCP<MOERTEL::Function> >::iterator curr = functions_.find(id);
  if (curr != functions_.end())
  {
    curr->second = func;
    return true;
  }
  Teuchos::RCP<MOERTEL::Function> newfunc = func;
  functions_.insert(pair<int,Teuchos::RCP<MOERTEL::Function> >(id,newfunc));
  return true;
}
#endif

/*----------------------------------------------------------------------*
 | evaluate the shape function id at the point xi            mwgee 06/05|
 | id     (in)   id of the function to evaluate                         |
 | xi     (in)   natural coordinates -1<xi<1 where to eval the function |
 | val    (out)  function values, if NULL on input, no evaluation       |
 | valdim (in)   dimension of val                                       |
 | deriv  (out)  derivatives of functions at xi, if NULL on input,      |
 |               no evaluation                                          | 
 *----------------------------------------------------------------------*/
bool MOERTEL::Segment::EvaluateFunction(int id, const double* xi, double* val, 
                                        int valdim, double* deriv)
{ 
  std::map<int,Teuchos::RCP<MOERTEL::Function> >::iterator curr = functions_.find(id);
  if (curr == functions_.end())
  {
	std::stringstream oss;
		oss << "***ERR*** MOERTEL::Segment::EvaluateFunction:\n"
         << "***ERR*** function id " << id << " does not exist on this segment\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }
  curr->second->EvaluateFunction(*this,xi,val,valdim,deriv);
  
  return true;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 06/05|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const MOERTEL::Segment& seg)
{
  seg.Print(); 
  return os;
}

/*----------------------------------------------------------------------*
 | get local numbering id for global node id on this segment mwgee 07/05|
 | return -1 of nid is not adjacent to this segment                     |
 *----------------------------------------------------------------------*/
int MOERTEL::Segment::GetLocalNodeId(int nid)
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
	std::stringstream oss;
		oss << "***ERR*** MOERTEL::Segment::GetLocalNodeId:\n"
         << "***ERR*** cannot find node " << nid << " in segment " << this->Id() << " list of nodes\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
  }
  return lid;
}

/*----------------------------------------------------------------------*
 | build an outward normal at a node adjacent to this        mwgee 07/05|
 *----------------------------------------------------------------------*/
double* MOERTEL::Segment::BuildNormalAtNode(int nid)
{ 
  // find this node in my list of nodes and get local numbering for it
  int lid = GetLocalNodeId(nid);
  
  // depending on what type of segment I am get local coordinates
  double xi[2];
  LocalCoordinatesOfNode(lid,xi);
  
  // build an outward unit normal at xi and return it
  return BuildNormal(xi);
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 | construct ptrs to redundant nodes from my node id list               |
 *----------------------------------------------------------------------*/
bool MOERTEL::Segment::GetPtrstoNodes(MOERTEL::Interface& interface)
{ 
  if (!interface.IsComplete()) return false;
  if (!interface.lComm()) return true;
  if (!nodeId_.size()) return false;
  
  // vector nodeptr_ might already exist, recreate it
  nodeptr_.clear();
  nodeptr_.resize(nodeId_.size());
  
  for (int i=0; i<(int)nodeId_.size(); ++i)
  {
    nodeptr_[i] = interface.GetNodeView(nodeId_[i]).get();
    if (!nodeptr_[i])
    {
		std::stringstream oss;
		oss << "***ERR*** MOERTEL::Segment::GetPtrstoNodes:\n"
           << "***ERR*** interface " << interface.Id() << " GetNodeView failed\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 10/05|
 | construct ptrs to nodes from vector                                  |
 *----------------------------------------------------------------------*/
bool MOERTEL::Segment::GetPtrstoNodes(std::vector<MOERTEL::Node*>& nodes)
{ 
  if (!nodeId_.size()) return false;
  
  // vector nodeptr_ might already exist, recreate it
  nodeptr_.clear();
  nodeptr_.resize(nodeId_.size());
  
  for (int i=0; i<(int)nodeId_.size(); ++i)
  {
    bool foundit = true;
    for (int j=0; j<(int)nodes.size(); ++j)
      if (nodes[j]->Id() == nodeId_[i])
      {
        foundit = true;
        nodeptr_[i] = nodes[j];
        break;
      }
    if (!foundit)
    {
		std::stringstream oss;
		oss << "***ERR*** MOERTEL::Segment::GetPtrstoNodes:\n"
           << "***ERR*** cannot find node " << nodeId_[i] << " in vector\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
	  throw ReportError(oss);
    }
  }
  return true;
}

/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 | return type of function with certain id                              |
 *----------------------------------------------------------------------*/
MOERTEL::Function::FunctionType MOERTEL::Segment::FunctionType(int id)
{ 
  // find the function with id id
  std::map<int,Teuchos::RCP<MOERTEL::Function> >::iterator curr = functions_.find(id);
  if (curr==functions_.end()) 
    return MOERTEL::Function::func_none;
  else
    return curr->second->Type();
}

#if 0
/*----------------------------------------------------------------------*
 |                                                           mwgee 07/05|
 | get ptr to function with id id                                       |
 *----------------------------------------------------------------------*/
MOERTEL::Function* MOERTEL::Segment::GetFunction(int id)
{ 
  // find the function with id id
  std::map<int,Teuchos::RCP<MOERTEL::Function> >::iterator curr = functions_.find(id);
  if (curr==functions_.end()) 
    return NULL;
  else
    return curr->second;
}
#endif
