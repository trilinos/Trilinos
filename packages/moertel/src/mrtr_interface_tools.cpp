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
#include "mrtr_interface.H"
#include "mrtr_projector.H"
#include "mrtr_utils.H"
#include "mrtr_pnode.H"
#include "mrtr_segment_bilineartri.H"


/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
MOERTEL::Interface::Interface(int Id,  bool oneD, Epetra_Comm& comm, int outlevel) :
Id_(Id),
outlevel_(outlevel),
oneD_(oneD),
isComplete_(false),
isIntegrated_(false),
gcomm_(comm),
lcomm_(Teuchos::null),
mortarside_(-1),
ptype_(MOERTEL::Interface::proj_continousnormalfield),
primal_(MOERTEL::Function::func_none),
dual_(MOERTEL::Function::func_none)
{
  return;
}

/*----------------------------------------------------------------------*
 |  copy-ctor (public)                                       mwgee 06/05|
 *----------------------------------------------------------------------*/
MOERTEL::Interface::Interface(MOERTEL::Interface& old) :
Id_(old.Id_),
outlevel_(old.outlevel_),
oneD_(old.oneD_),
isComplete_(old.isComplete_),
isIntegrated_(old.isIntegrated_),
gcomm_(old.gcomm_),
mortarside_(old.mortarside_),
ptype_(old.ptype_),
primal_(old.primal_),
dual_(old.dual_)
{
  // copy the nodes and segments
  for (int i=0; i<2; ++i)
  {
    // the local segment map
	std::map<int,Teuchos::RCP<MOERTEL::Segment> >::const_iterator seg_curr;
    for (seg_curr=old.seg_[i].begin(); seg_curr != old.seg_[i].end(); ++seg_curr)
    {
	  Teuchos::RCP<MOERTEL::Segment>  tmpseg = Teuchos::rcp(seg_curr->second->Clone());
      seg_[i].insert(std::pair<int,Teuchos::RCP<MOERTEL::Segment> >(tmpseg->Id(),tmpseg));
    }
    // the global segment map
    for (seg_curr=old.rseg_[i].begin(); seg_curr != old.rseg_[i].end(); ++seg_curr)
    {
	  Teuchos::RCP<MOERTEL::Segment> tmpseg = Teuchos::rcp(seg_curr->second->Clone());
      rseg_[i].insert(std::pair<int,Teuchos::RCP<MOERTEL::Segment> >(tmpseg->Id(),tmpseg));
    }
    // the local node map
	std::map<int,Teuchos::RCP<MOERTEL::Node> >::const_iterator node_curr;
    for (node_curr=old.node_[i].begin(); node_curr != old.node_[i].end(); ++node_curr)
    {
	  Teuchos::RCP<MOERTEL::Node> tmpnode = Teuchos::rcp( new MOERTEL::Node(*(node_curr->second)));
      node_[i].insert(std::pair<int,Teuchos::RCP<MOERTEL::Node> >(tmpnode->Id(),tmpnode));
    }
    // the global node map
    for (node_curr=old.rnode_[i].begin(); node_curr != old.rnode_[i].end(); ++node_curr)
    {
	  Teuchos::RCP<MOERTEL::Node> tmpnode = Teuchos::rcp(new MOERTEL::Node(*(node_curr->second)));
      rnode_[i].insert(std::pair<int,Teuchos::RCP<MOERTEL::Node> >(tmpnode->Id(),tmpnode));
    }
  }
  // copy the PID maps
  segPID_  = old.segPID_;
  nodePID_ = old.nodePID_;
  
  // copy the local communicator of this interface
  lcomm_ = old.lcomm_;
    
  // rebuild the node-segment topology on this new interface
  BuildNodeSegmentTopology(); 
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
MOERTEL::Interface::~Interface()
{ 
  // delete segments
  for (int i=0; i<2; ++i)
  {
    seg_[i].clear();
    rseg_[i].clear();
  } 
  
  // delete nodes
  for (int i=0; i<2; ++i)
  {
    node_[i].clear();
    rnode_[i].clear();
  } 
  
  // delete PID maps
  segPID_.clear();
  nodePID_.clear();
}

/*----------------------------------------------------------------------*
 |  print segments of this interface to cout                  (public)  |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::PrintSegments() const
{ 
  if (!lComm()) return true;
  
  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::const_iterator curr;
  for (int j=0; j<2; ++j) 
  {
    for (int k=0; k<lComm()->NumProc(); ++k) 
    {
      if (lComm()->MyPID()==k)
      {
        cout << "---global/local Proc " << gcomm_.MyPID() << "/" << k 
             << ":\t Segments Side " << j << endl;
        for (curr=rseg_[j].begin(); curr!=rseg_[j].end(); ++curr)
        {
		  Teuchos::RCP<MOERTEL::Segment> seg = curr->second;
          if (SegPID(seg->Id()) == k)
          {
            if (seg == Teuchos::null)
            {
              cout << "***ERR*** MOERTEL::Interface::PrintSegments:\n"
                   << "***ERR*** found NULL entry in map of segments\n"
                   << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
              return false;
            }
            cout << *seg;
          }
        }
      }
      lComm()->Barrier();
    }
  }
  lComm()->Barrier();
  return true;
}

/*----------------------------------------------------------------------*
 |  print nodes of this interface to cout                     (public)  |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::PrintNodes() const
{ 
  if (!lComm()) return true;
  
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::const_iterator curr;
  
  for (int j=0; j<2; ++j)
  {
    for (int k=0; k<lComm()->NumProc(); ++k)
    {
      if (lComm()->MyPID()==k)
      {
        cout << "---global/local Proc " << gcomm_.MyPID() << "/" << k  
             << ":\t Nodes Side " << j << endl;
        for (curr=rnode_[j].begin(); curr!=rnode_[j].end(); ++curr)
        {
		  Teuchos::RCP<MOERTEL::Node> node = curr->second;
          if (NodePID(node->Id()) == k)
          {
            if (node == Teuchos::null)
            {
              cout << "***ERR*** MOERTEL::Interface::PrintNodes:\n"
                   << "***ERR*** found NULL entry in map of nodes\n"
                   << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
              return false;
            }
            cout << *node;
          }
        }
        lComm()->Barrier();
      }
    }
  }
  lComm()->Barrier();
  return true;
}

/*----------------------------------------------------------------------*
 |  print interface to cout                                   (public)  |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Print() const
{ 
  
  if (!IsComplete())
  {
    if (gcomm_.MyPID()==0)
    {
      cout << "Complete() was NOT called\n";
      cout << "Cannot print node/segment info\n";
    }
    return true;
  }
  if (!lComm()) return true;
  if (gcomm_.MyPID()==0)
  {
    cout << "===== MOERTEL Interface # " << Id() << " =====\n";
    if (oneD_)
      cout << "Dimension: 1D\n";
    else
      cout << "Dimension: 2D\n";
    if (GetProjectionType()==MOERTEL::Interface::proj_none)
      cout << "ProjectionType: none\n";
    else if (GetProjectionType()==MOERTEL::Interface::proj_continousnormalfield)
      cout << "ProjectionType: continousnormalfield\n";
    else if (GetProjectionType()==MOERTEL::Interface::proj_orthogonal)
      cout << "ProjectionType: orthogonal\n";
    int mside = MortarSide();
    int sside = MortarSide();
    if (mside==1 || mside==0)
      sside = OtherSide(mside);
    cout << "Mortar Side " << mside << endl;
    cout << "Slave  Side " << sside << endl;
  }
  
  if (gcomm_.MyPID()==0)
    cout << "----- Segments -----\n";
  fflush(stdout);
  PrintSegments();

  if (gcomm_.MyPID()==0)
    cout << "----- Nodes    -----\n";
  fflush(stdout);
  PrintNodes();

  return true;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 06/05|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const MOERTEL::Interface& inter)
{ 
  inter.Print();
  return (os);
}

/*----------------------------------------------------------------------*
 |  add a single segment to a specified side of the interface (public)  |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::AddSegment(MOERTEL::Segment& seg, int side)
{ 
  // check whether this interface has been finalized before
  if (IsComplete())
  {
    if (OutLevel()>0)
      cout << "***ERR*** MOERTEL::Interface::AddSegment:\n"
           << "***ERR*** Cannot add segment as Complete() was called before\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  // check side
  if (side != 0 && side != 1)
  {
    if (OutLevel()>0)
      cout << "***ERR*** MOERTEL::Interface::AddSegment:\n"
           << "***ERR*** parameter side: " << side << " has to be 0 or 1\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  if (seg.Type()==MOERTEL::Segment::seg_BiLinearQuad)
  {
    if (seg.Nnode() != 4)
    {
      cout << "***ERR*** MOERTEL::Interface::AddSegment:\n"
           << "***ERR*** Unknown number of nodes " << seg.Nnode() << "for BilinearQuad\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
#if 0 // splitting the quad into 2 triangles
    // split the quad into 2 triangles
    int ids1[3];
    ids1[0] = seg.NodeIds()[0];
    ids1[1] = seg.NodeIds()[1];
    ids1[2] = seg.NodeIds()[2];
    int ids2[3];
    ids2[0] = seg.NodeIds()[0];
    ids2[1] = seg.NodeIds()[2];
    ids2[2] = seg.NodeIds()[3];

    // create 2 triangles, give second one the negative id
	Teuchos::RCP<MOERTEL::Segment> tmp1 = Teuchos::rcp( new MOERTEL::Segment_BiLinearTri(seg.Id(),3,ids1,seg.OutLevel()));
	Teuchos::RCP<MOERTEL::Segment> tmp2 = Teuchos::rcp( new MOERTEL::Segment_BiLinearTri(-seg.Id(),3,ids2,seg.OutLevel()));
    
    // add 2 triangles
	std::map<int,Teuchos::RCP<MOERTEL::Segment> >* s = 0;
    if (side==0) s = &(seg_[0]);
    else         s = &(seg_[1]);
    s->insert(std::pair<int,Teuchos::RCP<MOERTEL::Segment> >(tmp1->Id(),tmp1));    
    s->insert(std::pair<int,Teuchos::RCP<MOERTEL::Segment> >(tmp2->Id(),tmp2));    
#else // using a real quad
    // copy the segment
	Teuchos::RCP<MOERTEL::Segment> tmp = Teuchos::rcp( seg.Clone());
    // add segment
	std::map<int,Teuchos::RCP<MOERTEL::Segment> >* s = 0;
    if (side==0) s = &(seg_[0]);
    else         s = &(seg_[1]);
    s->insert(std::pair<int,Teuchos::RCP<MOERTEL::Segment> >(tmp->Id(),tmp));
#endif
  }
  else // all other types of segments
  {
    // copy the segment
	Teuchos::RCP<MOERTEL::Segment> tmp = Teuchos::rcp( seg.Clone());
    // add segment
	std::map<int,Teuchos::RCP<MOERTEL::Segment> >* s = 0;
    if (side==0) s = &(seg_[0]);
    else         s = &(seg_[1]);
    s->insert(std::pair<int,Teuchos::RCP<MOERTEL::Segment> >(tmp->Id(),tmp));
  }

  return true;
}

/*----------------------------------------------------------------------*
 |  add a single node to a specified side of the interface (public)     |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::AddNode(MOERTEL::Node& node, int side)
{ 
  // check whether this interface has been finalized before
  if (IsComplete())
  {
    if (OutLevel()>0)
      cout << "***ERR*** MOERTEL::Interface::AddNode:\n"
           << "***ERR*** Cannot add node as Complete() was called before\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  // check side
  if (side != 0 && side != 1)
  {
    if (OutLevel()>0)
      cout << "***ERR*** MOERTEL::Interface::AddNode:\n"
           << "***ERR*** parameter side: " << side << " has to be 0 or 1\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  // copy the node
  Teuchos::RCP<MOERTEL::Node>  tmp = Teuchos::rcp( new MOERTEL::Node(node));
  
  // add node
  std::map<int,Teuchos::RCP<MOERTEL::Node> >* n = 0;
  if (side==0) n = &(node_[0]);
  else         n = &(node_[1]);
  n->insert(std::pair<int,Teuchos::RCP<MOERTEL::Node> >(tmp->Id(),tmp));

  return true;
}

/*----------------------------------------------------------------------*
 |  set a MOERTEL::function derived class to all segments                  |
 |  on a specified side                                                 |
 |  side      (in)    side of interface to set function to              |
 |  id        (in)    id of the function                                |
 |  func      (in)    ptr to function class to set to segments          |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::SetFunctionAllSegmentsSide(int side, 
                                                 int id, MOERTEL::Function* func)
{ 
  if (side!=0 && side!=1)
  {
    cout << "***ERR*** MOERTEL::Interface::SetFunctionAllSegmentsSide:\n"
         << "***ERR*** side = " << side << " not equal 0 or 1\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (id<0)
  {
    cout << "***ERR*** MOERTEL::Interface::SetFunctionAllSegmentsSide:\n"
         << "***ERR*** id = " << id << " < 0 (out of range)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!func)
  {
    cout << "***ERR*** MOERTEL::Interface::SetFunctionAllSegmentsSide:\n"
         << "***ERR*** func = NULL on input\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  // set the function to my own segments
  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::iterator scurr;
  for (scurr=seg_[side].begin(); scurr!=seg_[side].end(); ++scurr)
    scurr->second->SetFunction(id,func);

  // if redundant segments are already build, set function there as well
  for (scurr=rseg_[side].begin(); scurr!=rseg_[side].end(); ++scurr)
    scurr->second->SetFunction(id,func);
  
  
  return true;
}

/*----------------------------------------------------------------------*
 |  set the Mortar (Master) side of the interface                       |
 |                                                                      |
 |  NOTE: This is a collective call that returns global number of       |
 |        segments of a side over all procs                             |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::SetMortarSide(int side)
{ 
  if (side!=0 && side!=1 && side!=-2)
  {
    cout << "***ERR*** MOERTEL::Interface::SetMortarSide:\n"
         << "***ERR*** side = " << side << " not equal 0 or 1\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    mortarside_=-1;
    return false;
  }
  mortarside_ = side;
  return true;
}

/*----------------------------------------------------------------------*
 |  get number of segments on interface side 0 or 1                     |
 |  side (in) side of which the total number of segments to return      |
 |                                                                      |
 |  NOTE: This is a collective call that returns global number of       |
 |        segments of a side over all procs                             |
 |        participating in the intra-communcicator of this interface.   |
 |        It returns 0 for procs not part of that intra-communcicator   |
 |        Complete() needs to be called before using this method        |
 *----------------------------------------------------------------------*/
int MOERTEL::Interface::GlobalNsegment(int side)
{ 
  if (!IsComplete())
  {
    cout << "***WRN*** MOERTEL::Interface::GlobalNsegment:\n"
         << "***WRN*** Complete() was not called on interface " << Id_ << "\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return -1;
  }
  if (side!=0 && side!=1)
  {
    cout << "***WRN*** MOERTEL::Interface::GlobalNsegment:\n"
         << "***WRN*** side = " << side << " not equal 0 or 1\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return (-1);
  }
  if (!lComm()) 
    return 0;
  int lnsegment = seg_[side].size();
  int gnsegment;
  lcomm_->SumAll(&lnsegment,&gnsegment,1);
  return(gnsegment);
}

/*----------------------------------------------------------------------*
 |  get number of segments on interface (both sides)                    |
 |  NOTE: This is a collective call that returns global number of       |
 |        segments over all procs                                       |
 |        participating in the intra-communcicator of this interface.   |
 |        It returns 0 for procs not part of that intra-communcicator   |
 |        Complete() needs to be called before using this method        |
 *----------------------------------------------------------------------*/
int MOERTEL::Interface::GlobalNsegment()
{ 
  if (!IsComplete())
  {
    cout << "***WRN*** MOERTEL::Interface::GlobalNsegment:\n"
         << "***WRN*** Complete() was not called on interface " << Id_ << "\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return -1;
  }
  if (!lComm()) 
    return 0;
  int nsegment = rseg_[0].size() + rseg_[1].size();
  return(nsegment);
}

/*----------------------------------------------------------------------*
 |  get number of nodes on interface side 0 or 1                     |
 |  side (in) side of which the total number of segments to return      |
 |                                                                      |
 |  NOTE: This is a collective call that returns global number of       |
 |        segments of a side over all procs                             |
 |        participating in the intra-communcicator of this interface.   |
 |        It returns 0 for procs not part of that intra-communcicator   |
 |        Complete() needs to be called before using this method        |
 *----------------------------------------------------------------------*/
int MOERTEL::Interface::GlobalNnode(int side)
{ 
  if (!IsComplete())
  {
    cout << "***WRN*** MOERTEL::Interface::GlobalNnode:\n"
         << "***WRN*** Complete() was not called on interface " << Id_ << "\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return -1;
  }
  if (side!=0 && side!=1)
  {
    cout << "***WRN*** MOERTEL::Interface::GlobalNnode:\n"
         << "***WRN*** side = " << side << " not equal 0 or 1\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return (-1);
  }
  if (!lComm()) 
    return 0;
  int gnnode = rnode_[side].size();
  return(gnnode);
}

/*----------------------------------------------------------------------*
 |  get number of nodes on interface (both sides)                       |
 |  NOTE: This is a collective call that returns global number of       |
 |        segments over all procs                                       |
 |        participating in the intra-communcicator of this interface.   |
 |        It returns 0 for procs not part of that intra-communcicator   |
 |        Complete() needs to be called before using this method        |
 *----------------------------------------------------------------------*/
int MOERTEL::Interface::GlobalNnode()
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::GlobalNnode:\n"
         << "***ERR*** Complete() was not called on interface " << Id_ << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return -1;
  }
  if (!lComm()) 
    return 0;
  int gnnode = rnode_[0].size() + rnode_[1].size();;
  return(gnnode);
}

/*----------------------------------------------------------------------*
 |  find PID (process id) for given nodal id nid                        |
 |  NOTE:                                                               |
 |  this method returns the PID from the intra-communicator of this     |
 |  interface. If called from a proc that is not part of this           |
 |  intra-communicator, the method returns -1                           |
 *----------------------------------------------------------------------*/
int MOERTEL::Interface::NodePID(int nid) const
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::NodePID:\n"
         << "***ERR*** Cannot search node, Complete() not called on interface " << Id_ << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return (-1);
  }
  if (!lComm())
  {
    if (OutLevel()>0)
    cout << "MOERTEL: ***WRN*** MOERTEL::Interface::NodePID:\n"
         << "MOERTEL: ***WRN*** Proc " << gcomm_.MyPID() << " not part of intra-communicator of interface " << Id() << "\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return (-1);
  }
  
  std::map<int,int>::const_iterator curr = nodePID_.find(nid);
  if (curr != nodePID_.end())
    return(curr->second);
  else
  {
    cout << "***ERR*** MOERTEL::Interface::NodePID:\n"
         << "***ERR*** Proc/Intra-Proc " << gcomm_.MyPID() << "/" << lcomm_->MyPID() << ": Cannot find node " << nid << " on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return(-1);
  }
}

/*----------------------------------------------------------------------*
 |  find PID (process id) for given segment id sid                      |
 *----------------------------------------------------------------------*/
int MOERTEL::Interface::SegPID(int sid) const
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::SegPID:\n"
         << "***ERR*** Cannot search segment, Complete() not called on interface " << Id_ << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return (-1);
  }
  if (!lComm())
  {
    if (OutLevel()>0)
    cout << "MOERTEL: ***WRN*** MOERTEL::Interface::NodePID:\n"
         << "MOERTEL: ***WRN*** Proc " << gcomm_.MyPID() << " not part of intra-communicator of interface " << Id() << "\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return (-1);
  }
  
  std::map<int,int>::const_iterator curr = segPID_.find(sid);
  if (curr != segPID_.end())
    return(curr->second);
  else
  {
    cout << "***ERR*** MOERTEL::Interface::SegPID:\n"
         << "***ERR*** Proc/Intra-Proc " << gcomm_.MyPID() << "/"<< lcomm_->MyPID() << ": Cannot find segment " << sid << "on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return(-1);
  }
}

/*----------------------------------------------------------------------*
 |  find PID (process id) for given segment id sid                      |
 *----------------------------------------------------------------------*/
int MOERTEL::Interface::OtherSide(int side) const
{ 
  if (side==0) return 1;
  else if (side==1) return 0;
  else
  {
    cout << "***WRN*** MOERTEL::Interface::OtherSide:\n"
         << "***WRN*** side " << side << " out of range (0 or 1)\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return -1;
  }
}

/*----------------------------------------------------------------------*
 |  get view of a local node with node id nid                           |
 |  if sid is not a local node will return NULL                         |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MOERTEL::Node> MOERTEL::Interface::GetNodeViewLocal(int nid)
{ 
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator curr = node_[0].find(nid);
  if (curr != node_[0].end())
    return(curr->second);
  curr = node_[1].find(nid);
  if (curr != node_[1].end())
    return(curr->second);
  return (Teuchos::null);
}

/*----------------------------------------------------------------------*
 |  get view of a node with node id nid                                 |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MOERTEL::Node> MOERTEL::Interface::GetNodeView(int nid)
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::GetNodeView:\n"
         << "***ERR*** Interface " << Id() << ": Complete() not called\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  if (!lComm()) return Teuchos::null;
  
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator curr = rnode_[0].find(nid);
  if (curr != rnode_[0].end())
    return(curr->second);
  curr = rnode_[1].find(nid);
  if (curr != rnode_[1].end())
    return(curr->second);
  return (Teuchos::null);
}

/*----------------------------------------------------------------------*
 |  get view of ALL nodes on this interface                             |
 |  method allocates a ptr vector. The calling method is in charge of    |
 | destroying it                                                        |
 | returns NULL if proc is not part of the local communicator            |
 *----------------------------------------------------------------------*/
MOERTEL::Node** MOERTEL::Interface::GetNodeView()
{ 
  if (!IsComplete())
  {
    cout << "***WRN*** MOERTEL::Interface::GetNodeView:\n"
         << "***WRN*** Interface " << Id() << ": Complete() not called\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return NULL;
  }
  if (!lComm()) return NULL;
  
  MOERTEL::Node** view = new MOERTEL::Node*[GlobalNnode()];
  int count=0;
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator curr;
  for (int i=0; i<2; ++i)
    for (curr=rnode_[i].begin(); curr != rnode_[i].end(); ++curr)
    {
      view[count] = curr->second.get();
      ++count;
    }
  return view;
}

/*----------------------------------------------------------------------*
 |  get view of ALL nodes on this interface                             |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::GetNodeView(std::vector<MOERTEL::Node*>& nodes)
{ 
  if (!IsComplete())
  {
    cout << "***WRN*** MOERTEL::Interface::GetNodeView:\n"
         << "***WRN*** Interface " << Id() << ": Complete() not called\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!lComm()) return false;
  
  nodes.resize(GlobalNnode());
  int count=0;
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator curr;
  for (int i=0; i<2; ++i)
    for (curr=rnode_[i].begin(); curr != rnode_[i].end(); ++curr)
    {
      nodes[count] = curr->second.get();
      ++count;
    }
  return true;
}

/*----------------------------------------------------------------------*
 |  get view of a local segment with id sid                             |
 *----------------------------------------------------------------------*/
Teuchos::RCP<MOERTEL::Segment>  MOERTEL::Interface::GetSegmentView(int sid)
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::GetSegmentView:\n"
         << "***ERR*** Interface " << Id() << ": Complete() not called\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  if (!lComm()) return Teuchos::null;
  
  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::iterator curr = rseg_[0].find(sid);
  if (curr != rseg_[0].end())
    return(curr->second);
  curr = rseg_[1].find(sid);
  if (curr != rseg_[1].end())
    return(curr->second);
  return (Teuchos::null);
}

/*----------------------------------------------------------------------*
 |  get view of ALL segments on this interface                          |
 |  method returns a ptr to a vector, calling method is in              |
 | charge of deleteing it                                               |
 *----------------------------------------------------------------------*/
MOERTEL::Segment** MOERTEL::Interface::GetSegmentView()
{ 
  if (!IsComplete())
  {
    cout << "***WRN*** MOERTEL::Interface::GetSegmentView:\n"
         << "***WRN*** Interface " << Id() << ": Complete() not called\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return NULL;
  }
  if (!lComm()) return NULL;
  
  MOERTEL::Segment** segs = new MOERTEL::Segment*[GlobalNsegment()];
  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::iterator curr;
  int count=0;
  for (int i=0; i<2; ++i)
    for (curr=rseg_[i].begin(); curr != rseg_[i].end(); ++curr)
    {
      segs[count] = curr->second.get();
      ++count;
    }
  return segs;
}

/*----------------------------------------------------------------------*
 |  find out which side a segment is on                                 |
 | returns -1 if it can't find the segment on either side               |
 *----------------------------------------------------------------------*/
int MOERTEL::Interface::GetSide(MOERTEL::Segment* seg)
{ 
  if (!lComm()) return -1;
  if (!IsComplete())
  {
    cout << "***WRN*** MOERTEL::Interface::GetSide:\n"
         << "***WRN*** Interface " << Id() << ": Complete() not called\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return -1;
  }
  if (!lComm())
  {
    cout << "***WRN*** MOERTEL::Interface::GetSide:\n"
         << "***WRN*** Interface " << Id() << ": Proc " << gcomm_.MyPID() << "not in intra-comm\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return -1;
  }
  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::iterator curr = rseg_[0].find(seg->Id());
  if (curr != rseg_[0].end())
    return(0);
  curr = rseg_[1].find(seg->Id());
  if (curr != rseg_[1].end())
    return(1);
  return (-1);
}

/*----------------------------------------------------------------------*
 |  find out which side a node is on                                    |
 | returns -1 if it can't find the node on either side                  |
 *----------------------------------------------------------------------*/
int MOERTEL::Interface::GetSide(MOERTEL::Node* node)
{ 
  if (!IsComplete())
  {
    cout << "***WRN*** MOERTEL::Interface::GetSide:\n"
         << "***WRN*** Interface " << Id() << ": Complete() not called\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return -1;
  }
  if (!lComm())
  {
    cout << "***WRN*** MOERTEL::Interface::GetSide:\n"
         << "***WRN*** Interface " << Id() << ": Proc " << gcomm_.MyPID() << "not in intra-comm\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return -1;
  }
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator curr = rnode_[0].find(node->Id());
  if (curr != rnode_[0].end())
    return(0);
  curr = rnode_[1].find(node->Id());
  if (curr != rnode_[1].end())
    return(1);
  return (-1);
}

/*----------------------------------------------------------------------*
 |  find out which side a node is on                                    |
 | returns -1 if it can't find the node on either side                  |
 *----------------------------------------------------------------------*/
int MOERTEL::Interface::GetSide(int nodeid)
{ 
  if (!IsComplete())
  {
    cout << "***WRN*** MOERTEL::Interface::GetSide:\n"
         << "***WRN*** Interface " << Id() << ": Complete() not called\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return -1;
  }
  if (!lComm())
  {
    cout << "***WRN*** MOERTEL::Interface::GetSide:\n"
         << "***WRN*** Interface " << Id() << ": Proc " << gcomm_.MyPID() << "not in intra-comm\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return -1;
  }
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator curr = rnode_[0].find(nodeid);
  if (curr != rnode_[0].end())
    return(0);
  curr = rnode_[1].find(nodeid);
  if (curr != rnode_[1].end())
    return(1);
  return (-1);
}

/*----------------------------------------------------------------------*
 | allreduce all segments of a side                                     |
 |                                                                      |
 | NOTE: this is a collective call of all procs in the intra-comm       |
 |       After call to RedundantNodes and RedundantSegments             |
 |       a call to BuildNodeSegmentTopology is necessary to complete    |
 |       the construction of redundant nodes/segments
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::RedundantSegments(int side)
{ 
  if (side != 0 && side != 1)
  {
    cout << "***ERR*** MOERTEL::Interface::RedundantSegments:\n"
         << "***ERR*** side=" << side << " out of range (0 or 1)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::RedundantSegments:\n"
         << "***ERR*** Complete() not called on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return (-1);
  }
  
  // send everybody who doesn't belong here out of here
  if (!lComm())
    return true;
  
  std::map<int,Teuchos::RCP<MOERTEL::Segment> >* rmap = &(rseg_[side]);
  // check whether redundant map has been build before
  if (rmap->size() != 0)
    return true;

  // add my own segments to the redundant map
  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::const_iterator curr;
  for (curr=seg_[side].begin(); curr != seg_[side].end(); ++curr)
  {
    //MOERTEL::Segment* tmp = curr->second->Clone();
    // FIXME: is this ok? it's not a deep copy anymore.....
	Teuchos::RCP<MOERTEL::Segment> tmp = curr->second;
    rmap->insert(std::pair<int,Teuchos::RCP<MOERTEL::Segment> >(curr->first,tmp));
  }
  
  // loop over all procs and broadcast proc's segments
  for (int proc=0; proc<lcomm_->NumProc(); ++proc)
  {
    int nseg = 0;
    if (proc==lcomm_->MyPID())
      nseg = MyNsegment(side);
    lcomm_->Broadcast(&nseg,1,proc);

    int  bsize = nseg*12;
	std::vector<int> bcast; 
    
    // pack proc's segments
    if (proc==lcomm_->MyPID())
    {
       int  count = 0;
       bcast.resize(bsize);
       for (curr=seg_[side].begin(); curr != seg_[side].end(); ++curr)
       {
         int  numint;
         int* spack = curr->second->Pack(&numint);
         if (count+numint>=bsize)
         {
           bsize += 5* numint;
           bcast.resize(bsize);
         }
         for (int i=0; i<numint; ++i)
           bcast[count++] = spack[i];
         delete [] spack;
       }
       bsize = count;
    }
    
    // broadcast proc's segments
	
    lcomm_->Broadcast(&bsize,1,proc); // Communicate the size of data needed to store the segments

    if (lcomm_->MyPID() != proc)
      bcast.resize(bsize);

	if(bsize > 0){ // Only send the segment(s) if there are one or more

		lcomm_->Broadcast(&bcast[0],bsize,proc);
    
		// Unpack proc's segments
		if (lcomm_->MyPID() != proc) {

			int count=0;
			for (int i=0; i<nseg; ++i) {

				// the type of segment is stored second in the pack
				MOERTEL::Segment* tmp = AllocateSegment(bcast[count+1],OutLevel());
				tmp->UnPack(&(bcast[count]));
				Teuchos::RCP<MOERTEL::Segment> tmp2 = Teuchos::rcp(tmp);
				count += bcast[count];
				rmap->insert(std::pair<int,Teuchos::RCP<MOERTEL::Segment> >(tmp2->Id(),tmp2));
			}
		}
	}

    bcast.clear();

  } // for (int proc=0; proc<lcomm_->NumProc(); ++proc)
  return true;
}


/*----------------------------------------------------------------------*
 | allreduce all nodes of a side                                        |
 |                                                                      |
 | NOTE: this is a collective call of all procs in the intra-comm       |
 |       After call to RedundantNodes and RedundantSegments             |
 |       a call to BuildNodeSegmentTopology is necessary to complete    |
 |       the construction of redundant nodes/segments
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::RedundantNodes(int side)
{ 
  if (side != 0 && side != 1)
  {
    cout << "***ERR*** MOERTEL::Interface::RedundantNodes:\n"
         << "***ERR*** side=" << side << " out of range (0 or 1)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::RedundantNodes:\n"
         << "***ERR*** Complete() not called on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return (-1);
  }
  
  // send everybody who doesn't belong here out of here
  if (!lComm())
    return true;

  std::map<int,Teuchos::RCP<MOERTEL::Node> >* rmap = &(rnode_[side]);
  // check whether redundant map has been build before
  if (rmap->size() != 0)
    return true;

  // add my own nodes to the redundant map
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::const_iterator curr;
  for (curr=node_[side].begin(); curr != node_[side].end(); ++curr)
  {
    //MOERTEL::Node* tmp = new MOERTEL::Node(*(curr->second));
    //FIXME: this is not a deep copy anymore. Is this ok?
	Teuchos::RCP<MOERTEL::Node> tmp = curr->second;
    rmap->insert(std::pair<int,Teuchos::RCP<MOERTEL::Node> >(curr->first,tmp));
  }
  
  // loop all procs and broadcast proc's nodes
  for (int proc=0; proc<lcomm_->NumProc(); ++proc)
  {
    int nnode = 0;
    if (proc==lcomm_->MyPID())
      nnode = MyNnode(side);
    lcomm_->Broadcast(&nnode,1,proc);
    
    int bsize = nnode*3;
	std::vector<double> bcast;
    
    // pack proc's nodes
    if (proc==lcomm_->MyPID())
    {
      int count = 0;
      bcast.resize(bsize);
      for (curr=node_[side].begin(); curr != node_[side].end(); ++curr)
      {
        int numdouble;
        double* npack = curr->second->Pack(&numdouble);
        if (count+numdouble>=bsize)
        {
          bsize += 3*numdouble;
          bcast.resize(bsize);
        }
        for (int i=0; i<numdouble; ++i)
          bcast[count++] = npack[i];
        delete [] npack;
      }
      bsize = count;
    }
    
    // bcast proc's nodes
	
    lcomm_->Broadcast(&bsize,1,proc);

    if (lcomm_->MyPID() != proc)
      bcast.resize(bsize);

	if(bsize > 0){ // Only send the segment(s) if there are one or more

		lcomm_->Broadcast(&bcast[0],bsize,proc);
    
		// Unpack proc's nodes
		if (lcomm_->MyPID() != proc) {

			int count=0;
			for (int i=0; i<nnode; ++i) {

				Teuchos::RCP<MOERTEL::Node> tmp = Teuchos::rcp(new MOERTEL::Node(OutLevel()));
				tmp->UnPack(&(bcast[count]));
				count += (int)bcast[count];
				rmap->insert(std::pair<int,Teuchos::RCP<MOERTEL::Node> >(tmp->Id(),tmp));
			}
		}    
	}

    bcast.clear();
  } // for (int proc=0; proc<lcomm_->NumProc(); ++proc)
  return true;
}

/*----------------------------------------------------------------------*
 | (re)build the topology info between nodes and segments               |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::BuildNodeSegmentTopology()
{ 
  if (!IsComplete())
  {
    if (OutLevel()>1)
    cout << "MOERTEL: ***WRN*** MOERTEL::Interface::BuildNodeSegmentTopology:\n"
         << "MOERTEL: ***WRN*** Complete() not called on interface " << Id() << "\n"
         << "MOERTEL: ***WRN*** Cannot build node<->segment topology\n"
         << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  if (!lComm()) return true;
  
  // loop nodes and find their adjacent segments
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator ncurr;
  for (int side=0; side<2; ++side)
  {
    for (ncurr=rnode_[side].begin(); ncurr != rnode_[side].end(); ++ncurr)
      ncurr->second->GetPtrstoSegments(*this);
  }
  
  // loop segments and find their adjacent nodes
  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::iterator scurr;
  for (int side=0; side<2; ++side)
  {
    for (scurr=rseg_[side].begin(); scurr != rseg_[side].end(); ++scurr)
      scurr->second->GetPtrstoNodes(*this);
  }
  return true;
}

/*----------------------------------------------------------------------*
 | set lagrange multiplier dofs starting from minLMGID                  |
 | to all slave nodes or segments that have a projection                |
 | On exit, return maxLMGID+1, where maxLMGID is the last LM dof number |
 | used here                                                            |
 | Note that this is collective for ALL procs                           |
 *----------------------------------------------------------------------*/
int MOERTEL::Interface::SetLMDofs(int minLMGID)
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::SetLMDofs:\n"
         << "***ERR*** Complete() was not called on interface " << Id_ << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return (0);
  }

  if (lComm())
  {
    int mside = MortarSide();
    int sside = OtherSide(mside);
    
    // create redundant flags
	std::vector<int> lhavelm(rnode_[sside].size());
	std::vector<int> ghavelm(rnode_[sside].size());
    for (int i=0; i<(int)rnode_[sside].size(); ++i) lhavelm[i] = 0;
    
    // loop through redundant nodes and add my flags
    int count=0;
	std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator curr;
    for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)
    {
      if (NodePID(curr->second->Id()) != lComm()->MyPID())
      {
        ++count;
        continue;
      }
	  Teuchos::RCP<std::map<int,double> > D = curr->second->GetD();
      if (D==Teuchos::null)
      {
        if (curr->second->GetM() != Teuchos::null)
          cout << *curr->second << "has no D but has M!!!\n";
        ++count;
        continue;
      }
      lhavelm[count] = 1;
      ++count;
    }
    if (count != (int)rnode_[sside].size())
    {
      cout << "***ERR*** MOERTEL::Interface::SetLMDofs:\n"
           << "***ERR*** number of redundant nodes wrong\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return (0);
    }

    // make the flags redundant
    lComm()->MaxAll(&lhavelm[0],&ghavelm[0],lhavelm.size());
    lhavelm.clear();

    // loop through nodes again and set lm dofs according to ghavelm
    count=0;
    for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)
    {
      if (!ghavelm[count])
      {
        ++count;
        continue;
      }

      // get number of dofs on this node to choose the same number of dofs
      // for the LM
      int ndof = curr->second->Ndof();
      
      // set LM dofs to this node and it's projection
      for (int i=0; i<ndof; ++i)
      {
        curr->second->SetLagrangeMultiplierId(minLMGID+i);
		Teuchos::RCP<MOERTEL::ProjectedNode> pnode = curr->second->GetProjectedNode();
        if (pnode!=Teuchos::null)
          pnode->SetLagrangeMultiplierId(minLMGID+i);
      }
      minLMGID += ndof;
      ++count;
    }
      ghavelm.clear();
  } // if (lComm())
  
  // broadcast minLMGID to all procs including those not in intra-comm
  int lbcaster = 0;
  int gbcaster = 0;
  if (lComm())
    if (lComm()->MyPID()==0)
      lbcaster = gcomm_.MyPID();
  gcomm_.MaxAll(&lbcaster,&gbcaster,1);
  gcomm_.Broadcast(&minLMGID,1,gbcaster);
  return(minLMGID);
}



#if 0
/*----------------------------------------------------------------------*
 | set lagrange multiplier dofs starting from minLMGID                  |
 | to all slave nodes or segments that have a projection                |
 | On exit, return maxLMGID+1, where maxLMGID is the last LM dof number |
 | used here                                                            |
 | Note that this is collective for ALL procs                           |
 *----------------------------------------------------------------------*/
int MOERTEL::Interface::SetLMDofs(int minLMGID)
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::SetLMDofs:\n"
         << "***ERR*** Complete() was not called on interface " << Id_ << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return (0);
  }

  if (lComm())
  {
    int mside = MortarSide();
    int sside = OtherSide(mside);
    
    if (IsOneDimensional())
    {
      // loop nodes on slave side and set LMdofs for those who have a projection
	  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator curr;
      for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)
      {
        // check whether this node has a projection
		Teuchos::RCP<MOERTEL::ProjectedNode> pnode = curr->second->GetProjectedNode();
        if (pnode==Teuchos::null) continue;
        
        
        // get number of dofs on this node to choose the same number of dofs
        // for the LM
        int ndof = curr->second->Ndof();
        
        // set LM dofs to this node and it's projection
        for (int i=0; i<ndof; ++i)
        {
          curr->second->SetLagrangeMultiplierId(minLMGID+i);
          pnode->SetLagrangeMultiplierId(minLMGID+i);
        }
        minLMGID += ndof;
      } // for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)
    } // if (IsOneDimensional())
    else
    {
      // create redundant flags
	  std::vector<int> lhavelm(rnode_[sside].size());
	  std::vector<int> ghavelm(rnode_[sside].size());
      for (int i=0; i<(int)rnode_[sside].size(); ++i) lhavelm[i] = 0;
      
      // loop through redundant nodes and add my flags
      int count=0;
	  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator curr;
      for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)
      {
        if (NodePID(curr->second->Id()) != lComm()->MyPID())
        {
          ++count;
          continue;
        }
		Teuchos::RCP<std::map<int,double> > D = curr->second->GetD();
        if (D==Teuchos::null)
        {
          if (curr->second->GetM() != Teuchos::null)
            cout << *curr->second << "has no D but has M!!!\n";
          ++count;
          continue;
        }
        lhavelm[count] = 1;
        ++count;
      }
      if (count != (int)rnode_[sside].size())
      {
        cout << "***ERR*** MOERTEL::Interface::SetLMDofs:\n"
             << "***ERR*** number of redundant nodes wrong\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return (0);
      }

      // make the flags redundant
      lComm()->MaxAll(&lhavelm[0],&ghavelm[0],lhavelm.size());
      lhavelm.clear();

      // loop through nodes again and set lm dofs according to ghavelm
      count=0;
      for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)
      {
        if (!ghavelm[count])
        {
          ++count;
          continue;
        }

        // get number of dofs on this node to choose the same number of dofs
        // for the LM
        int ndof = curr->second->Ndof();
        
        // set LM dofs to this node and it's projection
        for (int i=0; i<ndof; ++i)
        {
          curr->second->SetLagrangeMultiplierId(minLMGID+i);
		  Teuchos::RCP<MOERTEL::ProjectedNode> pnode = curr->second->GetProjectedNode();
          if (pnode!=Teuchos::null)
            pnode->SetLagrangeMultiplierId(minLMGID+i);
        }
        minLMGID += ndof;
        ++count;
      }
      ghavelm.clear();
    }
  } // if (lComm())
  
  // broadcast minLMGID to all procs including those not in intra-comm
  int lbcaster = 0;
  int gbcaster = 0;
  if (lComm())
    if (lComm()->MyPID()==0)
      lbcaster = gcomm_.MyPID();
  gcomm_.MaxAll(&lbcaster,&gbcaster,1);
  gcomm_.Broadcast(&minLMGID,1,gbcaster);
  return(minLMGID);
}
#endif

/*----------------------------------------------------------------------*
 | retrieve a vector containing a list of lm ids owned by this processor|
 | The calling routine is responsible for destroying this list          |
 *----------------------------------------------------------------------*/
std::vector<int>* MOERTEL::Interface::MyLMIds()
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::MyLMIds:\n"
         << "***ERR*** Complete() was not called on interface " << Id_ << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return (0);
  }

  int mside = MortarSide();
  int sside = OtherSide(mside);

  // allocate a vector with a guess
  std::vector<int>* lmids = new std::vector<int>;
  
  // procs not in intra-comm return an empty vector
  if (!lComm())
  {
    lmids->resize(0);
    return lmids;
  }
  
  lmids->resize(rnode_[sside].size()*10);
  int count=0;
    
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator curr;
  for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)
  {
	Teuchos::RCP<MOERTEL::Node> node = curr->second;
    if (NodePID(node->Id()) != lComm()->MyPID()) 
      continue;
    int  nlmdof = node->Nlmdof();
    if (!nlmdof) 
      continue; 
    const int* ids = node->LMDof();
    if (count+nlmdof>(int)lmids->size())
      lmids->resize(lmids->size()+50*nlmdof);
    for (int i=0; i<nlmdof; ++i)
      (*lmids)[count++] = ids[i];
  }
  lmids->resize(count);


  return lmids;
}

/*----------------------------------------------------------------------*
 | detect end segments and reduce the order of the lm shape functions   |
 | on these end segments                                                |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::DetectEndSegmentsandReduceOrder()
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::DetectEndSegmentsandReduceOrder:\n"
         << "***ERR*** Complete() was not called on interface " << Id_ << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!lComm()) return true;

  if (IsOneDimensional())
    return DetectEndSegmentsandReduceOrder_2D();
  else
    return DetectEndSegmentsandReduceOrder_3D();
  return false;
}

#if 0 // old version
/*----------------------------------------------------------------------*
 | detect end segments and reduce the order of the lm shape functions   |
 | on these end segments                                                |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::DetectEndSegmentsandReduceOrder_2D()
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::DetectEndSegmentsandReduceOrder:\n"
         << "***ERR*** Complete() was not called on interface " << Id_ << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!lComm()) return true;

  if (!IsOneDimensional())
    return false;

  int mside = MortarSide();
  int sside = OtherSide(mside);

  /*
  in the 1D interface case, an end segment is detected as follows:
  - an end segments is attached to a node that has a projection and is the 
    ONLY segment attached to that node
  - an end segment is attached to a node that has a pseudo projection (that is
    a node carrying lagrange mutlipliers but not having a projection) 
  */
  
  // loop all nodes on the slave side and find those with only one segment
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator curr;
  for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)
  {
	Teuchos::RCP<MOERTEL::Node> node = curr->second;
    bool foundit = false;
    if (node->Nseg()<2)
      foundit = true;
    if (node->GetProjectedNode() != Teuchos::null)
      if (!(node->GetProjectedNode()->Segment()))
        foundit = true;
    if (!foundit)
      continue;
    
    MOERTEL::Segment** segs = node->Segments();

    for (int i=0; i<node->Nseg(); ++i)
    {
      MOERTEL::Function::FunctionType type = 
        segs[i]->FunctionType(1);
        
      MOERTEL::Function_Constant1D* tmp1 = NULL;  
      switch (type)
      {
        // for linear and dual linear reduce function order to constant
        case MOERTEL::Function::func_Constant1D:
        case MOERTEL::Function::func_Linear1D:
        case MOERTEL::Function::func_DualLinear1D:
          tmp1 = new Function_Constant1D(OutLevel());
          segs[i]->SetFunction(1,tmp1);
        break;
        case MOERTEL::Function::func_none:
          cout << "***ERR*** MOERTEL::Interface::DetectEndSegmentsandReduceOrder:\n"
               << "***ERR*** interface " << Id() << " function type of function 1 on segment " << segs[0]->Id() << " is func_none\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          exit(EXIT_FAILURE);     
        break;
        default:
          cout << "***ERR*** MOERTEL::Interface::DetectEndSegmentsandReduceOrder:\n"
               << "***ERR*** interface " << Id() << " function type of function 1 on segment " << segs[0]->Id() << " is unknown\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          exit(EXIT_FAILURE);     
        break;
      } // switch (type)
      if (tmp1) delete tmp1; tmp1 = NULL;
    } // for (int i=0; i<node->Nseg(); ++i)
  } // for (curr=rnode_[sside].begin(); curr!=rnode_[sside].end(); ++curr)
  return true;
}
#endif

/*----------------------------------------------------------------------*
 | detect end segments and reduce the order of the lm shape functions   |
 | on these end segments                                                |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::DetectEndSegmentsandReduceOrder_2D()
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::DetectEndSegmentsandReduceOrder:\n"
         << "***ERR*** Complete() was not called on interface " << Id_ << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!lComm()) return true;

  if (!IsOneDimensional())
    return false;

  int mside = MortarSide();
  if (mside != 1 && mside != 0)
  {
    cout << "***ERR*** MOERTEL::Interface::DetectEndSegmentsandReduceOrder_2D:\n"
         << "***ERR*** Mortar side not set on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  int sside = OtherSide(mside);

  // A node attached to only one element AND on the boundary is
  // considered a corner node and is member of ONE support set
  // It is in the modified support psi tilde of the closest internal node
  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator ncurr;
  for (ncurr=rnode_[sside].begin(); ncurr != rnode_[sside].end(); ++ncurr)
  {
    if (!(ncurr->second->IsOnBoundary())) continue;
    MOERTEL::Segment** seg   = ncurr->second->Segments();
    MOERTEL::Node**    nodes = seg[0]->Nodes();
    // find the supporting node for this potential corner node on same element
    for (int i=0; i<seg[0]->Nnode(); ++i)
    {
      if (nodes[i]->Id() == ncurr->second->Id()) continue;
      if (!nodes[i]->IsOnBoundary())
      {
        //cout << "Supporting neighbor node on same element is \n" << *nodes[i];
        ncurr->second->AddSupportedByNode(nodes[i]);
      }
    }
    if (!ncurr->second->NSupportSet())
    {
      cout << "***ERR*** MOERTEL::Interface::DetectEndSegmentsandReduceOrder_2D:\n"
           << "***ERR*** Cannot find a supporting internal node for corner node " << ncurr->second->Id() << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }
  return true;
}


/*----------------------------------------------------------------------*
  // prepare the boundary modification for 3D interfaces
  // Nodes on the edge of an interface will not carry LMs so they
  // do not conflict with other interfaces
  // The choice of the Mortar side in the case of several interfaces
  // is then arbitrary
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::DetectEndSegmentsandReduceOrder_3D()
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::DetectEndSegmentsandReduceOrder_3D:\n"
         << "***ERR*** Complete() was not called on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  if (!lComm()) return true;

  if (IsOneDimensional())
    return false;

  int mside = MortarSide();
  if (mside != 1 && mside != 0)
  {
    cout << "***ERR*** MOERTEL::Interface::DetectEndSegmentsandReduceOrder_3D:\n"
         << "***ERR*** Mortar side not set on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  int sside = OtherSide(mside);
  
  // 1
  // A node attached to only one element AND on the boundary is
  // considered a corner node and is member of ONE support set
  // It is in the modified support psi tilde of the closest internal node
  
  // 2
  // A node on the boundary that is not a corner node is member of 
  // so many support sets as the # of internal nodes it is topologically 
  // connected to
  
  // See B.Wohlmuth:"Discretization Methods and Iterative Solvers
  //                 Based on Domain Decomposition", pp 33/34, Springer 2001.

  std::map<int,Teuchos::RCP<MOERTEL::Node> >::iterator ncurr;

  // do 1
  for (ncurr=rnode_[sside].begin(); ncurr != rnode_[sside].end(); ++ncurr)
  {
    if (!(ncurr->second->IsOnBoundary())) continue;
    if (ncurr->second->Nseg() != 1) continue;
    MOERTEL::Segment** seg   = ncurr->second->Segments();
    MOERTEL::Node**    nodes = seg[0]->Nodes();
    // find the supporting node for this potential corner node on same element
    for (int i=0; i<seg[0]->Nnode(); ++i)
    {
      if (nodes[i]->Id() == ncurr->second->Id()) continue;
      if (!nodes[i]->IsOnBoundary())
      {
        //cout << "Supporting neighbor node on same element is \n" << *nodes[i];
        ncurr->second->AddSupportedByNode(nodes[i]);
      }
    }
    if (!ncurr->second->NSupportSet()) // we've not found a supporting node yet
    {
      for (int i=0; i<seg[0]->Nnode(); ++i)
      {
        if (nodes[i]->Id() == ncurr->second->Id()) continue;
        MOERTEL::Segment** neighborsegs = nodes[i]->Segments();
        for (int j=0; j<nodes[i]->Nseg(); ++j)
        {
          MOERTEL::Node** neighborneighbornodes = neighborsegs[j]->Nodes();
          for (int k=0; k<neighborsegs[j]->Nnode(); ++k)
            if (!neighborneighbornodes[k]->IsOnBoundary())
            {
              ncurr->second->AddSupportedByNode(neighborneighbornodes[k]);
            }
        }
      }
    }
    if (!ncurr->second->NSupportSet())
    {
      cout << "***ERR*** MOERTEL::Interface::DetectEndSegmentsandReduceOrder_3D:\n"
           << "***ERR*** Cannot find a supporting internal node for corner node " << ncurr->second->Id() << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }
  
  // do 2
  for (ncurr=rnode_[sside].begin(); ncurr != rnode_[sside].end(); ++ncurr)
  {
    // do only nodes on boundary that don't have a supporting node yet
    if (!(ncurr->second->IsOnBoundary()) || ncurr->second->NSupportSet()) continue;
    //cout << "Looking at boundary node  " << *(ncurr->second);
    // loop all segments adjacent to this node
    MOERTEL::Segment** segs = ncurr->second->Segments();
    for (int i=0; i<ncurr->second->Nseg(); ++i)
    {
      MOERTEL::Node** neighbornodes = segs[i]->Nodes();
      for (int j=0; j<segs[i]->Nnode(); ++j)
      {
        if (neighbornodes[j]->IsOnBoundary()) continue;
        //cout << "Supporting neighbor node on same element is \n" << *neighbornodes[j];
        ncurr->second->AddSupportedByNode(neighbornodes[j]);
      }
    }
    if (!ncurr->second->NSupportSet())
    {
      cout << "***ERR*** MOERTEL::Interface::DetectEndSegmentsandReduceOrder_3D:\n"
           << "***ERR*** Cannot find a supporting internal node for boundary node " << ncurr->second->Id() << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }


#if 0
  MOERTEL::Function_ConstantTri* tmp = new Function_ConstantTri(OutLevel());
  for (ncurr=rnode_[sside].begin(); ncurr != rnode_[sside].end(); ++ncurr)
  {
    if (!(ncurr->second->IsOnBoundary())) continue;
    MOERTEL::Segment** seg   = ncurr->second->Segments();
    for (int i=0; i<ncurr->second->Nseg(); ++i)
      seg[i]->SetFunction(1,tmp);    
  }  
  if (tmp) delete tmp;  
#endif
  
  
  return true;
}

/*----------------------------------------------------------------------*
 | set the type of shape functions that are supposed to be used         |
 | primal: type of shape function for the trace space                   |
 | dual:   type of shape function for the LM space                      |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::SetFunctionTypes(MOERTEL::Function::FunctionType primal,
                                       MOERTEL::Function::FunctionType dual)
{ 
  primal_ = primal;
  dual_   = dual;
  return true;
}

/*----------------------------------------------------------------------*
 | set functions to all segments depending on the variables primal_     |
 | and dual_
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::SetFunctionsFromFunctionTypes()
{ 
  if (!lComm()) return true;
  if (!IsComplete())
  {
    cout << "***ERR*** MOERTEL::Interface::SetFunctionsFromFunctionTypes:\n"
         << "***ERR*** interface " << Id() << " : Complete() was not called\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);     
  }
  if (dual_==MOERTEL::Function::func_none)
  {
    cout << "***ERR*** MOERTEL::Interface::SetFunctionsFromFunctionTypes:\n"
         << "***ERR*** interface " << Id() << " : no dual function type set\n"
         << "***ERR*** use SetFunctionTypes(..) to set function types\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);     
  }
  
#if 0
  if (primal_==MOERTEL::Function::func_BiLinearQuad) primal_ = MOERTEL::Function::func_LinearTri;
  if (dual_==MOERTEL::Function::func_BiLinearQuad) dual_ = MOERTEL::Function::func_LinearTri;
  if (dual_==MOERTEL::Function::func_DualBiLinearQuad) dual_ = MOERTEL::Function::func_DualLinearTri;
#endif

  // set the primal shape functions
  MOERTEL::Function_Linear1D*         func1 = NULL;
  MOERTEL::Function_Constant1D*       func2 = NULL;
  MOERTEL::Function_DualLinear1D*     func3 = NULL;
  MOERTEL::Function_LinearTri*        func4 = NULL;
  MOERTEL::Function_DualLinearTri*    func5 = NULL;
  MOERTEL::Function_ConstantTri*      func6 = NULL;
  MOERTEL::Function_BiLinearQuad*     func7 = NULL;
  MOERTEL::Function_DualBiLinearQuad* func8 = NULL;
  switch(primal_)
  {
    case MOERTEL::Function::func_Linear1D:
      func1 = new MOERTEL::Function_Linear1D(OutLevel());
      SetFunctionAllSegmentsSide(0,0,func1);
      SetFunctionAllSegmentsSide(1,0,func1);
      delete func1; func1 = NULL;
    break;
    case MOERTEL::Function::func_DualLinear1D:
      cout << "MOERTEL: ***ERR*** MOERTEL::Interface::SetFunctionsFromFunctionTypes:\n"
           << "MOERTEL: ***ERR*** interface " << Id() << " : setting discontious dual shape functions as\n"
           << "MOERTEL: ***ERR*** primal isoparametric trace space function is probably a bad idea...\n"
           << "MOERTEL: ***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);     
    break;
    case MOERTEL::Function::func_Constant1D:
      cout << "MOERTEL: ***ERR*** MOERTEL::Interface::SetFunctionsFromFunctionTypes:\n"
           << "MOERTEL: ***ERR*** interface " << Id() << " : setting constant shape functions as\n"
           << "MOERTEL: ***ERR*** primal isoparametric trace space function is probably a bad idea...\n"
           << "MOERTEL: ***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);     
    break;
    case MOERTEL::Function::func_LinearTri:
      func4 = new MOERTEL::Function_LinearTri(OutLevel());
      SetFunctionAllSegmentsSide(0,0,func4);
      SetFunctionAllSegmentsSide(1,0,func4);
      delete func4; func4 = NULL;
    break;
    case MOERTEL::Function::func_BiLinearQuad:
      func7 = new MOERTEL::Function_BiLinearQuad(OutLevel());
      SetFunctionAllSegmentsSide(0,0,func7);
      SetFunctionAllSegmentsSide(1,0,func7);
      delete func7; func7 = NULL;
    break;
    case MOERTEL::Function::func_none:
      cout << "***ERR*** MOERTEL::Interface::SetFunctionsFromFunctionTypes:\n"
           << "***ERR*** interface " << Id() << " : no primal function type set\n"
           << "***ERR*** use SetFunctionTypes(..) to set function types\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);     
    break;
    default:
      cout << "***ERR*** MOERTEL::Interface::SetFunctionsFromFunctionTypes:\n"
           << "***ERR*** interface " << Id() << " : Unknown function type: " << primal_ << endl
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);     
    break;
  }
  
  int side = MortarSide();
  if (side != 1 && side != 0)
  {
    cout << "***ERR*** MOERTEL::Interface::SetFunctionsFromFunctionTypes:\n"
         << "***ERR*** interface " << Id() << " : Mortar Side not set set\n"
         << "***ERR*** use SetMortarSide(int side) to choose mortar side first\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);     
  }
  side = OtherSide(side); // we want the slave side here

  switch(dual_)
  {
    case MOERTEL::Function::func_Linear1D:
      func1 = new MOERTEL::Function_Linear1D(OutLevel());
      SetFunctionAllSegmentsSide(side,1,func1);
      delete func1; func1 = NULL;
    break;
    case MOERTEL::Function::func_DualLinear1D:
      func3 = new MOERTEL::Function_DualLinear1D(OutLevel());
      SetFunctionAllSegmentsSide(side,1,func3);
      delete func3; func3 = NULL;
    break;
    case MOERTEL::Function::func_Constant1D:
      func2 = new MOERTEL::Function_Constant1D(OutLevel());
      SetFunctionAllSegmentsSide(side,1,func2);
      delete func2; func2 = NULL;
    break;
    case MOERTEL::Function::func_LinearTri:
      func4 = new MOERTEL::Function_LinearTri(OutLevel());
      SetFunctionAllSegmentsSide(side,1,func4);
      delete func4; func4 = NULL;
    break;
    case MOERTEL::Function::func_DualLinearTri:
      func5 = new MOERTEL::Function_DualLinearTri(OutLevel());
      SetFunctionAllSegmentsSide(side,1,func5);
      delete func5; func5 = NULL;
    break;
    case MOERTEL::Function::func_ConstantTri:
      func6 = new MOERTEL::Function_ConstantTri(OutLevel());
      SetFunctionAllSegmentsSide(side,1,func6);
      delete func6; func6 = NULL;
    break;
    case MOERTEL::Function::func_BiLinearQuad:
      func7 = new MOERTEL::Function_BiLinearQuad(OutLevel());
      SetFunctionAllSegmentsSide(side,1,func7);
      delete func7; func7 = NULL;
    break;
    case MOERTEL::Function::func_DualBiLinearQuad:
      func8 = new MOERTEL::Function_DualBiLinearQuad(OutLevel());
      SetFunctionAllSegmentsSide(side,1,func8);
      delete func8; func8 = NULL;
    break;
    case MOERTEL::Function::func_none:
      cout << "***ERR*** MOERTEL::Interface::SetFunctionsFromFunctionTypes:\n"
           << "***ERR*** interface " << Id() << " : no dual function type set\n"
           << "***ERR*** use SetFunctionTypes(..) to set function types\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);     
    break;
    default:
      cout << "***ERR*** MOERTEL::Interface::SetFunctionsFromFunctionTypes:\n"
           << "***ERR*** interface " << Id() << " : Unknown function type: " << dual_ << endl
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);     
    break;
  }

  return true;
}
