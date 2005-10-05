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


#include "mrtr_interface.H"
#include "mrtr_projector.H"
#include "mrtr_utils.H"
#include "mrtr_pnode.H"


/*----------------------------------------------------------------------*
 |  build averaged normals and make projection of nodes                 |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::Project()
{ 
  bool ok = false;
  
  //-------------------------------------------------------------------
  // interface needs to be complete
  if (!IsComplete())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MRTR::Interface::Project:\n"
           << "***ERR*** Complete() not called on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // send all procs not member of this interface's intra-comm out of here
  if (!lComm()) return true;

  //-------------------------------------------------------------------
  // interface segments need to have at least one function on each side
  map<int,MRTR::Segment*>::iterator curr;
  for (int side=0; side<2; ++side)
    for (curr=rseg_[side].begin(); curr!=rseg_[side].end(); ++curr)
      if (curr->second->Nfunctions() < 1)
      {
        cout << "***ERR*** MRTR::Interface::Project:\n"
             << "***ERR*** interface " << Id_ << ", mortar side\n"
             << "***ERR*** segment " << curr->second->Id() << " needs at least 1 function set\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    
  //-------------------------------------------------------------------
  // build nodal normals on both sides
  map<int,MRTR::Node*>::iterator ncurr;
  for (int side=0; side<2; ++side)
    for (ncurr=rnode_[side].begin(); ncurr!=rnode_[side].end(); ++ncurr)
    {
#if 1
      cout << *(ncurr->second);
#endif
      ncurr->second->BuildAveragedNormal();
#if 1
      cout << *(ncurr->second);
#endif
    }

  //-------------------------------------------------------------------
  // check the type of projection to be used and project nodes
  // projection along the normal field:
  // uses the slave side's interpolated normal field to
  // project slave nodes onto the master surfaces.
  // Then projects master nodes along the same normal field on slave
  // surfaces. Overrides the normals in the master nodes by the negative
  // of the normal field of their projection point
  if (GetProjectionType() == MRTR::Interface::proj_continousnormalfield)
  {
    ok = ProjectNodes_NormalField();
    if (!ok) return false;
  }
  else if (GetProjectionType() == MRTR::Interface::proj_orthogonal)
  {
    ok = ProjectNodes_Orthogonal();
    if (!ok) return false;
  }
  else
  {
    cout << "***ERR*** MRTR::Interface::Project:\n"
         << "***ERR*** interface " << Id() << "\n"
         << "***ERR*** unknown type of nodal projection\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  return true;
}

/*----------------------------------------------------------------------*
 | do projection of nodes on master and slave side                      |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::ProjectNodes_NormalField()
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MRTR::Interface::ProjectNodes_NormalField:\n"
         << "***ERR*** Complete() not called on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  if (!lComm()) return true;

  // project the slave nodes onto the master surface along the slave normal field
  ProjectNodes_SlavetoMaster_NormalField();  
  
  // project the master nodes ontp the slave surface along slave's normal field
  ProjectNodes_MastertoSlave_NormalField();

  return true;
}

/*----------------------------------------------------------------------*
 | project the slave nodes onto master segments along slave normal field|
 *----------------------------------------------------------------------*/
bool MRTR::Interface::ProjectNodes_SlavetoMaster_NormalField()
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MRTR::Interface::ProjectNodes_SlavetoMaster_NormalField:\n"
         << "***ERR*** Complete() not called on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  if (!lComm()) return true;
  
  int mside = MortarSide();
  int sside = OtherSide(mside);

  // iterate over all nodes of the slave side and project those belonging to me
  map<int,MRTR::Node*>::iterator scurr;
  for (scurr=rnode_[sside].begin(); scurr!=rnode_[sside].end(); ++scurr)
  {
    MRTR::Node* snode = scurr->second;
    if (NodePID(snode->Id()) != lComm()->MyPID())
      continue;
    
    const double* sx = snode->X();
    double mindist = 1.0e+20;
    MRTR::Node* closenode = NULL;
    
    // find a node on the master side, that is closest to me
    map<int,MRTR::Node*>::iterator mcurr;
    for (mcurr=rnode_[mside].begin(); mcurr!=rnode_[mside].end(); ++mcurr)
    {
      MRTR::Node* mnode = mcurr->second;
      const double* mx = mnode->X();
      
      // build distance | mnode->X() - snode->X() |
      double dist = 0.0;
      for (int i=0; i<3; ++i) dist += (mx[i]-sx[i])*(mx[i]-sx[i]);
      dist = sqrt(dist);
      if (dist <= mindist)
      {
        mindist = dist;
	closenode = mnode;
      }
      //cout << "snode " << snode->Id() << " mnode " << mnode->Id() << " mindist " << mindist  << " dist " << dist << endl;
    }
    if (!closenode)
    {
      cout << "***ERR*** MRTR::Interface::ProjectNodes_SlavetoMaster_NormalField:\n"
           << "***ERR*** Weired: for slave node " << snode->Id() << " no closest master node found\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }

#if 0
    cout << "snode " << snode->Id() << " closenode " << closenode->Id() << endl;
    cout << "snode\n" << *snode;
    cout << "closenode\n" << *closenode;
#endif    

    // get segments attached to closest node cnode
    int  nseg = closenode->Nseg();
    MRTR::Segment** segs = closenode->Segments(); 
    
    // create a projection-iterator
    MRTR::Projector projector(IsOneDimensional());

    // finding a good geometric projection is somehow 
    // critical. We work with some tolerance here and pick the 'best'
    // out of all acceptable projections made    
    // loop these segments and project onto them along snode's normal vector
    double bestdist[2];
    bestdist[0] = bestdist[1] = 1.0e+20;
    MRTR::Segment* bestseg = NULL;
    for (int i=0; i<nseg; ++i)
    {
      // project the slave node onto that master segment
      double xi[2]; xi[0] = xi[1] = 0.0;
      projector.ProjectNodetoSegment_NodalNormal(*snode,*(segs[i]),xi);
      
      // check whether xi is better then previous projections
      if (IsOneDimensional())
      {
        if (abs(xi[0]) < abs(bestdist[0])) 
	{ 
	  bestdist[0] = xi[0];
	  bestdist[1] = xi[1];
	  bestseg = segs[i];
	}
      }
      else
      {
	double dist  = sqrt(xi[0]*xi[0]+xi[1]*xi[1]);
	double bdist = sqrt(bestdist[0]*bestdist[0]+bestdist[1]*bestdist[1]); 
        // it's better in both directions
	if (abs(xi[0])<abs(bestdist[0]) && abs(xi[1])<abs(bestdist[1]))
	{
	  bestdist[0] = xi[0];
	  bestdist[1] = xi[1];
	  bestseg = segs[i];
	}
	//it's better in one direction and 'in' in the other
	else if ((abs(xi[0])<abs(bestdist[0]) && abs(xi[1]) < 1.2) ||
	         (abs(xi[1])<abs(bestdist[1]) && abs(xi[0]) < 1.2) )
	{
	  bestdist[0] = xi[0];
	  bestdist[1] = xi[1];
	  bestseg = segs[i];
	}
	// it's better by radius
	else if (dist<bdist)
	{
	  bestdist[0] = xi[0];
	  bestdist[1] = xi[1];
	  bestseg = segs[i];
	}
      }
    } // for (int i=0; i<nseg; ++i)
    
    // check whether the bestseg and bestdist are inside the segment
    // (with some tolerance of 10%)
    bool ok = false;
    if (IsOneDimensional())
      if (abs(bestdist[0]) < 1.2) ok = true;
    else
      if (abs(bestdist[0])<1.2 && abs(bestdist[1])<1.2) ok = true;
    
    if (ok)  // the projection is good
    {
      // create a projected node and store it in snode
      MRTR::ProjectedNode* pnode 
        = new MRTR::ProjectedNode(*snode,bestdist,bestseg);
      snode->SetProjectedNode(pnode);
    }
  } // for (scurr=rnode_[sside].begin(); scurr!=rnode_[sside].end(); ++scurr)
  lComm()->Barrier();
  
  // loop all slave nodes again and make the projections redundant
  double* bcast = new double[4*rnode_[sside].size()]; // that's the max
  for (int proc=0; proc<lComm()->NumProc(); ++proc)
  {
    int blength = 0;
    if (proc==lComm()->MyPID())
    {
      for (scurr=rnode_[sside].begin(); scurr!=rnode_[sside].end(); ++scurr)
      {
        MRTR::Node* snode = scurr->second;
        if (proc != NodePID(snode->Id())) continue; // I cannot have a projection on a node not owned by me
        MRTR::ProjectedNode* pnode = snode->GetProjectedNode();
        if (!pnode) continue; // this node does not have a projection
        const double* xi = pnode->Xi();
        bcast[blength] = (double)pnode->Id();            
        ++blength;
        if (pnode->Segment())
          bcast[blength] = (double)pnode->Segment()->Id(); 
        else
          bcast[blength] = -1.0; // indicating this node does not have projection but lagrange multipliers
        ++blength;
        bcast[blength] = xi[0];
        ++blength;
        bcast[blength] = xi[1];
        ++blength;
      } 
      if (blength > 4*rnode_[sside].size())
      {
        cout << "***ERR*** MRTR::Interface::ProjectNodes_SlavetoMaster_NormalField:\n"
             << "***ERR*** Overflow in communication buffer occured\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
    }
    lComm()->Broadcast(&blength,1,proc);
    lComm()->Broadcast(bcast,blength,proc);
    if (proc!=lComm()->MyPID())
    {
      int i;
      for (i=0; i<blength;)
      {
        int     nid = (int)bcast[i]; ++i;
        int     sid = (int)bcast[i]; ++i;
        double* xi  = &bcast[i];     ++i; ++i;
        MRTR::Node* snode = GetNodeView(nid);
        MRTR::Segment* seg = NULL;
        if (sid!=-1)
          seg = GetSegmentView(sid);
        if (!snode)
        {
          cout << "***ERR*** MRTR::Interface::ProjectNodes_SlavetoMaster_NormalField:\n"
               << "***ERR*** Cannot get view of node or segment\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          exit(EXIT_FAILURE);
        }
        MRTR::ProjectedNode* pnode = new MRTR::ProjectedNode(*snode,xi,seg);
        snode->SetProjectedNode(pnode);
      }
      if (i != blength)
      {
        cout << "***ERR*** MRTR::Interface::ProjectNodes_SlavetoMaster_NormalField:\n"
             << "***ERR*** Mismatch in dimension of recv buffer\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
    }
  } // for (int proc=0; proc<lComm()->NumProc(); ++proc)
  delete [] bcast; bcast = NULL;
  lComm()->Barrier();

#if 1
  // Postprocess the projections
  // The slave side of the interface might be larger then the master side
  // of the interface so not all slave nodes have a projection.
  // For those slave nodes without a projection attached to a slave segment
  // which overlaps with the master side, lagrange mutlipliers have to be
  // introduced. This is done by checking all nodes without a projection 
  // whether they are attached to some slave segment on which another node
  // HAS a projection. If this case is found, a pseudo ProjectedNode is 
  // introduced for that node.
  for (scurr=rnode_[sside].begin(); scurr!=rnode_[sside].end(); ++scurr)
  {
    MRTR::Node* snode = scurr->second;
    
    // don't do anything on nodes that already have a projection
    if (snode->GetProjectedNode())
      continue;

    // get segments adjacent to this node  
    int nseg             = snode->Nseg();
    MRTR::Segment** segs = snode->Segments();
    
    // loop segments and check for other nodes with projection
    bool foundit = false;
    for (int i=0; i<nseg; ++i)
    {
      int nnode = segs[i]->Nnode();
      MRTR::Node** nodes = segs[i]->Nodes();
      for (int j=0; j<nnode; ++j)
        if (nodes[j]->GetProjectedNode())
          if (nodes[j]->GetProjectedNode()->Segment())
          {
            foundit = true;
            break;
          }
      if (foundit) break;
    }
    
    if (foundit)
    {
#if 0
      cout << "Node without projection:\n" << *snode;        
      cout << "...get's lagrange multipliers\n\n";
#endif
      MRTR::ProjectedNode* pnode = new MRTR::ProjectedNode(*snode,NULL,NULL);
      snode->SetProjectedNode(pnode);
    }
  } // for (scurr=rnode_[sside].begin(); scurr!=rnode_[sside].end(); ++scurr)
#endif  

  return true;
}


/*----------------------------------------------------------------------*
 | project nodes master to slave along slave cont. normal field         |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::ProjectNodes_MastertoSlave_NormalField()
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MRTR::Interface::ProjectNodes_MastertoSlave_NormalField:\n"
         << "***ERR*** Complete() not called on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  if (!lComm()) return true;
  
  int mside = MortarSide();
  int sside = OtherSide(mside);

  // iterate over all nodes of the master side and project those belonging to me
  map<int,MRTR::Node*>::iterator mcurr;
  for (mcurr=rnode_[mside].begin(); mcurr!=rnode_[mside].end(); ++mcurr)
  {
    MRTR::Node* mnode = mcurr->second;
    if (NodePID(mnode->Id()) != lComm()->MyPID())
      continue;
    
    const double* mx = mnode->X();
    double mindist = 1.0e+20;
    MRTR::Node* closenode = NULL;
    
    // find a node on the slave side that is closest to me
    map<int,MRTR::Node*>::iterator scurr;
    for (scurr=rnode_[sside].begin(); scurr!=rnode_[sside].end(); ++scurr)
    {
      MRTR::Node* snode = scurr->second;
      const double* sx = snode->X();
      
      // build distance | snode->X() - mnode->X() |
      double dist = 0.0;
      for (int i=0; i<3; ++i) dist += (mx[i]-sx[i])*(mx[i]-sx[i]);
      dist = sqrt(dist);
      if (dist < mindist)
      {
        mindist = dist;
        closenode = snode;
      }
    } 
    if (!closenode)
    {
      cout << "***ERR*** MRTR::Interface::ProjectNodes_MastertoSlave_NormalField:\n"
           << "***ERR*** Weired: for master node " << mnode->Id() << " no closest master node found\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    //cout << "mnode " << mnode->Id() << " closenode " << closenode->Id() << endl;
    
    // get segments attached to closest node closenode
    int  nseg = closenode->Nseg();
    MRTR::Segment** segs = closenode->Segments(); 
    
    // create a projection operator
    MRTR::Projector projector(IsOneDimensional());
    
    // loop these segments and find best projection
    double bestdist[2];
    bestdist[0] = bestdist[1] = 1.0e+20;
    MRTR::Segment* bestseg = NULL;
    for (int i=0; i<nseg; ++i)
    {
      // project the master node on the slave segment along the segments interpolated normal field
      double xi[2]; xi[0] = xi[1] = 0.0;
      projector.ProjectNodetoSegment_SegmentNormal(*mnode,*(segs[i]),xi);
      
      // check whether xi is better then previous projections
      if (IsOneDimensional())
      {
        if (abs(xi[0]) < abs(bestdist[0])) 
	{ 
	  bestdist[0] = xi[0];
	  bestdist[1] = xi[1];
	  bestseg = segs[i];
	}
      }
      else
      {
	double dist  = sqrt(xi[0]*xi[0]+xi[1]*xi[1]);
	double bdist = sqrt(bestdist[0]*bestdist[0]+bestdist[1]*bestdist[1]); 
        // it's better in both directions
	if (abs(xi[0])<abs(bestdist[0]) && abs(xi[1])<abs(bestdist[1]))
	{
	  bestdist[0] = xi[0];
	  bestdist[1] = xi[1];
	  bestseg = segs[i];
	}
	//it's better in one direction and 'in' in the other
	else if ((abs(xi[0])<abs(bestdist[0]) && abs(xi[1]) < 1.2) ||
	         (abs(xi[1])<abs(bestdist[1]) && abs(xi[0]) < 1.2) )
	{
	  bestdist[0] = xi[0];
	  bestdist[1] = xi[1];
	  bestseg = segs[i];
	}
	// it's better by radius
	else if (dist<bdist)
	{
	  bestdist[0] = xi[0];
	  bestdist[1] = xi[1];
	  bestseg = segs[i];
	}
      }
    } // for (int i=0; i<nseg; ++i)
    
    // check whether the bestseg/bestdist are inside that segment
    // (with some tolerance of 20%)
    bool ok = false;
    if (IsOneDimensional())
      if (abs(bestdist[0]) < 1.1) ok = true;
    else
      if (abs(bestdist[0])<1.1 && abs(bestdist[1])<1.1) ok = true;
    
    if (ok) // the projection is good
    {
      // build the interpolated normal and overwrite the mnode normal with -n
      int          nsnode = bestseg->Nnode();
      MRTR::Node** snodes = bestseg->Nodes();
      double* val   = new double[nsnode];
      bestseg->EvaluateFunction(0,bestdist,val,nsnode,NULL);
      double NN[3]; NN[0] = NN[1] = NN[2] = 0.0;
      for (int i=0; i<nsnode; ++i)
      {
        const double* N = snodes[i]->N();
        for (int j=0; j<3; ++j)
          NN[j] -= val[i]*N[j];
      }
      delete [] val; val = NULL;
      mnode->SetN(NN);

      // create projected node and store it in mnode
      MRTR::ProjectedNode* pnode
        = new MRTR::ProjectedNode(*mnode,bestdist,bestseg);
      mnode->SetProjectedNode(pnode);
    }
    else // this mnode does not have a valid projection
    {
      if (OutLevel()>5)
      cout << "***WRN***: Node " << mnode->Id() << " does not have projection\n\n";
      mnode->SetProjectedNode(NULL);
    }
  } // for (scurr=rnode_[mside].begin(); scurr!=rnode_[mside].end(); ++scurr)

  // loop all master nodes again and make the projection and the new normal redundant
  int bsize = 7*rnode_[mside].size();
  double* bcast = new double[bsize]; // that's the max
  for (int proc=0; proc<lComm()->NumProc(); ++proc)
  {
    int blength = 0;
    if (proc==lComm()->MyPID())
    {
      for (mcurr=rnode_[mside].begin(); mcurr!=rnode_[mside].end(); ++mcurr)
      {
        MRTR::Node* mnode = mcurr->second;
        if (proc != NodePID(mnode->Id())) continue; // cannot have a projection on a node i don't own
        MRTR::ProjectedNode* pnode = mnode->GetProjectedNode();
        if (!pnode) continue; // this node does not have a projection
        const double* xi = pnode->Xi();
        const double* N  = mnode->N();
        bcast[blength] = (double)pnode->Id();
        ++blength;
        bcast[blength] = (double)pnode->Segment()->Id();
        ++blength;
        bcast[blength] = xi[0];
        ++blength;
        bcast[blength] = xi[1];
        ++blength;
        bcast[blength] = N[0];
        ++blength;
        bcast[blength] = N[1];
        ++blength;
        bcast[blength] = N[2];
        ++blength;
      }
      if (blength > bsize)
      {
        cout << "***ERR*** MRTR::Interface::ProjectNodes_MastertoSlave_NormalField:\n"
             << "***ERR*** Overflow in communication buffer occured\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
    }
    lComm()->Broadcast(&blength,1,proc);
    lComm()->Broadcast(bcast,blength,proc);
    if (proc!=lComm()->MyPID())
    {
      int i;
      for (i=0; i<blength;)
      {
        int     nid = (int)bcast[i];  ++i;
        int     sid = (int)bcast[i];  ++i;
        double* xi  =      &bcast[i]; ++i; ++i; 
        double* n   =      &bcast[i]; ++i; ++i; ++i;
        MRTR::Node*    mnode = GetNodeView(nid);
        MRTR::Segment* seg   = GetSegmentView(sid);
        if (!mnode || !seg)
        {
          cout << "***ERR*** MRTR::Interface::ProjectNodes_MastertoSlave_NormalField:\n"
               << "***ERR*** Cannot get view of node or segment\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          exit(EXIT_FAILURE);
        }
        mnode->SetN(n);
        MRTR::ProjectedNode* pnode = new MRTR::ProjectedNode(*mnode,xi,seg);
        mnode->SetProjectedNode(pnode);
      }
      if (i != blength)
      {
        cout << "***ERR*** MRTR::Interface::ProjectNodes_MastertoSlave_NormalField:\n"
             << "***ERR*** Mismatch in dimension of recv buffer\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
    }
  } // for (int proc=0; proc<lComm()->NumProc(); ++proc)
  delete [] bcast; bcast = NULL;

  return true;
}

/*----------------------------------------------------------------------*
 | do projection of nodes on master and slave side                      |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::ProjectNodes_Orthogonal()
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MRTR::Interface::ProjectNodes_Orthogonal:\n"
         << "***ERR*** Complete() not called on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  if (!lComm()) return true;

  // project the master nodes onto the slave surface orthogonaly
  ProjectNodes_MastertoSlave_Orthogonal();

  // project the slave nodes onto the master surface orthogonal to adjacent slave segment
  ProjectNodes_SlavetoMaster_Orthogonal();  
  
  return true;
}

/*----------------------------------------------------------------------*
 | project nodes master to slave along slave orthogonal         |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::ProjectNodes_MastertoSlave_Orthogonal()
{ 
  if (!IsComplete())
  {
    cout << "***ERR*** MRTR::Interface::ProjectNodes_MastertoSlave_Orthogonal:\n"
         << "***ERR*** Complete() not called on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  if (!lComm()) return true;
  
  int mside = MortarSide();
  int sside = OtherSide(mside);

  // iterate over all master nodes and project those belonging to me
    map<int,MRTR::Node*>::iterator mcurr;
  for (mcurr=rnode_[mside].begin(); mcurr!=rnode_[mside].end(); ++mcurr)
  {
    MRTR::Node* mnode = mcurr->second;
    if (NodePID(mnode->Id()) != lComm()->MyPID())
      continue;
      
    const double* mx = mnode->X();
    double mindist = 1.0e+20;
    MRTR::Node* closenode = NULL;
    
    // find a node on the slave side that is closest to me
    map<int,MRTR::Node*>::iterator scurr;
    for (scurr=rnode_[sside].begin(); scurr!=rnode_[sside].end(); ++scurr)
    {
      MRTR::Node* snode = scurr->second;
      const double* sx = snode->X();
      
      // build distance | snode->X() - mnode->X() |
      double dist = 0.0;
      for (int i=0; i<3; ++i) dist += (mx[i]-sx[i])*(mx[i]-sx[i]);
      dist = sqrt(dist);
      if (dist < mindist)
      {
        mindist = dist;
        closenode = snode;
      }
    } 
    if (!closenode)
    {
      cout << "***ERR*** MRTR::Interface::ProjectNodes_MastertoSlave_Orthogonal:\n"
           << "***ERR*** Weired: for master node " << mnode->Id() << " no closest master node found\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }


    // get segments attached to closest node closenode
    int  nseg = closenode->Nseg();
    MRTR::Segment** segs = closenode->Segments(); 
    
    // create a projection operator
    MRTR::Projector projector(IsOneDimensional());
    
    // loop these segments and find best projection
    double bestdist[2];
    bestdist[0] = bestdist[1] = 1.0e+20;
    MRTR::Segment* bestseg = NULL;
    for (int i=0; i<nseg; ++i)
    {
      // project the master node orthogonally on the slave segment
      double xi[2]; xi[0] = xi[1] = 0.0;
      projector.ProjectNodetoSegment_SegmentOrthogonal(*mnode,*(segs[i]),xi);
      
      // check whether xi is better than previous projection
      if (IsOneDimensional())
      {
        if (abs(xi[0]) < abs(bestdist[0])) 
        {
	  bestdist[0] = xi[0];
	  bestdist[1] = xi[1];
	  bestseg = segs[i];
        }
      }
      else
      {
        cout << "***ERR*** MRTR::Interface::ProjectNodes_MastertoSlave_Orthogonal:\n"
             << "***ERR*** not impl. for 3D\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
      
    } // for (int i=0; i<nseg; ++i)
      
    // check whether this best projection is good
    bool ok = false;
    if (IsOneDimensional())
    {
      if (abs(bestdist[0]) < 1.01) 
        ok = true;
    }
    else
    {
      cout << "***ERR*** MRTR::Interface::ProjectNodes_MastertoSlave_Orthogonal:\n"
           << "***ERR*** not impl. for 3D\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    
    if (ok) // the projection is good
    {
      // create a projected node and store it in mnode
      MRTR::ProjectedNode* pnode = 
        new MRTR::ProjectedNode(*mnode,bestdist,bestseg);
      mnode->SetProjectedNode(pnode);
    } 
    else // this mnode does not have a valid projection
    {
      if (OutLevel()>5)
      cout << "***WRN***: Node " << mnode->Id() << " does not have projection\n\n";
      //mnode->SetProjectedNode(NULL);
    }
  } // for (mcurr=rnode_[mside].begin(); mcurr!=rnode_[mside].end(); ++mcurr)

  // loop all master nodes again and make projection redundant
  int bsize = 4*rnode_[mside].size();
  double* bcast = new double[bsize]; // that's the max
  for (int proc=0; proc<lComm()->NumProc(); ++proc)
  {
    int blength = 0;
    if (proc==lComm()->MyPID())
    {
      for (mcurr=rnode_[mside].begin(); mcurr!=rnode_[mside].end(); ++mcurr)
      {
        MRTR::Node* mnode = mcurr->second;
        if (proc != NodePID(mnode->Id())) continue; // cannot have a projection on a node i don't own
        MRTR::ProjectedNode* pnode = mnode->GetProjectedNode();
        if (!pnode) continue; // this node does not have a projection
        const double* xi = pnode->Xi();
        bcast[blength] = (double)pnode->Id();
        ++blength;
        bcast[blength] = (double)pnode->Segment()->Id();
        ++blength;
        bcast[blength] = xi[0];
        ++blength;
        bcast[blength] = xi[1];
        ++blength;
      } // for (mcurr=rnode_[mside].begin(); mcurr!=rnode_[mside].end(); ++mcurr)
      if (blength>bsize)
      {
        cout << "***ERR*** MRTR::Interface::ProjectNodes_MastertoSlave_Orthogonal:\n"
             << "***ERR*** Overflow in communication buffer occured\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
    } // if (proc==lComm()->MyPID())
    lComm()->Broadcast(&blength,1,proc);
    lComm()->Broadcast(bcast,blength,proc);
    if (proc!=lComm()->MyPID())
    {
      int i;
      for (i=0; i<blength;)
      {
        int     nid = (int)bcast[i];  ++i;
        int     sid = (int)bcast[i];  ++i;
        double* xi  =      &bcast[i]; ++i; ++i; 
        MRTR::Node*    mnode = GetNodeView(nid);
        MRTR::Segment* seg   = GetSegmentView(sid);
        if (!mnode || !seg)
        {
          cout << "***ERR*** MRTR::Interface::ProjectNodes_MastertoSlave_Orthogonal:\n"
               << "***ERR*** Cannot get view of node or segment\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          exit(EXIT_FAILURE);
        }
        MRTR::ProjectedNode* pnode = new MRTR::ProjectedNode(*mnode,xi,seg);
        mnode->SetProjectedNode(pnode);
      }
      if (i != blength)
      {
        cout << "***ERR*** MRTR::Interface::ProjectNodes_MastertoSlave_Orthogonal:\n"
             << "***ERR*** Mismatch in dimension of recv buffer: " << i << " != " << blength << "\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
    } // if (proc!=lComm()->MyPID())
  }  // for (int proc=0; proc<lComm()->NumProc(); ++proc)
  delete [] bcast; bcast = NULL;

  return true;
}

/*----------------------------------------------------------------------*
 | project the slave nodes onto master segments orthogonal              |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::ProjectNodes_SlavetoMaster_Orthogonal()
{
  if (!IsComplete())
  {
    cout << "***ERR*** MRTR::Interface::ProjectNodes_SlavetoMaster_Orthogonal:\n"
         << "***ERR*** Complete() not called on interface " << Id() << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  if (!lComm()) return true;
  
  int mside = MortarSide();
  int sside = OtherSide(mside);

  // iterate over all nodes of the slave side and project those belonging to me
  map<int,MRTR::Node*>::iterator scurr;
  for (scurr=rnode_[sside].begin(); scurr!=rnode_[sside].end(); ++scurr)
  {
    MRTR::Node* snode = scurr->second;

#if 0
    cout << "now projecting\n " << *snode;
#endif    
    
    if (NodePID(snode->Id()) != lComm()->MyPID())
      continue;
    
    const double* sx = snode->X();
    double mindist = 1.0e+20;
    MRTR::Node* closenode = NULL;
    
    // find a node on the master side, that is closest to me
    map<int,MRTR::Node*>::iterator mcurr;
    for (mcurr=rnode_[mside].begin(); mcurr!=rnode_[mside].end(); ++mcurr)
    {
      MRTR::Node* mnode = mcurr->second;
      const double* mx = mnode->X();
      
      // build distance | mnode->X() - snode->X() |
      double dist = 0.0;
      for (int i=0; i<3; ++i) dist += (mx[i]-sx[i])*(mx[i]-sx[i]);
      dist = sqrt(dist);
      if (dist <= mindist)
      {
        mindist = dist;
	closenode = mnode;
      }
      //cout << "snode " << snode->Id() << " mnode " << mnode->Id() << " mindist " << mindist  << " dist " << dist << endl;
    }
    if (!closenode)
    {
      cout << "***ERR*** MRTR::Interface::ProjectNodes_SlavetoMaster_Orthogonal:\n"
           << "***ERR*** Weired: for slave node " << snode->Id() << " no closest master node found\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }

    // get segments attached to closest node cnode
    int  nmseg = closenode->Nseg();
    MRTR::Segment** msegs = closenode->Segments(); 
    
    // create a projection-iterator
    MRTR::Projector projector(IsOneDimensional());

    // loop segments and find all projections onto them
    int nsseg = snode->Nseg();
    MRTR::Segment** ssegs = snode->Segments(); 
    for (int i=0; i<nmseg; ++i)
    {
      // loop all segments that are adjacent to the slave node
      for (int j=0; j<nsseg; ++j)
      {
        // project the slave node onto that master segment
        double xi[2]; xi[0] = xi[1] = 0.0;
        projector.ProjectNodetoSegment_Orthogonal_to_Slave(*snode,*(msegs[i]),xi,*(ssegs[j]));
        
        // check whether this projection is good
        bool ok = false;
        if (IsOneDimensional())
        {
          if (abs(xi[0])<1.01)
            ok = true;
        }
        else
        {
          cout << "***ERR*** MRTR::Interface::ProjectNodes_SlavetoMaster_Orthogonal:\n"
               << "***ERR*** not impl. for 3D\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          exit(EXIT_FAILURE);
        }
        
        if (ok) // the projection is good
        {
          // create a projected node and store it in snode
          MRTR::ProjectedNode* pnode = 
            new MRTR::ProjectedNode(*snode,xi,msegs[i],ssegs[j]->Id());
          snode->SetProjectedNode(pnode);
#if 0
          cout << " snode id: " << pnode->Id()
               << " projects on mseg: " << msegs[i]->Id()
               << " orth to sseg " << ssegs[j]->Id() << endl;
#endif          
        }
      } // for (int j=0; j<nsseg; ++j)
    } // for (int i=0; i<nmseg; ++i)
  } // for (scurr=rnode_[sside].begin(); scurr!=rnode_[sside].end(); ++scurr)  



  // loop all slave nodes again and make projections redundant
  if (lComm()->NumProc()>1)
  {
    vector<double> bcast(10*rnode_[sside].size());
    for (int proc=0; proc<lComm()->NumProc(); ++proc)
    {
      int blength=0;
      if (proc==lComm()->MyPID())
      {
        for (scurr=rnode_[sside].begin(); scurr!=rnode_[sside].end(); ++scurr)
        {
          MRTR::Node* snode = scurr->second;
          if (proc != NodePID(snode->Id())) continue; // cannot have a projection on a node i don't own
          int npnode=0;
          MRTR::ProjectedNode** pnode = snode->GetProjectedNode(npnode);        
          if (!pnode) continue; // no projection on this one
          bcast[blength] = (double)snode->Id();
          ++blength;
          bcast[blength] = (double)npnode;
          ++blength;
          for (int j=0; j<npnode; ++j)
          {
            bcast[blength] = (double)pnode[j]->Segment()->Id();
            ++blength;
            const double* xi = pnode[j]->Xi();
            bcast[blength] = xi[0];
            ++blength;
            bcast[blength] = xi[1];
            ++blength;
            bcast[blength] = pnode[j]->OrthoSegment();
            ++blength;
          }
          if (bcast.size() < blength+20) 
            bcast.resize(bcast.size()+40);
        } // for (mcurr=rnode_[mside].begin(); mcurr!=rnode_[mside].end(); ++mcurr)
        if (blength>=bcast.size())
        {
          cout << "***ERR*** MRTR::Interface::ProjectNodes_SlavetoMaster_Orthogonal:\n"
               << "***ERR*** Overflow in communication buffer occured\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          exit(EXIT_FAILURE);
        }
      } // if (proc==lComm()->MyPID())
      
      lComm()->Broadcast(&blength,1,proc);
      if (proc!=lComm()->MyPID()) bcast.resize(blength);
      lComm()->Broadcast(&bcast[0],blength,proc);
      
      if (proc!=lComm()->MyPID())
      {
        int i;
        for (i=0; i<blength;)
        {
          int nid    = (int)bcast[i] ; ++i;
          MRTR::Node* snode = GetNodeView(nid);
          int npnode = (int) bcast[i]; ++i;
          for (int j=0; j<npnode; ++j)
          {
            int sid     = (int)bcast[i]; ++i;
            double* xi  = &bcast[i];     ++i; ++i;
            int orthseg = (int)bcast[i]; ++i;
            MRTR::Segment* seg   = GetSegmentView(sid);
            MRTR::ProjectedNode* pnode = new MRTR::ProjectedNode(*snode,xi,seg,orthseg);
            snode->SetProjectedNode(pnode);
          }
        }
        if (i != blength)
        {
          cout << "***ERR*** MRTR::Interface::ProjectNodes_SlavetoMaster_Orthogonal:\n"
               << "***ERR*** Mismatch in dimension of recv buffer: " << i << " != " << blength << "\n"
               << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
          exit(EXIT_FAILURE);
        }
      } // if (proc!=lComm()->MyPID())
    } // for (int proc=0; proc<lComm()->NumProc(); ++proc)
    bcast.clear();
  } // if (lComm()->NumProc()>1)

#if 1
  // Postprocess the projections
  // The slave side of the interface might be larger then the master side
  // of the interface so not all slave nodes have a projection.
  // For those slave nodes without a projection attached to a slave segment
  // which overlaps with the master side, lagrange mutlipliers have to be
  // introduced. This is done by checking all nodes without a projection 
  // whether they are attached to some slave segment on which another node
  // HAS a projection. If this case is found, a pseudo ProjectedNode is 
  // introduced for that node.
  for (scurr=rnode_[sside].begin(); scurr!=rnode_[sside].end(); ++scurr)
  {
    MRTR::Node* snode = scurr->second;
    // do only my own nodes

    // don't do anything on nodes that already have a projection
    if (snode->GetProjectedNode())
      continue;

    // get segments adjacent to this node  
    int nseg             = snode->Nseg();
    MRTR::Segment** segs = snode->Segments();
    // loop segments and check for other nodes with projection
    bool foundit = false;
    for (int i=0; i<nseg; ++i)
    {
      int nnode = segs[i]->Nnode();
      MRTR::Node** nodes = segs[i]->Nodes();
      for (int j=0; j<nnode; ++j)
        if (nodes[j]->GetProjectedNode())
          if (nodes[j]->GetProjectedNode()->Segment())
          {
            foundit = true;
            break;
          }
      if (foundit) break;
    }
    if (foundit)
    {
#if 0
      cout << "Node without projection:\n" << *snode;        
      cout << "...get's lagrange multipliers\n\n";
#endif
      MRTR::ProjectedNode* pnode = new MRTR::ProjectedNode(*snode,NULL,NULL);
      snode->SetProjectedNode(pnode);
    }
  }
#endif


  return true; 
}
#endif // TRILINOS_PACKAGE
