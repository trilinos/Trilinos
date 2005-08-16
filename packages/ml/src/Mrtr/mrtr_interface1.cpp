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
#include "mrtr_utils.H"

/*----------------------------------------------------------------------*
 |  finalize construction of this interface                             |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::Complete()
{ 
  if (IsComplete())
  {
    if (OutLevel()>0)
      cout << "***WRN*** MRTR::Interface::InterfaceComplete:\n"
           << "***WRN*** InterfaceComplete() was called before, do nothing\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return true;
  }
  
  //-------------------------------------------------------------------
  // check for NULL entries in maps
  bool ok = true;
  for (int i=0; i<2; ++i)
  {
    map<int,MRTR::Node*>::const_iterator curr;
    for (curr=node_[i].begin(); curr!=node_[i].end(); ++curr)
    {
      MRTR::Node* node = curr->second;
      if (!node)
      {
        cout << "***ERR*** MRTR::Interface::Complete:\n"
             << "***ERR*** Interface # " << Id_ << ":\n"
             << "***ERR*** found NULL entry in map of nodes\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        ok = false;
      }
    }
  }
  for (int i=0; i<2; ++i)
  {
    map<int,MRTR::Segment*>::const_iterator curr;
    for (curr=seg_[i].begin(); curr!=seg_[i].end(); ++curr)
    {
      MRTR::Segment* seg = curr->second;
      if (!seg)
      {
        cout << "***ERR*** MRTR::Interface::Complete:\n"
             << "***ERR*** Interface # " << Id_ << ":\n"
             << "***ERR*** found NULL entry in map of segments\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        ok = false;
      }
    }
  }
  int lok = ok;
  int gok = 1;
  gcomm_.MinAll(&lok,&gok,1);  
  if (!gok) return false;
  
  //-------------------------------------------------------------------
  // check whether all nodes for segments are present
  // (take in account that node might be on different processor)
  // this test is expensive and does not scale. It is therefore only performed
  // when user requests a high output level
  if (OutLevel()>1)
  {
    for (int proc=0; proc<gcomm_.NumProc(); ++proc)
    {
      for (int side=0; side<2; ++side)
      {
        // create length of list of all nodes adjacent to segments on proc
        int sendsize =  0;
        if (proc==gcomm_.MyPID())
        {
          map<int,MRTR::Segment*>::const_iterator curr;
          for (curr=seg_[side].begin(); curr!=seg_[side].end(); ++curr)
            sendsize += curr->second->Nnode();
        }
        gcomm_.Broadcast(&sendsize,1,proc);
        
        // create list of all nodes adjacent to segments on proc
        int* ids = new int[sendsize];
        if (proc==gcomm_.MyPID())
        {
          map<int,MRTR::Segment*>::const_iterator curr;
          int counter=0;
          for (curr=seg_[side].begin(); curr!=seg_[side].end(); ++curr)
          {
            const int* segids = curr->second->NodeIds();
            for (int i=0; i<curr->second->Nnode(); ++i)
              ids[counter++] = segids[i];
          }
        }
        gcomm_.Broadcast(ids,sendsize,proc);
        
        // check on all processors for nodes in ids
        int* foundit  = new int[sendsize];
        int* gfoundit = new int[sendsize];
        for (int i=0; i<sendsize; ++i) 
        {
          foundit[i] = 0;
          if (node_[side].find(ids[i]) != node_[side].end()) 
            foundit[i] = 1;
        }
        gcomm_.MaxAll(foundit,gfoundit,sendsize);
        for (int i=0; i<sendsize; ++i)
        {
          if (gfoundit[i]!=1)
          {
            if (gcomm_.MyPID()==proc)
            cout << "***ERR*** MRTR::Interface::Complete:\n"
                 << "***ERR*** cannot find segment's node # " << ids[i] << "\n"
                 << "***ERR*** in map of all nodes on all procs\n"
                 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
            delete [] ids;
            delete [] foundit;
            delete [] gfoundit;
            gcomm_.Barrier();
            return false;
          }
        }
        
        // tidy up
        delete [] ids;
        delete [] foundit;
        delete [] gfoundit;
      } // for (int size=0; side<2; ++side)
    } // for (int proc=0; proc<gcomm_.NumProc(); ++proc)
  }
  
  //-------------------------------------------------------------------
  // find all procs that have business on this interface (own nodes/segments)
  // build a Epetra_comm that contains only those procs
  // this intra-communicator will be used to handle most stuff on this 
  // interface so the interface will not block all other procs
  {
#ifdef PARALLEL
    int* lin = new int[gcomm_.NumProc()];
    int* gin = new int[gcomm_.NumProc()];
    for (int i=0; i<gcomm_.NumProc(); ++i) lin[i] = 0;
    
    // check ownership of any segments
    for (int i=0; i<2; ++i)
      if (seg_[i].size() != 0)
      {
        lin[gcomm_.MyPID()] = 1;
        break;
      }
    // check ownership of any nodes
    for (int i=0; i<2; ++i)
      if (node_[i].size() != 0)
      {
        lin[gcomm_.MyPID()] = 1;
        break;
      }
    gcomm_.MaxAll(lin,gin,gcomm_.NumProc());
    delete [] lin; lin = NULL;
    
    // typecast the Epetra_Comm to Epetra_MpiComm
    Epetra_MpiComm* epetrampicomm = dynamic_cast<Epetra_MpiComm*>(&gcomm_);
    if (!epetrampicomm)
    {
      cout << "***ERR*** MRTR::Interface::Complete:\n"
           << "***ERR*** Interface " << Id() << ": Epetra_Comm is not an Epetra_MpiComm\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }

    // split the communicator into participating and none-participating procs
    int color;
    int key = gcomm_.MyPID();
    // I am taking part in the new comm if I have any ownership 
    if (gin[gcomm_.MyPID()]) 
      color = 0; 
    // I am not taking part in the new comm
    else                    
      color = MPI_UNDEFINED;
      
    // tidy up
    delete [] gin; gin = NULL;

    // create the local communicator   
    MPI_Comm  mpi_global_comm = epetrampicomm->GetMpiComm();
    MPI_Comm* mpi_local_comm  = new MPI_Comm();
    MPI_Comm_split(mpi_global_comm,color,key,mpi_local_comm);

    // create the new Epetra_MpiComm
    if (*mpi_local_comm == MPI_COMM_NULL)
      lcomm_ = NULL;
    else
      lcomm_ = new Epetra_MpiComm(*mpi_local_comm); // FIXME: who destroys the MPI_Comm inside?

#if 0
    // test this stuff on the mpi level
    int grank,lrank;
    MPI_Comm_rank(mpi_global_comm,&grank);
    if (*mpi_local_comm != MPI_COMM_NULL)
      MPI_Comm_rank(*mpi_local_comm,&lrank);
    else
      lrank = -1;
    for (int proc=0; proc<gcomm_.NumProc(); ++proc)
    {
      if (proc==gcomm_.MyPID())
      cout << "using mpi    comms: I am global rank " << grank << " and local rank " << lrank << endl;
      gcomm_.Barrier();
    }
    // test this stuff on the epetra level
    if (lComm())
      for (int proc=0; proc<lcomm_->NumProc(); ++proc)
      {
        if (proc==lcomm_->MyPID())
        cout << "using epetra comms: I am global rank " << gcomm_.MyPID() << " and local rank " << lcomm_->MyPID() << endl;
        lcomm_->Barrier();
      }
    gcomm_.Barrier();
#endif


    
#else  // the easy serial case
    Epetra_SerialComm* serialcomm = dynamic_cast<Epetra_SerialComm*>(&gcomm_);
    if (!serialcomm)
    {
      cout << "***ERR*** MRTR::Interface::Complete:\n"
           << "***ERR*** Interface " << Id() << ": Epetra_Comm is not an Epetra_SerialComm\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    lcomm_ = new Epetra_SerialComm(*serialcomm);
#endif // end of #ifdef PARALLEL    
  }
  
  //-------------------------------------------------------------------
  // create a map of all nodes to there PID (process id)
  if (lComm())
    for (int proc=0; proc<lcomm_->NumProc(); ++proc)
    {
      int lnnodes = 0;
      if (proc==lcomm_->MyPID())
        lnnodes = node_[0].size() + node_[1].size();
      lcomm_->Broadcast(&lnnodes,1,proc);
      int* ids = new int[lnnodes];
      if (proc==lcomm_->MyPID())
      {
        map<int,MRTR::Node*>::const_iterator curr;
        int counter=0;
        for (int side=0; side<2; ++side)
          for (curr=node_[side].begin(); curr!=node_[side].end(); ++curr)
            ids[counter++] = curr->first;
      }
      lcomm_->Broadcast(ids,lnnodes,proc);
      for (int i=0; i<lnnodes; ++i)
        nodePID_.insert(pair<int,int>(ids[i],proc));
      delete [] ids;
    }
  
  //-------------------------------------------------------------------
  // create a map of all segments to there PID (process id)
  if (lComm())
    for (int proc=0; proc<lcomm_->NumProc(); ++proc)
    {
      int lnsegs = 0;
      if (proc==lcomm_->MyPID())
        lnsegs = seg_[0].size() + seg_[1].size();
      lcomm_->Broadcast(&lnsegs,1,proc);
      int* ids = new int[lnsegs];
      if (proc==lcomm_->MyPID())
      {
        map<int,MRTR::Segment*>::const_iterator curr;
        int counter=0;
        for (int side=0; side<2; ++side)
          for (curr=seg_[side].begin(); curr!=seg_[side].end(); ++curr)
            ids[counter++] = curr->first;
      }
      lcomm_->Broadcast(ids,lnsegs,proc);
      for (int i=0; i<lnsegs; ++i)
        segPID_.insert(pair<int,int>(ids[i],proc));
      delete [] ids;
    }
  
  //-------------------------------------------------------------------
  // set isComplete_ flag
  // we set it here already as we will be using some methods that require it
  // from now on
  isComplete_ = true;
  
  //-------------------------------------------------------------------
  // make the nodes know there adjacent segments
  // find max number of nodes to a segment
  if (lComm())
  {
    int lmaxnnode = 0;
    int gmaxnnode = 0;
    for (int side=0; side<2; ++side)
    {
      map<int,MRTR::Segment*>::const_iterator scurr;
      for (scurr=seg_[side].begin(); scurr!=seg_[side].end(); ++scurr)
        if (lmaxnnode < scurr->second->Nnode())
          lmaxnnode = scurr->second->Nnode();
    }
    lcomm_->MaxAll(&lmaxnnode,&gmaxnnode,1);
    
    // loop all procs and broadcast their adjacency
    for (int proc=0; proc<lcomm_->NumProc(); ++proc)
    {
      // local number of segments
      int lnseg = 0;
      if (proc==lcomm_->MyPID())
        lnseg = seg_[0].size() + seg_[1].size();
      lcomm_->Broadcast(&lnseg,1,proc);
      
      // allocate vector to hold adjacency
      int offset = gmaxnnode+2;
      int size   = lnseg*offset;
      int* adj   = new int[size];
      
      // proc fills adjacency vector adj and broadcasts
      if (proc==lcomm_->MyPID())
      {
        int count = 0;
        for (int side=0; side<2; ++side)
        {
          map<int,MRTR::Segment*>::const_iterator scurr;
          for (scurr=seg_[side].begin(); scurr!=seg_[side].end(); ++scurr)
          {
            Segment* seg = scurr->second;
            adj[count]   = seg->Id();
            adj[count+1] = seg->Nnode();
            const int* ids = seg->NodeIds();
            for (int i=0; i<seg->Nnode(); ++i)
              adj[count+2+i] = ids[i];
            count += offset;
          }
        }
      }
      lcomm_->Broadcast(adj,size,proc);
      
      // all procs read adj and add segment to the nodes they own
      int count = 0;
      for (int i=0; i<lnseg; ++i)
      {
        int segid = adj[count];
        int nnode = adj[count+1];
        for (int j=0; j<nnode; ++j)
        {
          int nid = adj[count+2+j];
          if (lcomm_->MyPID() == NodePID(nid))
          {
            // I own this node, so set the segment segid in it
            MRTR::Node* node = GetNodeViewLocal(nid);
            if (!node)
            {
              cout << "***ERR*** MRTR::Interface::Complete:\n"
                   << "***ERR*** cannot find node " << nid << endl
                   << "***ERR*** in map of all nodes on this proc\n"
                   << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
              exit(EXIT_FAILURE);
            }
            node->AddSegment(segid);
          }
          else
            continue;
        }
        count += offset;
      }
      delete [] adj;
    } // for (int proc=0; proc<lcomm_->NumProc(); ++proc)
  } // if (lComm())
  
  //-------------------------------------------------------------------
  // build redundant segments and nodes
  if (lComm())
  {
    int ok = 0;
    ok += RedundantSegments(0);  
    ok += RedundantSegments(1);  
    ok += RedundantNodes(0);
    ok += RedundantNodes(1);
    if (ok != 4)
    {
      cout << "***ERR*** MRTR::Interface::Complete:\n"
           << "***ERR*** building of redundant information failed\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
  }

  //-------------------------------------------------------------------
  // make topology segments <-> nodes for each side
  if (lComm())
    BuildNodeSegmentTopology(); 

  //-------------------------------------------------------------------
  // delete distributed nodes and segments
  for (int i=0; i<2; ++i)
  {
    MRTR::DestroyMap(seg_[i]);
    MRTR::DestroyMap(node_[i]);
  }


  //-------------------------------------------------------------------
  // we are done
  // note that there might not be any functions on the interface yet
  // they still have to be set
  
  return ok;
}

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
      cout << "***ERR*** MRTR::Interface::Mortar_Integrate:\n"
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
        cout << "***ERR*** MRTR::Interface::Mortar_Integrate:\n"
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
      ncurr->second->BuildAveragedNormal();

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
  else
  {
    cout << "***ERR*** MRTR::Interface::Mortar_Integrate:\n"
         << "***ERR*** interface " << Id() << "\n"
         << "***ERR*** currently only projection type MRTR::Interface::proj_continousnormalfield\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  return true;
}



/*----------------------------------------------------------------------*
 |  make mortar integration of this interface                           |
 *----------------------------------------------------------------------*/
bool MRTR::Interface::Mortar_Integrate(Epetra_CrsMatrix& D, 
                                       Epetra_CrsMatrix& M)
{ 
  bool ok = false;
  
  //-------------------------------------------------------------------
  // interface needs to be complete
  if (!IsComplete())
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MRTR::Interface::Mortar_Integrate:\n"
           << "***ERR*** Complete() not called on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // send all procs not member of this interface's intra-comm out of here
  // FIXME for testing, leave them in
  //if (!lComm()) return true;

  //-------------------------------------------------------------------
  // interface needs to have a mortar side assigned
  if (MortarSide()==-1)
  {
    if (gcomm_.MyPID()==0)
      cout << "***ERR*** MRTR::Interface::Mortar_Integrate:\n"
           << "***ERR*** mortar side was not assigned on interface " << Id_ << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //-------------------------------------------------------------------
  // interface segments need to have at least one function on the mortar side
  // and two functions on the slave side
  int mside = MortarSide();
  int sside = OtherSide(mside);
  map<int,MRTR::Segment*>::iterator scurr;
  for (scurr=seg_[mside].begin(); scurr!=seg_[mside].end(); ++scurr)
    if (scurr->second->Nfunctions() < 1)
    {
      cout << "***ERR*** MRTR::Interface::Mortar_Integrate:\n"
           << "***ERR*** interface " << Id_ << ", mortar side\n"
           << "***ERR*** segment " << scurr->second->Id() << " needs at least 1 function set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  for (scurr=seg_[sside].begin(); scurr!=seg_[sside].end(); ++scurr)
    if (scurr->second->Nfunctions() < 2)
    {
      cout << "***ERR*** MRTR::Interface::Mortar_Integrate:\n"
           << "***ERR*** interface " << Id_ << ", slave side\n"
           << "***ERR*** segment " << scurr->second->Id() << " needs at least 2 function set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    
  //-------------------------------------------------------------------
  // do the integration of the master side
  if (IsOneDimensional())
  {
    ok = Integrate_MasterSide_2D(M);
    if (!ok) return false;
  }
  else
  {
    cout << "***ERR*** MRTR::Interface::Mortar_Integrate:\n"
         << "***ERR*** Interface " << Id() << " 2D interface integration not yet impl.\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }


  //-------------------------------------------------------------------
  // do the integration of the slave side
  if (IsOneDimensional())
  {
    ok = Integrate_SlaveSide_2D(D);
    if (!ok) return false;
  }
  else
  {
    cout << "***ERR*** MRTR::Interface::Mortar_Integrate:\n"
         << "***ERR*** Interface " << Id() << " 2D interface integraion not yet impl.\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }

  //-------------------------------------------------------------------
  // set the flag that this interface has been successfully integrated
  isIntegrated_ = true;
  
  return true;
}


#endif // TRILINOS_PACKAGE
