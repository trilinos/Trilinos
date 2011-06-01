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
#include "mrtr_utils.H"

/*----------------------------------------------------------------------*
 |  finalize construction of this interface                             |
 *----------------------------------------------------------------------*/
bool MOERTEL::Interface::Complete()
{ 
  if (IsComplete())
  {
    if (OutLevel()>0)
      cout << "MOERTEL: ***WRN*** MOERTEL::Interface::InterfaceComplete:\n"
           << "MOERTEL: ***WRN*** InterfaceComplete() was called before, do nothing\n"
           << "MOERTEL: ***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return true;
  }
  
  //-------------------------------------------------------------------
  // check for NULL entries in maps
  bool ok = true;
  for (int i=0; i<2; ++i)
  {
	std::map<int,Teuchos::RCP<MOERTEL::Node> >::const_iterator curr;
    for (curr=node_[i].begin(); curr!=node_[i].end(); ++curr)
    {
      if (curr->second == Teuchos::null)
      {
        cout << "***ERR*** MOERTEL::Interface::Complete:\n"
             << "***ERR*** Interface # " << Id_ << ":\n"
             << "***ERR*** found NULL entry in map of nodes\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        ok = false;
      }
    }
  }
  for (int i=0; i<2; ++i)
  {
	std::map<int,Teuchos::RCP<MOERTEL::Segment> >::const_iterator curr;
    for (curr=seg_[i].begin(); curr!=seg_[i].end(); ++curr)
    {
      if (curr->second == Teuchos::null)
      {
        cout << "***ERR*** MOERTEL::Interface::Complete:\n"
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
#if 1
  if (OutLevel()>9)
  {
    for (int proc=0; proc<gcomm_.NumProc(); ++proc)
    {
      for (int side=0; side<2; ++side)
      {
        // create length of list of all nodes adjacent to segments on proc
        int sendsize =  0;
        if (proc==gcomm_.MyPID())
        {
		  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::const_iterator curr;
          for (curr=seg_[side].begin(); curr!=seg_[side].end(); ++curr)
            sendsize += curr->second->Nnode();
        }
        gcomm_.Broadcast(&sendsize,1,proc);
        
        // create list of all nodes adjacent to segments on proc
		std::vector<int> ids(sendsize);
        if (proc==gcomm_.MyPID())
        {
		  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::const_iterator curr;
          int counter=0;
          for (curr=seg_[side].begin(); curr!=seg_[side].end(); ++curr)
          {
            const int* segids = curr->second->NodeIds();
            for (int i=0; i<curr->second->Nnode(); ++i)
              ids[counter++] = segids[i];
          }
        }
        gcomm_.Broadcast(&ids[0],sendsize,proc);
        
        // check on all processors for nodes in ids
		std::vector<int> foundit(sendsize);
		std::vector<int> gfoundit(sendsize);
        for (int i=0; i<sendsize; ++i) 
        {
          foundit[i] = 0;
          if (node_[side].find(ids[i]) != node_[side].end()) 
            foundit[i] = 1;
        }
        gcomm_.MaxAll(&foundit[0],&gfoundit[0],sendsize);
        for (int i=0; i<sendsize; ++i)
        {
          if (gfoundit[i]!=1)
          {
            if (gcomm_.MyPID()==proc)
            cout << "***ERR*** MOERTEL::Interface::Complete:\n"
                 << "***ERR*** cannot find segment's node # " << ids[i] << "\n"
                 << "***ERR*** in map of all nodes on all procs\n"
                 << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
            ids.clear();
            foundit.clear();
            gfoundit.clear();
            gcomm_.Barrier();
            return false;
          }
        }
        
        // tidy up
        ids.clear();
        foundit.clear();
        gfoundit.clear();
      } // for (int size=0; side<2; ++side)
    } // for (int proc=0; proc<gcomm_.NumProc(); ++proc)
  }
#endif  
  //-------------------------------------------------------------------
  // find all procs that have business on this interface (own nodes/segments)
  // build a Epetra_comm that contains only those procs
  // this intra-communicator will be used to handle most stuff on this 
  // interface so the interface will not block all other procs
  {
#ifdef EPETRA_MPI
	std::vector<int> lin(gcomm_.NumProc());
	std::vector<int> gin(gcomm_.NumProc());
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
    gcomm_.MaxAll(&lin[0],&gin[0],gcomm_.NumProc());
    lin.clear();
    
    // typecast the Epetra_Comm to Epetra_MpiComm
    Epetra_MpiComm* epetrampicomm = dynamic_cast<Epetra_MpiComm*>(&gcomm_);
    if (!epetrampicomm)
    {
	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Interface::Complete:\n"
           << "***ERR*** Interface " << Id() << ": Epetra_Comm is not an Epetra_MpiComm\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      throw ReportError(oss);
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
    gin.clear();

    // create the local communicator   
    MPI_Comm  mpi_global_comm = epetrampicomm->GetMpiComm();
    MPI_Comm* mpi_local_comm  = new MPI_Comm();
    MPI_Comm_split(mpi_global_comm,color,key,mpi_local_comm);

    // create the new Epetra_MpiComm
    if (*mpi_local_comm == MPI_COMM_NULL)
      lcomm_ = Teuchos::null;
    else
      lcomm_ = Teuchos::rcp(new Epetra_MpiComm(*mpi_local_comm)); // FIXME: who destroys the MPI_Comm inside?

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
	  std::stringstream oss;
			oss << "***ERR*** MOERTEL::Interface::Complete:\n"
           << "***ERR*** Interface " << Id() << ": Epetra_Comm is not an Epetra_SerialComm\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      throw ReportError(oss);
    }
    lcomm_ = Teuchos::rcp(new Epetra_SerialComm(*serialcomm));
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
	  std::vector<int> ids(lnnodes);
      if (proc==lcomm_->MyPID())
      {
		std::map<int,Teuchos::RCP<MOERTEL::Node> >::const_iterator curr;
        int counter=0;
        for (int side=0; side<2; ++side)
          for (curr=node_[side].begin(); curr!=node_[side].end(); ++curr)
            ids[counter++] = curr->first;
      }
      lcomm_->Broadcast(&ids[0],lnnodes,proc);
      for (int i=0; i<lnnodes; ++i)
        nodePID_.insert(std::pair<int,int>(ids[i],proc));
      ids.clear();
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
	  std::vector<int> ids(lnsegs);
      if (proc==lcomm_->MyPID())
      {
		std::map<int,Teuchos::RCP<MOERTEL::Segment> >::const_iterator curr;
        int counter=0;
        for (int side=0; side<2; ++side)
          for (curr=seg_[side].begin(); curr!=seg_[side].end(); ++curr)
            ids[counter++] = curr->first;
      }
      lcomm_->Broadcast(&ids[0],lnsegs,proc);
      for (int i=0; i<lnsegs; ++i)
        segPID_.insert(std::pair<int,int>(ids[i],proc));
      ids.clear();
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
	  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::const_iterator scurr;
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
	  std::vector<int> adj(size);
      
      // proc fills adjacency vector adj and broadcasts
      if (proc==lcomm_->MyPID())
      {
        int count = 0;
        for (int side=0; side<2; ++side)
        {
		  std::map<int,Teuchos::RCP<MOERTEL::Segment> >::const_iterator scurr;
          for (scurr=seg_[side].begin(); scurr!=seg_[side].end(); ++scurr)
          {
			Teuchos::RCP<MOERTEL::Segment> seg = scurr->second;
            adj[count]   = seg->Id();
            adj[count+1] = seg->Nnode();
            const int* ids = seg->NodeIds();
            for (int i=0; i<seg->Nnode(); ++i)
              adj[count+2+i] = ids[i];
            count += offset;
          }
        }
      }
      lcomm_->Broadcast(&adj[0],size,proc);
      
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
			Teuchos::RCP<MOERTEL::Node> node = GetNodeViewLocal(nid);
            if (node == Teuchos::null)
            {
				std::stringstream oss;
					oss << "***ERR*** MOERTEL::Interface::Complete:\n"
						<< "***ERR*** cannot find node " << nid << "\n"
                   << "***ERR*** in map of all nodes on this proc\n"
                   << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
				throw ReportError(oss);
            }
            node->AddSegment(segid);
          }
          else
            continue;
        }
        count += offset;
      }
      adj.clear();
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
		std::stringstream oss;
			oss << "***ERR*** MOERTEL::Interface::Complete:\n"
           << "***ERR*** building of redundant information failed\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
		throw ReportError(oss);
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
    seg_[i].clear();
    node_[i].clear();
  }


  //-------------------------------------------------------------------
  // we are done
  // note that there might not be any functions on the interface yet
  // they still have to be set
  
  return ok;
}
