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

#include "mrtr_manager.H"
#include "EpetraExt_MatrixMatrix.h"  // for adding matrices
#include <EpetraExt_Transpose_RowMatrix.h>

/*----------------------------------------------------------------------*
 |  ctor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
MRTR::Manager::Manager(Epetra_Comm& comm, int outlevel) :
outlevel_(outlevel),
comm_(comm),
inputmap_(NULL),
inputmatrixisdeep_(false),
inputmatrix_(NULL),
constraintsmap_(NULL),
D_(NULL),
M_(NULL),
saddlemap_(NULL),
saddlematrix_(NULL)
{
}

/*----------------------------------------------------------------------*
 |  dtor (public)                                            mwgee 06/05|
 *----------------------------------------------------------------------*/
MRTR::Manager::~Manager()
{
  // delete interfaces
  map<int,MRTR::Interface*>::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    if (curr->second)
    {
      delete curr->second;
      curr->second = NULL;
    }
    else if (OutLevel()>0)
      cout << "***WRN*** MRTR::Manager::~Manager:\n"
           << "***WRN*** found NULL entry in map of interfaces\n"
           << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
  }
  interface_.clear();
  
  if (inputmap_) delete inputmap_; inputmap_ = NULL;
  if (inputmatrix_) 
  {
    if (inputmatrixisdeep_)
      delete inputmatrix_; 
    inputmatrix_ = NULL;
    inputmatrixisdeep_ = false;
  }
  if (constraintsmap_) delete constraintsmap_; constraintsmap_ = NULL;
  if (D_) delete D_; D_ = NULL;
  if (M_) delete M_; M_ = NULL;
  if (saddlemap_) delete saddlemap_; saddlemap_ = NULL;
  if (saddlematrix_) delete saddlematrix_; saddlematrix_ = NULL;
}

/*----------------------------------------------------------------------*
 |  << operator                                              mwgee 06/05|
 *----------------------------------------------------------------------*/
ostream& operator << (ostream& os, const MRTR::Manager& man)
{ 
  man.Print();
  return (os);
}

/*----------------------------------------------------------------------*
 |  print all data                                           mwgee 06/05|
 *----------------------------------------------------------------------*/
bool MRTR::Manager::Print() const
{ 
  comm_.Barrier();

  if (Comm().MyPID()==0)
  cout << "\n========================= Mortar Manager =========================\n\n";
  fflush(stdout);
  
  comm_.Barrier();

  map<int,MRTR::Interface*>::const_iterator curr;
  for (curr=interface_.begin(); curr!=interface_.end(); ++curr)
  {
    MRTR::Interface* inter = curr->second;
    if (!inter)
    {
      cout << "***ERR*** MRTR::Manager::Print:\n"
           << "***ERR*** found NULL entry in map of interfaces\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    cout << *inter;
  }
  comm_.Barrier();
  if (inputmap_)
  {
    if (comm_.MyPID() == 0)
    cout << "\n------------------ Input RowMap of original Problem ---------------\n";
    comm_.Barrier();
    cout << *inputmap_;
  }
  comm_.Barrier();
  if (constraintsmap_)
  {
    if (comm_.MyPID() == 0)
    cout << "\n------------------ RowMap of Constraints ---------------\n";
    comm_.Barrier();
    cout << *constraintsmap_;
  }
  comm_.Barrier();
  if (D_)
  {
    if (comm_.MyPID() == 0)
    cout << "\n------------------ Coupling Matrix D ---------------\n";
    comm_.Barrier();
    cout << *D_;
  }
  comm_.Barrier();
  if (M_)
  {
    if (comm_.MyPID() == 0)
    cout << "\n------------------ Coupling Matrix M ---------------\n";
    comm_.Barrier();
    cout << *M_;
  }
  comm_.Barrier();

  if (Comm().MyPID()==0)
  cout << "\n========================= End Mortar Manager =========================\n\n";
  fflush(stdout);

  return true;
}

/*----------------------------------------------------------------------*
 |  Add an interface (public)                                mwgee 06/05|
 *----------------------------------------------------------------------*/
bool MRTR::Manager::AddInterface(MRTR::Interface& interface) 
{
  if (!interface.IsComplete())
  {
    cout << "***ERR*** MRTR::Manager::AddInterface:\n"
         << "***ERR*** Cannot add segment as Complete() was not called before\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  MRTR::Interface* tmp = new MRTR::Interface(interface);
  interface_.insert(pair<int,MRTR::Interface*>(tmp->Id(),tmp));
  
  return true;
}

/*----------------------------------------------------------------------*
 |  Store the rowmap of the underlying matrix from the application 07/05|
 *----------------------------------------------------------------------*/
bool MRTR::Manager::SetInputMap(Epetra_Map* map)
{
  if (inputmap_)
    delete inputmap_;
  inputmap_ = new Epetra_Map(*map);
  return true;
}

/*----------------------------------------------------------------------*
 |  Set ptr to the uncoupled matrix on input                       07/05|
 |  Note that MRTR::Manager does not make a deep copy here by default   |
 *----------------------------------------------------------------------*/
bool MRTR::Manager::SetInputMatrix(Epetra_CrsMatrix* inputmatrix, bool DeepCopy)
{
  if (DeepCopy)
  {
    inputmatrix_ = new Epetra_CrsMatrix(*inputmatrix);
    inputmatrixisdeep_ = true;
    return true;
  }
  else
  {
    inputmatrix_ = inputmatrix;
    inputmatrixisdeep_ = false;
    return true;
  }
  return false;
}

/*----------------------------------------------------------------------*
 |                                                                 08/05|
 |  modified version of the epetraext matrixmatrixadd                   |
 |  NOTE:                                                               |
 |  - A has to be FillComplete, B must NOT be FillComplete()            |
 *----------------------------------------------------------------------*/
int MRTR::Manager::MatrixMatrixAdd(const Epetra_CrsMatrix& A, bool transposeA,double scalarA,
                      Epetra_CrsMatrix& B,double scalarB )
{
  //
  //This method forms the matrix-matrix sum B = scalarA * op(A) + scalarB * B, where

  if (!A.Filled())
  {
     cout << "***ERR*** MRTR::Manager::MatrixMatrixAdd:\n"
          << "***ERR*** FillComplete was not called on A\n"
          << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
     exit(EXIT_FAILURE);
  }
  if (B.Filled())
  {
     cout << "***ERR*** MRTR::Manager::MatrixMatrixAdd:\n"
          << "***ERR*** FillComplete was called on B\n"
          << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
     exit(EXIT_FAILURE);
  }
  
  //explicit tranpose A formed as necessary
  Epetra_CrsMatrix* Aprime = 0;
  EpetraExt::RowMatrix_Transpose* Atrans = 0;
  if( transposeA )
  {
    Atrans = new EpetraExt::RowMatrix_Transpose(false);
    Aprime = &(dynamic_cast<Epetra_CrsMatrix&>(((*Atrans)(const_cast<Epetra_CrsMatrix&>(A)))));
  }
  else
    Aprime = const_cast<Epetra_CrsMatrix*>(&A);
    
  //Loop over B's rows and sum into
  int MaxNumEntries = EPETRA_MAX( Aprime->MaxNumEntries(), B.MaxNumEntries() );
  int NumEntries;
  int * Indices = new int[MaxNumEntries];
  double * Values = new double[MaxNumEntries];

  int NumMyRows = A.NumMyRows();
  int Row, err;

  if( scalarA )
  {
    for( int i = 0; i < NumMyRows; ++i )
    {
      Row = Aprime->GRID(i);
      EPETRA_CHK_ERR(Aprime->ExtractGlobalRowCopy(Row,MaxNumEntries,NumEntries,Values,Indices));
      if( scalarA != 1.0 )
        for( int j = 0; j < NumEntries; ++j ) Values[j] *= scalarA;
      if( B.Filled() ) {//Sum In Values
        err = B.SumIntoGlobalValues( Row, NumEntries, Values, Indices );
        assert( err == 0 );
      }
      else {
        err = B.InsertGlobalValues( Row, NumEntries, Values, Indices );
        assert( err == 0 || err == 1 );
      }
    }
  }

  delete [] Indices;
  delete [] Values;

  if( Atrans ) delete Atrans;

  return(0);
}

/*----------------------------------------------------------------------*
 |  Choose dofs for lagrange multipliers (private)           mwgee 07/05|
 | Note that this is collective for ALL procs                           |
 *----------------------------------------------------------------------*/
Epetra_Map* MRTR::Manager::LagrangeMultiplierDofs()
{
  if (!inputmap_)
  {
    cout << "***ERR*** MRTR::Manager::LagrangeMultiplierDofs:\n"
         << "***ERR*** inputmap==NULL, Need to set an input-rowmap first\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  // find the largest row number in inputmap
  int maxinputGID = inputmap_->MaxAllGID();

  // start numbering lagrange multipliers from maxinputGID+1
  int minLMGID = maxinputGID+1;
  int maxLMGID = 0;
  
  int length = 0;
  
  // loop interfaces and set LM dofs on the slave side for all slave nodes
  // that have a projection
  // start with minLMGID and return maxLMGID+1 on a specific interface
  // Note this is collective for ALL procs
  map<int,MRTR::Interface*>::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    length -= minLMGID;
    maxLMGID = curr->second->SetLMDofs(minLMGID);
    if (!maxLMGID && maxLMGID!=minLMGID)
    {
      cout << "***ERR*** MRTR::Manager::LagrangeMultiplierDofs:\n"
           << "***ERR*** interface " << curr->second->Id() << " returned false\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    minLMGID = maxLMGID;
    length += maxLMGID;
  }  

  vector<int> mylmids(length);
  int count=0;
  
  // get a vector of LM dof ids from each proc and each interface
  // and add it to the global one
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    MRTR::Interface* inter = curr->second;
    vector<int>* lmids = inter->MyLMIds();
    if (count+lmids->size() > mylmids.size())
      mylmids.resize(mylmids.size()+5*lmids->size());
    for (int i=0; i<lmids->size(); ++i)
      mylmids[count++] = (*lmids)[i];
  }  
  mylmids.resize(count);
  int lsize = count;
  int gsize = 0;
  comm_.SumAll(&lsize,&gsize,1);
  
  // create the rowmap for the constraints
  // Note that this map contains the global communicator from the MRTR::Manager
  // NOT any interface local one
  Epetra_Map* map = new Epetra_Map(gsize,lsize,&(mylmids[0]),0,comm_);

  // tidy up
  mylmids.clear();
  
  return map;
}






/*----------------------------------------------------------------------*
 |  integrate all mortar interfaces                       (public) 07/05|
 |  Note: All processors have to go in here!                            |
 *----------------------------------------------------------------------*/
bool MRTR::Manager::Mortar_Integrate()
{
  //-------------------------------------------------------------------
  // check whether we have interfaces
  int ninter = Ninterfaces();
  if (!ninter)
    return true;

  //-------------------------------------------------------------------
  // check whether we have an input map
  if (!inputmap_)
  {
    cout << "***ERR*** MRTR::Manager::Mortar_Integrate:\n"
         << "***ERR*** inputmap==NULL, Need to set an input-rowmap first\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    exit(EXIT_FAILURE);
  }
  
  //-------------------------------------------------------------------
  // check whether we have a mortar side chosen on each interface or 
  // whether we have to chose it automatically
  {
    map<int,MRTR::Interface*>::iterator curr;
    bool foundit = true;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      int mside = curr->second->MortarSide();
      if (mside==-2)
      {
        foundit = false;
        break;
      }
    }
    if (!foundit); // we have to chose mortar sides ourself
      ChooseMortarSide();
  }  

  //-------------------------------------------------------------------
  // build projections for all interfaces
  {
    map<int,MRTR::Interface*>::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      bool ok  = curr->second->Project();
      if (!ok)
      {
        cout << "***ERR*** MRTR::Manager::Mortar_Integrate:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false on projection\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
    }
  }  

  //-------------------------------------------------------------------
  // this is probably the place to put detection of end segments
  // for each end segment, the order of the lagrange multiplier shape
  // function will be reduced by one
#if 0
  {
    map<int,MRTR::Interface*>::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {
      bool ok = curr->second->DetectEndSegmentsandReduceOrder();
      if (!ok)
      {
        cout << "***ERR*** MRTR::Manager::Mortar_Integrate:\n"
             << "***ERR*** interface " << curr->second->Id() << " returned false from DetectEndSegmentsandReduceOrder\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
    }
  }
#endif
  //-------------------------------------------------------------------
  // choose dofs for lagrange multipliers and set them to slave nodes
  // build the rowmap for the coupling matrices M and D
  {
    if (constraintsmap_) delete constraintsmap_;
    constraintsmap_ = LagrangeMultiplierDofs();
    if (!constraintsmap_)
    {
      cout << "***ERR*** MRTR::Manager::Mortar_Integrate:\n"
           << "***ERR*** LagrangeMultiplierDofs() returned NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
  }
    
  //-------------------------------------------------------------------
  // build the map for the saddle point problem
  {
    // the saddle point problem rowmap is the inputmap + the constraintmap
    int numglobalelements = inputmap_->NumGlobalElements() +
                            constraintsmap_->NumGlobalElements();
    int nummyelements     = inputmap_->NumMyElements() +
                            constraintsmap_->NumMyElements();
    vector<int> myglobalelements(nummyelements);
    int count = 0;
    int* inputmyglobalelements = inputmap_->MyGlobalElements();
    for (int i=0; i<inputmap_->NumMyElements(); ++i)
      myglobalelements[count++] = inputmyglobalelements[i];
    int* constraintsmyglobalelements = constraintsmap_->MyGlobalElements();
    for (int i=0; i<constraintsmap_->NumMyElements(); ++i)
      myglobalelements[count++] = constraintsmyglobalelements[i];
    if (count != nummyelements)
    {
        cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
             << "***ERR*** Mismatch in dimensions\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
    }
    if (saddlemap_) delete saddlemap_;
    saddlemap_ = new Epetra_Map(numglobalelements,nummyelements,
                                &(myglobalelements[0]),0,comm_);
    myglobalelements.clear();
  }

  //-------------------------------------------------------------------
  // build the Epetra_CrsMatrix D and M
  {
    if (D_) delete D_;
    if (M_) delete M_;
    D_ = new Epetra_CrsMatrix(Copy,*saddlemap_,5,false);
    M_ = new Epetra_CrsMatrix(Copy,*saddlemap_,40,false);
  }


  //-------------------------------------------------------------------
  // integrate all interfaces
  {
    map<int,MRTR::Interface*>::iterator curr;
    for (curr=interface_.begin(); curr != interface_.end(); ++curr)
    {  
      MRTR::Interface* inter = curr->second;
      bool ok = inter->Mortar_Integrate(*D_,*M_);
      if (!ok)
      {
        cout << "***ERR*** MRTR::Manager::Mortar_Integrate:\n"
             << "***ERR*** interface " << inter->Id() << " returned false on integration\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
    }
  }
  
  //-------------------------------------------------------------------
  // call FillComplete() on M_ and D_ 
  D_->FillComplete(*saddlemap_,*saddlemap_);
  M_->FillComplete(*saddlemap_,*saddlemap_);
  
  return true;
}

/*----------------------------------------------------------------------*
 | Create the saddle point problem (public)                  mwgee 07/05|
 | Note that this is collective for ALL procs                           |
 *----------------------------------------------------------------------*/
Epetra_CrsMatrix* MRTR::Manager::MakeSaddleProblem()
{
  // check whether all interfaces are complete and integrated
  map<int,MRTR::Interface*>::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    if (curr->second->IsComplete() == false)
    {
      cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not Complete()\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    if (curr->second->IsIntegrated() == false)
    {
      cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not integrated yet\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
  }
  
  // check whether we have an inputmap
  if (!inputmap_)
  {
      cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
           << "***ERR*** No inputrowmap set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
  }
  
  // check whether we have a constraintsmap_
  if (!constraintsmap_)
  {
      cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
           << "***ERR*** onstraintsmap is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
  }
  
  // check for saddlemap_
  if (!saddlemap_)
  {
      cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
           << "***ERR*** saddlemap_==NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
  }

  // check for inputmatrix
  if (!inputmatrix_)
  {
      cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
           << "***ERR*** No inputmatrix set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
  }

  // check whether we have M and D matrices
  if (!D_ || !M_)
  {
      cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
           << "***ERR*** Matrix M or D is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
  }
  
  // create a matrix for the saddle problem and fill it
  if (saddlematrix_) delete saddlematrix_; 
  saddlematrix_ = new Epetra_CrsMatrix(Copy,*saddlemap_,90);

  // add values from inputmatrix
  MatrixMatrixAdd(*inputmatrix_,false,1.0,*saddlematrix_,0.0);
  
  // add values from D_
  MatrixMatrixAdd(*D_,false,1.0,*saddlematrix_,1.0);
  MatrixMatrixAdd(*D_,true,1.0,*saddlematrix_,1.0);  

  // add values from M_
  MatrixMatrixAdd(*M_,false,1.0,*saddlematrix_,1.0);
  MatrixMatrixAdd(*M_,true,1.0,*saddlematrix_,1.0);  

  saddlematrix_->FillComplete();

  return saddlematrix_;
}

/*----------------------------------------------------------------------*
 |                                                                 09/05|
 |  choose the mortar side                                              |
 *----------------------------------------------------------------------*/
bool MRTR::Manager::ChooseMortarSide()
{
  // find all 2D interfaces
  vector<MRTR::Interface*> inter(Ninterfaces());
  map<int,MRTR::Interface*>::iterator curr;
  curr=interface_.begin();
  int count = 0;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    if (curr->second->IsOneDimensional())
    {
      inter[count] = curr->second;
      ++count;
    }
  }
  inter.resize(count);
  
  // call choice of the mortar side for all 1D interfaces
  bool ok = ChooseMortarSide_2D(inter);

  // for 2D interfaces not yet impl.

  // tidy up
  inter.clear();  

  return ok;
}


/*----------------------------------------------------------------------*
 |                                                                 09/05|
 |  choose the mortar side for 1D interfaces                            |
 *----------------------------------------------------------------------*/
bool MRTR::Manager::ChooseMortarSide_2D(vector<MRTR::Interface*> inter)
{
  // number of interfaces
  int ninter = inter.size();
  if (ninter < 2) return true;
  
  if (OutLevel()>5)
  {
    fflush(stdout);
    if (Comm().MyPID()==0)
      cout << "---MRTR::Manager: finding common nodes on interfaces\n";
    fflush(stdout);
  }
  
  // get a view of all nodes from all interfaces
  vector<MRTR::Node**> nodes(ninter);
  vector<int>          nnodes(ninter);
  for (int i=0; i<ninter; ++i)
  {
    nnodes[i] = inter[i]->GlobalNnode(); // returns 0 for procs not in lComm()
    nodes[i]  = inter[i]->GetNodeView(); // get vector of ALL nodes on inter i
  }
  
  // loop all interfaces
  for (int i=0; i<ninter; ++i)
  {
    // loop all nodes on that interface 
    // (procs not part of inter[i]->lComm() don't loop because they have nnodes[i]=0)
    for (int j=0; j<nnodes[i]; ++j)
    {
      // do nothing for node that does not belong to me
      if (inter[i]->NodePID(nodes[i][j]->Id()) != inter[i]->lComm()->MyPID())
        continue;
        
      // do nothing for a node that has been flagged cornernode before
      if (nodes[i][j]->IsCorner())
        continue;
        
      // the node we are currently looking for
      int actnodeid       = nodes[i][j]->Id();
      MRTR::Node* inode   = nodes[i][j];
      
      // search all other interfaces for this node
      for (int k=0; k<ninter; ++k)
      {
        if (k==i) 
          continue; // don't do the current interface
        
        MRTR::Node* knode = NULL;
        for (int l=0; l<nnodes[k]; ++l)
          if (actnodeid==nodes[k][l]->Id())
          {
            knode = nodes[k][l];
            break;
          }
        if (!knode) // node inode is not on interface k 
          continue;
          
        // found node actnodeid on interface i and k
        if (OutLevel()>5)
        {
          cout << "Node " << actnodeid << " on interfaces " << inter[i]->Id() << " and " << inter[k]->Id() << endl;
          fflush(stdout);
        }
        
        // flag that node on interfaces i and k as cornernode
        inode->SetCorner();
        knode->SetCorner();
      } // for (int k=0; k<ninter; ++k)
    } // for (int j=0; j< nnodes[i]; ++j)
  } // for (int i=0; i<ninter; ++i)

  
  // make the cornernode information redundant
    // loop all interfaces
  for (int i=0; i<ninter; ++i)
  {
    // if I'm not part of inter[i], continue
    if (!inter[i]->lComm())
      continue;
    for (int proc=0; proc<inter[i]->lComm()->NumProc(); ++proc)
    {
      // get the corner nodes I have found
      vector<int> bcast(4);
      int bsize = 0;
      if (proc==inter[i]->lComm()->MyPID())
        for (int j=0; j<nnodes[i]; ++j)
        {
          if (nodes[i][j]->IsCorner())
          {
            if (bcast.size() <= bsize)
              bcast.resize(bsize+10);
            bcast[bsize] = j;
            ++bsize;  
          }
        } 
      inter[i]->lComm()->Broadcast(&bsize,1,proc);
      bcast.resize(bsize);
      inter[i]->lComm()->Broadcast(&(bcast[0]),bsize,proc);
      if (proc!=inter[i]->lComm()->MyPID())
        for (int k=0; k<bsize; ++k)
        {
          int j = bcast[k];
          nodes[i][j]->SetCorner();
        }
      bcast.clear();
    } // for (int proc=0; proc<Comm()->NumProc(); ++proc)
  } // for (i=0; i<ninter; ++i)

  
  // we have list of interfaces inter[i], now we need
  // to build list of cornernodes associated with each interface
  vector< vector<int> > cornernodes(ninter);
  for (int i=0; i<ninter; ++i) 
  {
    int bprocs = 0; int bprocr = 0;
    int bsize = 0;
    if (inter[i]->lComm())
      if (inter[i]->lComm()->MyPID()==0)
      {
        bprocs = Comm().MyPID();
        for (int j=0; j<nnodes[i]; ++j)
          if (nodes[i][j]->IsCorner())
          {
            if (cornernodes[i].size() <= bsize)
              cornernodes[i].resize(bsize+5);
            cornernodes[i][bsize] = nodes[i][j]->Id();
            ++bsize;
          }
      }
    Comm().MaxAll(&bprocs,&bprocr,1);
    Comm().Broadcast(&bsize,1,bprocr);
    cornernodes[i].resize(bsize);
    Comm().Broadcast(&(cornernodes[i][0]),bsize,bprocr);
  }
  
  // we now color the interfaces in such a way that every corner node
  // has lagrange multipliers either never or once but never more then once
  vector<int> flags((2*ninter));
  vector<int> flagr((2*ninter));
  for (int i=0; i<ninter; ++i)
  {
     // check whether inter[i] has mortar side assigned
     // if not, go through all cornernodes and check dependency on
     // other interfaces
     // if there is a dependency, choose mortar side of inter[i] approp., 
     // if there is no dependency, choose something

     // create a vector of flags as follows:
     // flags[0..ninter-1]        flags for side 0 of inter[i]
     // flags[ninter..2*ninter-1] flags for side 1 of inter[i]
     // flags[k]==  0 : no dependency
     // flags[k]==  1 : dependency but no side chosen yet on inter[k]
     // flags[k]==  2 : side should be chosen master/mortar side
     // flags[k]==  3 : side should be chosen slave/LM side
     for (int k=0; k<2*ninter; ++k) flags[k] = 0;
       
     if (inter[i]->MortarSide() == -2)
     {
       // actnode[0] holds the active node id
       // actnode[1] holds the side this node is on on inter[i]  
       int actnodes[2]; 
       int actnoder[2]; 
       // loop cornernodes on inter[i]
       for (int j=0; j<cornernodes[i].size(); ++j)
       {
         actnodes[0] = 0; actnodes[1] = 0;    
         if (inter[i]->lComm())
           if (inter[i]->lComm()->MyPID()==0)
           {
             actnodes[0] = cornernodes[i][j];
             actnodes[1] = inter[i]->GetSide(actnodes[0]);
           }
         Comm().MaxAll(actnodes,actnoder,2);
         
         // loop all other interfaces and check for actnoder[0]
         for (int k=0; k<ninter; ++k)
         {
           if (k==i) continue; // don't do inter[i]
           if (inter[k]->lComm())
             if (inter[k]->lComm()->MyPID()==0)
             {
               for (int l=0; l<cornernodes[k].size(); ++l)
               {
                 if (actnoder[0]==cornernodes[k][l])
                 {
                   flags[k]        = 1; // found a dependency
                   flags[ninter+k] = 1; // found a dependancy
                   int mside = inter[k]->MortarSide();
                   if (mside != -2) // side has been chosen before on inter[k]
                   {
                     int nodeside = inter[k]->GetSide(actnoder[0]);
                     // found node on master/mortar side
                     // so set flags such that node is on slave on inter[i]
                     if (nodeside == mside) 
                     {
                       // actnode[1] has to become slave side on inter[i]
                       int sside_i = actnoder[1];
                       int mside_i = inter[k]->OtherSide(sside_i);
                       flags[sside_i*ninter+k] = 3;
                       flags[mside_i*ninter+k] = 2;
                     }
                     // found node on slave/LM side
                     // so set flags that node has to be on master on inter[i]
                     else if (nodeside == inter[k]->OtherSide(mside))
                     {
                       int mside_i = actnoder[1];
                       int sside_i = inter[k]->OtherSide(mside_i);
                       flags[sside_i*ninter+k] = 3;
                       flags[mside_i*ninter+k] = 2;
                     }
                     else
                     {
                       cout << "***ERR*** MRTR::Manager::ChooseMortarSide_2D:\n"
                            << "***ERR*** Matrix M or D is NULL\n"
                            << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
                       exit(EXIT_FAILURE);
                     }
                   }
                   break;
                 }
               } // for (int l=0; l<cornernodes[k].size(); ++l)
             }
         } // for (int k=0; k<ninter; ++k)
       } // for (int j=0; j<cornernodes[i].size(); ++j)

       // sum all flags
       Comm().MaxAll(&(flags[0]),&(flagr[0]),(2*ninter));
       
       cout << "Flags side 0:\n";
       for (int j=0; j<ninter; ++j) cout << "inter " << j << " flag " << flagr[j] << endl;
       cout << "Flags side 1:\n";
       for (int j=ninter; j<2*ninter; ++j) cout << "inter " << j << " flag " << flagr[j] << endl;
       
       // loop through flags and make a decision what to chose on inter[i]
       for (int k=0; k<ninter; ++k)
       {
         if (i==k) continue; 
         // I don't care what is chosen
         if (flagr[k]==0 && flagr[ninter+k]==0); 
         // i don't care either
         else if (flagr[k]==1 && flagr[ninter+k]==1);
         // side 0 has to be master, side 1 has to be slave
         else if (flagr[k]==2 && flagr[ninter+k]==3)
         {
           // check whether a side has been set before and this side is differen?
           if (inter[i]->MortarSide()==1)
           {
             cout << "***ERR*** MRTR::Manager::ChooseMortarSide_2D:\n"
                  << "***ERR*** want to set mortar side to 0 but has been set to 1 before\n"
                  << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
             exit(EXIT_FAILURE);
           }
           inter[i]->SetMortarSide(0);
         }
         // side 0 has to be slave, side 1 has to be master
         else if (flagr[k]==3 && flagr[ninter+k]==2)
         {
           // check whether a side has been set before and this side is differen?
           if (inter[i]->MortarSide()==0)
           {
             cout << "***ERR*** MRTR::Manager::ChooseMortarSide_2D:\n"
                  << "***ERR*** want to set mortar side to 1 but has been set to 0 before\n"
                  << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
             exit(EXIT_FAILURE);
           }
           inter[i]->SetMortarSide(1);
         }
         else
         {
           cout << "***ERR*** MRTR::Manager::ChooseMortarSide_2D:\n"
                << "***ERR*** unknown type of coloring flag: " << flags[k] << "\n"
                << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
           exit(EXIT_FAILURE);
         }
       }
       // Obviously we don't need to care what side inter[i] is going to have
       if (inter[i]->MortarSide()==-2)
         inter[i]->SetMortarSide(0);
         



     } // if (inter[i]->MortarSide() == -2)
     
     

    
    
  } // for (int i=0; i<ninter; ++i)


  // tidy up
  nnodes.clear();
  for (int i=0; i<ninter; ++i) 
    if (nodes[i]) delete [] nodes[i];
  nodes.clear();

  return true;
}




















#endif // TRILINOS_PACKAGE
