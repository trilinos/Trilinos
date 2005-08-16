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
  // this is probably the place to put detection of end nodes and end segments
  

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
  // build the Epetra_CrsMatrix D and M
  {
    if (!constraintsmap_)
    {
      cout << "***ERR*** MRTR::Manager::Mortar_Integrate:\n"
           << "***ERR*** constraintsmap_ == NULL, cannot generate D and M\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    if (D_) delete D_;
    if (M_) delete M_;
    D_ = new Epetra_CrsMatrix(Copy,*constraintsmap_,5,false);
    M_ = new Epetra_CrsMatrix(Copy,*constraintsmap_,40,false);
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
  D_->FillComplete(*inputmap_,*constraintsmap_);
  M_->FillComplete(*inputmap_,*constraintsmap_);

  return true;
}

/*----------------------------------------------------------------------*
 |  Choose dofs for lagrange multipliers (private)           mwgee 07/05|
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

  // create a rowmap for the saddle point problem
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
  
  // create a matrix for the saddle problem and fill it
  if (saddlematrix_) delete saddlematrix_; 
  saddlematrix_ = new Epetra_CrsMatrix(Copy,*saddlemap_,90);


  
  // add values from inputmatrix
  // get a column map to convert column indices to global
  const Epetra_Map& inputcolmap = inputmatrix_->ColMap();
  for (int i=0; i<inputmatrix_->NumMyRows(); ++i)
  {
    // get a row in local indices (FillComplete() has been called on inputmatrix_))
    int numentries;
    double* values;
    int* indices;
    int err = inputmatrix_->ExtractMyRowView(i,numentries,values,indices);
    if (err)
    {
      cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
           << "***ERR*** Proc " << comm_.MyPID() << ": inputmatrix_->ExtractMyRowView returned nonzero\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    
    // get the global row index
    int grow = inputmatrix_->GRID(i);
    if (grow<-1)
    {
      cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
           << "***ERR*** Proc " << comm_.MyPID() << " does not have global row for local row " << i << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    
    
    // loop row and add values to saddlematrix_
    for (int j=0; j<numentries; ++j)
    {
      int gcol = inputcolmap.GID(indices[j]);
      if (gcol<0)
      {
        cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
             << "***ERR*** Proc " << comm_.MyPID() << " does not have global col for local col " << indices[j] << " in column map\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
      //cout << "Proc " << comm_.MyPID() << ": grow " << grow << " gcol " << gcol << " value " << values[j] << endl;
      err = saddlematrix_->InsertGlobalValues(grow,1,&values[j],&gcol);
      if (err)
      {
        cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
             << "***ERR*** Proc " << comm_.MyPID() << " saddlematrix_->InsertGlobalValues returned " << err << "\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
    } // for (int j=0; j<numentries; ++j)
  } // for (int i=0; i<inputmatrix_->NumMyRows(); ++i)




  
  // add values from matrix D_;
  const Epetra_Map& dcolmap = D_->ColMap();
  for (int i=0; i<D_->NumMyRows(); ++i)
  {
    // get a row in local indices (FillComplete() has been called on D_))
    int numentries;
    double* values;
    int* indices;
    int err = D_->ExtractMyRowView(i,numentries,values,indices);
    if (err)
    {
      cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
           << "***ERR*** Proc " << comm_.MyPID() << ": D_->ExtractMyRowView returned nonzero\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    
    // get the global row index
    int grow = D_->GRID(i);
    if (grow<-1)
    {
      cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
           << "***ERR*** Proc " << comm_.MyPID() << " does not have global row for local row " << i << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    
    // loop row and values to saddlematrix_
    for (int j=0; j<numentries; ++j)
    {
      int gcol = dcolmap.GID(indices[j]);
      if (gcol<0)
      {
        cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
             << "***ERR*** Proc " << comm_.MyPID() << " does not have global col for local col " << indices[j] << " in column map\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
      //cout << "Proc " << comm_.MyPID() << ": grow " << grow << " gcol " << gcol << " value " << values[j] << endl;
      err = saddlematrix_->InsertGlobalValues(grow,1,&values[j],&gcol);
      if (err)
      {
        cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
             << "***ERR*** Proc " << comm_.MyPID() << " saddlematrix_->InsertGlobalValues returned " << err << "\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
    } // for (int j=0; j<numentries; ++j)
  } // for (int i=0; i<D_->NumMyRows(); ++i)


  
  
  
  // add values from matrix M_;
  const Epetra_Map& mcolmap = M_->ColMap();
  for (int i=0; i<M_->NumMyRows(); ++i)
  {
    // get a row in local indices (FillComplete() has been called on M_))
    int numentries;
    double* values;
    int* indices;
    int err = M_->ExtractMyRowView(i,numentries,values,indices);
    if (err)
    {
      cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
           << "***ERR*** Proc " << comm_.MyPID() << ": M_->ExtractMyRowView returned nonzero\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    
    // get the global row index
    int grow = M_->GRID(i);
    if (grow<-1)
    {
      cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
           << "***ERR*** Proc " << comm_.MyPID() << " does not have global row for local row " << i << "\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      exit(EXIT_FAILURE);
    }
    
    // loop row and values to saddlematrix_
    for (int j=0; j<numentries; ++j)
    {
      int gcol = mcolmap.GID(indices[j]);
      if (gcol<0)
      {
        cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
             << "***ERR*** Proc " << comm_.MyPID() << " does not have global col for local col " << indices[j] << " in column map\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
      cout << "Proc " << comm_.MyPID() << ": grow " << grow << " gcol " << gcol << " value " << values[j] << endl;
      err = saddlematrix_->InsertGlobalValues(grow,1,&values[j],&gcol);
      if (err)
      {
        cout << "***ERR*** MRTR::Manager::MakeSaddleProblem:\n"
             << "***ERR*** Proc " << comm_.MyPID() << " saddlematrix_->InsertGlobalValues returned " << err << "\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        exit(EXIT_FAILURE);
      }
    } // for (int j=0; j<numentries; ++j)
  } // for (int i=0; i<M_->NumMyRows(); ++i)
  
  
  
  
  
  
  return saddlematrix_;
}

#endif // TRILINOS_PACKAGE
