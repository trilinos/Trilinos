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
# Questions? Contact Glen Hansen (gahanse@sandia.gov)
#
# ************************************************************************
#@HEADER
*/
/* ******************************************************************** */
/* See the file COPYRIGHT for a complete copyright notice, contact      */
/* person and disclaimer.                                               */
/* ******************************************************************** */
#include "mrtr_manager.H"
#include "mrtr_utils.H"
#include "EpetraExt_MatrixMatrix.h"  // for adding matrices
#include "EpetraExt_Transpose_RowMatrix.h"
#include "Epetra_Time.h"

/*----------------------------------------------------------------------*
 |  Store the rowmap of the underlying matrix from the application 07/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::SetProblemMap(const Epetra_Map* map)
{
  problemmap_ = Teuchos::rcp(new Epetra_Map(*map));
  return true;
}

/*----------------------------------------------------------------------*
 |  Set ptr to the uncoupled matrix on input                       07/05|
 |  Note that MOERTEL::Manager does not make a deep copy here by default   |
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::SetInputMatrix(Epetra_CrsMatrix* inputmatrix, bool DeepCopy)
{
  if (DeepCopy)
    inputmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(*inputmatrix));
  else
  {
    inputmatrix_ = Teuchos::rcp(inputmatrix);
    inputmatrix_.release();
  }
  ResetSolver();
  return true;
}


/*----------------------------------------------------------------------*
 | Create the saddle point problem map (private)             mwgee 11/05|
 | Note that this is collective for ALL procs                           |
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::BuildSaddleMap()
{
  // check whether all interfaces are complete and integrated
  std::map<int,Teuchos::RCP<MOERTEL::Interface> >::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    if (curr->second->IsComplete() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::BuildSaddleMap:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not Complete()\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }

  // check whether we have a problemmap_
  if (problemmap_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::BuildSaddleMap:\n"
           << "***ERR*** No problemmap_ set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
  }

  // check whether we have a constraintsmap_
  if (constraintsmap_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::BuildSaddleMap:\n"
           << "***ERR*** onstraintsmap is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
  }

  // the saddle point problem rowmap is the problemmap_ + the constraintmap
  int numglobalelements = problemmap_->NumGlobalElements() +
                          constraintsmap_->NumGlobalElements();
  int nummyelements     = problemmap_->NumMyElements() +
                          constraintsmap_->NumMyElements();
  std::vector<int> myglobalelements(nummyelements);
  int count = 0;
  int* inputmyglobalelements = problemmap_->MyGlobalElements();
  for (int i=0; i<problemmap_->NumMyElements(); ++i)
    myglobalelements[count++] = inputmyglobalelements[i];
  int* constraintsmyglobalelements = constraintsmap_->MyGlobalElements();
  for (int i=0; i<constraintsmap_->NumMyElements(); ++i)
    myglobalelements[count++] = constraintsmyglobalelements[i];
  if (count != nummyelements)
  {
    cout << "***ERR*** MOERTEL::Manager::BuildSaddleMap:\n"
         << "***ERR*** Mismatch in dimensions\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  saddlemap_ = Teuchos::rcp(new Epetra_Map(numglobalelements,nummyelements,
                              &(myglobalelements[0]),0,comm_));
  myglobalelements.clear();

  return true;
}

/*----------------------------------------------------------------------*
 | Create the saddle point problem (public)                  mwgee 07/05|
 | Note that this is collective for ALL procs                           |
 *----------------------------------------------------------------------*/
Epetra_CrsMatrix* MOERTEL::Manager::MakeSaddleProblem()
{
  // time this process
  Epetra_Time time(Comm());
  time.ResetStartTime();

  // check whether all interfaces are complete and integrated
  std::map<int,Teuchos::RCP<MOERTEL::Interface> >::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    if (curr->second->IsComplete() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not Complete()\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
    }
    if (curr->second->IsIntegrated() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not integrated yet\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
    }
  }

  // check whether we have a problemmap_
  if (problemmap_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** No problemmap_ set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check whether we have a constraintsmap_
  if (constraintsmap_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** onstraintsmap is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check for saddlemap_
  if (saddlemap_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** saddlemap_==NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check for inputmatrix
  if (inputmatrix_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** No inputmatrix set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check whether we have M and D matrices
  if (D_==Teuchos::null || M_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** Matrix M or D is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // create a matrix for the saddle problem and fill it
  saddlematrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*saddlemap_,90));

  // add values from inputmatrix
  MOERTEL::MatrixMatrixAdd(*inputmatrix_,false,1.0,*saddlematrix_,0.0);

  // add values from D_
  MOERTEL::MatrixMatrixAdd(*D_,false,1.0,*saddlematrix_,1.0);
  MOERTEL::MatrixMatrixAdd(*D_,true,1.0,*saddlematrix_,1.0);

  // add values from M_
  MOERTEL::MatrixMatrixAdd(*M_,false,1.0,*saddlematrix_,1.0);
  MOERTEL::MatrixMatrixAdd(*M_,true,1.0,*saddlematrix_,1.0);

  saddlematrix_->FillComplete();
  saddlematrix_->OptimizeStorage();

  // time this process
  double t = time.ElapsedTime();
  if (OutLevel()>5 && Comm().MyPID()==0)
    cout << "MOERTEL (Proc 0): Construct saddle system in " << t << " sec\n";


  return saddlematrix_.get();
}

#if 0 // old working version (slow)
/*----------------------------------------------------------------------*
 | Create the spd system (public)                            mwgee 12/05|
 | Note that this is collective for ALL procs                           |
 *----------------------------------------------------------------------*/
Epetra_CrsMatrix* MOERTEL::Manager::MakeSPDProblem()
{
  // time this process
  Epetra_Time time(Comm());
  time.ResetStartTime();

  // check whether all interfaces are complete and integrated
  std::map<int,Teuchos::RCP<MOERTEL::Interface> >::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    if (curr->second->IsComplete() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not Complete()\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
    }
    if (curr->second->IsIntegrated() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not integrated yet\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
    }
  }

  // check whether we have a problemmap_
  if (problemmap_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** No problemmap_ set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check whether we have a constraintsmap_
  if (constraintsmap_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** onstraintsmap is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check for saddlemap_
  if (saddlemap_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** saddlemap_==NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check for inputmatrix
  if (inputmatrix_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** No inputmatrix set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check whether we have M and D matrices
  if (D_==Teuchos::null || M_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** Matrix M or D is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // we need a map from lagrange multiplier dofs to primal dofs on the same node
  std::vector<MOERTEL::Node*> nodes(0);
  std::map<int,int> lm_to_dof;
  for (curr=interface_.begin(); curr!=interface_.end(); ++curr)
  {
	Teuchos::RCP<MOERTEL::Interface> inter = curr->second;
    inter->GetNodeView(nodes);
    for (int i=0; i<(int)nodes.size(); ++i)
    {
      if (!nodes[i]->Nlmdof())
        continue;
      const int* dof = nodes[i]->Dof();
      const int* lmdof = nodes[i]->LMDof();
      for (int j=0; j<nodes[i]->Nlmdof(); ++j)
      {
        //cout << "j " << j << " maps lmdof " << lmdof[j] << " to dof " << dof[j] << endl;
        lm_to_dof[lmdof[j]] = dof[j];
      }
    }
  }
  lm_to_dof_ = Teuchos::rcp(new std::map<int,int>(lm_to_dof)); // this is a very useful map for the moertel_ml_preconditioner
  /*
               _              _
              |               |
              |  Arr  Arn  Mr |
         S =  |               |
              |  Anr  Ann  D  |
              |
              |  MrT  D    0  |
              |_          _   |

               _           _
              |            |
              |  Arr  Arn  |
         A =  |            |
              |  Anr  Ann  |
              |_          _|

        1) Ann is square and we need it's Range/DomainMap annmap

               _         _
        WT =  |_ 0 Dinv _|

        2) Build WT (has rowmap/rangemap annmap and domainmap problemmap_)

               _    _
              |     |
              |  Mr |
         B =  |     |
              |  D  |
              |_   _|

        3) Build B (has rowmap/rangemap problemmap_ and domainmap annmap)

        4) Build I, the identity matrix with maps problemmap_,problemmap_);

        After constructing WT ,B and I we can start building Atilde (spdmatrix_)

           Atilde = A + ( B WT - I) A W B^T + B WT A (W B^T - I)

        5) Build BWT = B * WT

        6) Build BWTmI = BWT - I

        7) Build BWTmIAWBT = BWTmI * A * W * B^T

        8) Allocate spdmatrix_  = A + BWTmIAWBT

        9) Build WBTmI = WT^T * B^T - I

        10) Build BWTAWBTmI = BWT * A * WBTmI and add to spdmatrix_
            Call FillComplete on spdmatrix_

        11) Build ImBWT = I - BWT and store it

  */

  int err=0;
  //--------------------------------------------------------------------------
  // 1) create the rangemap of Ann
  std::vector<int> myanngids(problemmap_->NumMyElements());
  int count=0;
  std::map<int,int>::iterator intintcurr;
  for (intintcurr=lm_to_dof.begin(); intintcurr!=lm_to_dof.end(); ++intintcurr)
  {
    if (problemmap_->MyGID(intintcurr->second)==false)
      continue;
    if ((int)myanngids.size()<=count)
      myanngids.resize(myanngids.size()+50);
    myanngids[count] = intintcurr->second;
    ++count;
  }
  myanngids.resize(count);
    int numglobalelements;
  Comm().SumAll(&count,&numglobalelements,1);
  Epetra_Map* annmap = new Epetra_Map(numglobalelements,count,&myanngids[0],0,Comm());
  annmap_ = Teuchos::rcp(annmap);
  myanngids.clear();

#if 0
  //--------------------------------------------------------------------------
  // 1.5) split matrix into blocks Arr Arn Anr Ann
  Teuchos::RCP<Epetra_Map>       A11row = Teuchos::null;
  Teuchos::RCP<Epetra_Map>       A22row = annmap_;
  Teuchos::RCP<Epetra_CrsMatrix> A11    = Teuchos::null;
  Teuchos::RCP<Epetra_CrsMatrix> A12    = Teuchos::null;
  Teuchos::RCP<Epetra_CrsMatrix> A21    = Teuchos::null;
  Teuchos::RCP<Epetra_CrsMatrix> A22    = Teuchos::null;
  MOERTEL::SplitMatrix2x2(inputmatrix_,A11row,A22row,A11,A12,A21,A22);
#endif

#if 0
  //--------------------------------------------------------------------------
  // 1.7) create a shifted version of M and D
  Epetra_CrsMatrix* MTshifted = new Epetra_CrsMatrix(Copy,*annmap,1,false);
  Epetra_CrsMatrix* Dshifted  = new Epetra_CrsMatrix(Copy,*annmap,1,false);
  std::vector<int> gindices(500);
  for (intintcurr=lm_to_dof.begin(); intintcurr!=lm_to_dof.end(); ++intintcurr)
  {
    const int lmdof = intintcurr->first;
    const int dof   = intintcurr->second;
    if (D_->MyGRID(lmdof)==false)
      continue;
    const int lmlrid = D_->LRID(lmdof);
    int numentries;
    int* indices;
    double* values;

    // do D
    err = D_->ExtractMyRowView(lmlrid,numentries,values,indices);
    if (err) cout << "D_->ExtractMyRowView returned err=" << err << endl;
    if (numentries>(int)gindices.size()) gindices.resize(numentries);
    for (int j=0; j<numentries; ++j)
    {
      gindices[j] = D_->GCID(indices[j]);
      if (gindices[j]<0) cout << "Cannot find gcid for indices[j]\n";
    }
    err = Dshifted->InsertGlobalValues(dof,numentries,values,&gindices[0]);
    if (err<0) cout << "Dshifted->InsertGlobalValues returned err=" << err << endl;

    // do MT
    err = M_->ExtractMyRowView(lmlrid,numentries,values,indices);
    if (err) cout << "M_->ExtractMyRowView returned err=" << err << endl;
    if (numentries>(int)gindices.size()) gindices.resize(numentries);
    for (int j=0; j<numentries; ++j)
    {
      gindices[j] = M_->GCID(indices[j]);
      if (gindices[j]<0) cout << "Cannot find gcid for indices[j]\n";
    }
    err = MTshifted->InsertGlobalValues(dof,numentries,values,&gindices[0]);
    if (err<0) cout << "MTshifted->InsertGlobalValues returned err=" << err << endl;

  }
  gindices.clear();
  Dshifted->FillComplete(*problemmap_,*annmap);
  Dshifted->OptimizeStorage();
  MTshifted->FillComplete(*problemmap_,*annmap);
  MTshifted->OptimizeStorage();
  Dshifted_ = Teuchos::rcp(Dshifted);
  MTshifted_ = Teuchos::rcp(MTshifted);
#endif

  //--------------------------------------------------------------------------
  // 2) create WT, D and MT
  Epetra_CrsMatrix* WT        = new Epetra_CrsMatrix(Copy,*annmap,1,false);
  for (intintcurr=lm_to_dof.begin(); intintcurr!=lm_to_dof.end(); ++intintcurr)
  {
    int lmdof = intintcurr->first;
    int dof   = intintcurr->second;
    if (D_->MyGRID(lmdof)==false)
      continue;
    int lmlrid = D_->LRID(lmdof);
    int numentries;
    int* indices;
    double* values;
    err = D_->ExtractMyRowView(lmlrid,numentries,values,indices);
    if (err) cout << "D_->ExtractMyRowView returned err=" << err << endl;
    bool foundit = false;
    for (int j=0; j<numentries; ++j)
    {
      int gcid = D_->GCID(indices[j]);
      if (gcid<0) cout << "Cannot find gcid for indices[j]\n";
      //cout << "Proc " << Comm().MyPID() << " lmdof " << lmdof << " dof " << dof << " gcid " << gcid << " val " << values[j] << endl;
      if (gcid==dof)
      {
        double val = 1./values[j];
        err = WT->InsertGlobalValues(dof,1,&val,&dof);
        if (err<0) cout << "WT->InsertGlobalValues returned err=" << err << endl;
        foundit = true;
        break;
      }
    }
    if (!foundit)
    {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** Cannot compute inverse of D_\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      cout << "lmdof " << lmdof << " dof " << dof << endl;
      return NULL;
    }
  }
  WT->FillComplete(*problemmap_,*annmap);
  WT_ = Teuchos::rcp(WT);

  //--------------------------------------------------------------------------
  // 3) create B
  // create a temporary matrix with rowmap of the Ann block
  Epetra_CrsMatrix* tmp = new Epetra_CrsMatrix(Copy,*annmap,120);
  std::vector<int> newindices(100);
  for (intintcurr=lm_to_dof.begin(); intintcurr!=lm_to_dof.end(); ++intintcurr)
  {
    int lmdof = intintcurr->first;
    int dof   = intintcurr->second;
    if (D_->MyGRID(lmdof)==false)
      continue;
    int lmlrid = D_->LRID(lmdof);
    if (lmlrid<0) cout << "Cannot find lmlrid for lmdof\n";
    int numentries;
    int* indices;
    double* values;

    // extract and add values from D
    err = D_->ExtractMyRowView(lmlrid,numentries,values,indices);
    if (err) cout << "D_->ExtractMyRowView returned err=" << err << endl;
    if (numentries>(int)newindices.size()) newindices.resize(numentries);
    for (int j=0; j<numentries; ++j)
    {
      newindices[j] = D_->GCID(indices[j]);
      if (newindices[j]<0) cout << "Cannot find gcid for indices[j]\n";
    }
    //cout << "Inserting from D in row " << dof << " cols/val ";
    //for (int j=0; j<numentries; ++j) cout << newindices[j] << "/" << values[j] << " ";
    //cout << endl;
    err = tmp->InsertGlobalValues(dof,numentries,values,&newindices[0]);
    if (err) cout << "tmp->InsertGlobalValues returned err=" << err << endl;

    // extract and add values from M
    err = M_->ExtractMyRowView(lmlrid,numentries,values,indices);
    if (err) cout << "M_->ExtractMyRowView returned err=" << err << endl;
    if (numentries>(int)newindices.size()) newindices.resize(numentries);
    for (int j=0; j<numentries; ++j)
    {
      newindices[j] = M_->GCID(indices[j]);
      if (newindices[j]<0) cout << "Cannot find gcid for indices[j]\n";
    }
    //cout << "Inserting from M in row " << dof << " cols/val ";
    //for (int j=0; j<numentries; ++j) cout << newindices[j] << "/" << values[j] << " ";
    //cout << endl;
    err = tmp->InsertGlobalValues(dof,numentries,values,&newindices[0]);
    if (err) cout << "tmp->InsertGlobalValues returned err=" << err << endl;
  }
  tmp->FillComplete(*(problemmap_.get()),*annmap);

  // B is transposed of tmp
  EpetraExt::RowMatrix_Transpose* trans = new EpetraExt::RowMatrix_Transpose(false);
  Epetra_CrsMatrix* B = &(dynamic_cast<Epetra_CrsMatrix&>(((*trans)(const_cast<Epetra_CrsMatrix&>(*tmp)))));
  delete tmp; tmp = NULL;
  B_ = Teuchos::rcp(new Epetra_CrsMatrix(*B));
  newindices.clear();

  //--------------------------------------------------------------------------
  // 4) create I
  Epetra_CrsMatrix* I = new Epetra_CrsMatrix(Copy,*problemmap_,1,true);
  for (int i=0; i<I->NumMyRows(); ++i)
  {
    double one = 1.0;
    int grid = I->GRID(i);
    if (grid<0) cout << "Cannot find grid for i\n";
    err = I->InsertGlobalValues(grid,1,&one,&grid);
    if (err<0) cout << "I->InsertGlobalValues returned err=" << err << endl;
  }
  I->FillComplete(*problemmap_,*problemmap_);
  I_ = Teuchos::rcp(I);
  //--------------------------------------------------------------------------
  // 5) Build BWT = B * WT
  Epetra_CrsMatrix* BWT = MOERTEL::MatMatMult(*B,false,*WT,false,OutLevel());

  //--------------------------------------------------------------------------
  // 6) Build BWTmI = BWT - I
  Epetra_CrsMatrix* BWTmI = new Epetra_CrsMatrix(Copy,*problemmap_,10,false);
  MOERTEL::MatrixMatrixAdd(*BWT,false,1.0,*BWTmI,0.0);
  MOERTEL::MatrixMatrixAdd(*I,false,-1.0,*BWTmI,1.0);
  BWTmI->FillComplete();

  //--------------------------------------------------------------------------
  // 7) Build BWTmIAWBT = BWTmI * A * W * B^T
  Epetra_CrsMatrix* BWTmIA = MOERTEL::MatMatMult(*BWTmI,false,*inputmatrix_,false,OutLevel());
  Epetra_CrsMatrix* WBT    = MOERTEL::MatMatMult(*WT,true,*B,true,OutLevel());
  Epetra_CrsMatrix* BWTmIAWBT = MOERTEL::MatMatMult(*BWTmIA,false,*WBT,false,OutLevel());
  delete BWTmIA; BWTmIA = NULL;

  //--------------------------------------------------------------------------
  // 8) Allocate spdmatrix_ and add A and BWTmIAWBT
  spdmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*problemmap_,10,false));
  MOERTEL::MatrixMatrixAdd(*BWTmIAWBT,false,1.0,*spdmatrix_,0.0);
  delete BWTmIAWBT; BWTmIAWBT = NULL;
  MOERTEL::MatrixMatrixAdd(*inputmatrix_,false,1.0,*spdmatrix_,1.0);

  //--------------------------------------------------------------------------
  // 9) Build WBTmI = WT^T * B^T - I
  Epetra_CrsMatrix* WBTmI = new Epetra_CrsMatrix(Copy,*problemmap_,10,false);
  MOERTEL::MatrixMatrixAdd(*I,false,-1.0,*WBTmI,0.0);
  MOERTEL::MatrixMatrixAdd(*WBT,false,1.0,*WBTmI,1.0);
  WBTmI->FillComplete();
  delete WBT; WBT = NULL;

  //--------------------------------------------------------------------------
  // 10) Build BWTAWBTmI = BWT * A * WBTmI and add to spdmatrix_
  Epetra_CrsMatrix* BWTA      = MOERTEL::MatMatMult(*BWT,false,*inputmatrix_,false,OutLevel());
  Epetra_CrsMatrix* BWTAWBTmI = MOERTEL::MatMatMult(*BWTA,false,*WBTmI,false,OutLevel());
  delete BWTA; BWTA = NULL;
  delete WBTmI; WBTmI = NULL;
  MOERTEL::MatrixMatrixAdd(*BWTAWBTmI,false,1.0,*spdmatrix_,1.0);
  delete BWTAWBTmI; BWTAWBTmI = NULL;
  spdmatrix_->FillComplete();
  spdmatrix_->OptimizeStorage();

  //--------------------------------------------------------------------------
  // 11) Build ImBWT = I - BWT and store it as spdrhs_
  spdrhs_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*problemmap_,10,false));
  MOERTEL::MatrixMatrixAdd(*I,false,1.0,*spdrhs_,0.0);
  //delete I; I = NULL;
  MOERTEL::MatrixMatrixAdd(*BWT,false,-1.0,*spdrhs_,1.0);
  delete BWT; BWT = NULL;
  spdrhs_->FillComplete();
  spdrhs_->OptimizeStorage();

  //--------------------------------------------------------------------------
  // tidy up
  lm_to_dof.clear();
  //delete annmap;
  //delete WT; WT = NULL;
  delete trans; B = NULL;

  // time this process
  double t = time.ElapsedTime();
  if (OutLevel()>5 && Comm().MyPID()==0)
    cout << "MOERTEL (Proc 0): Construct spd system in " << t << " sec\n";

  return spdmatrix_.get();
}
#endif

#if 1 // new version (faster)
/*----------------------------------------------------------------------*
 | Create the spd system (public)                            mwgee 12/05|
 | Note that this is collective for ALL procs                           |
 *----------------------------------------------------------------------*/
Epetra_CrsMatrix* MOERTEL::Manager::MakeSPDProblem()
{
  // time this process
  Epetra_Time time(Comm());
  time.ResetStartTime();

  // check whether all interfaces are complete and integrated
  std::map<int,Teuchos::RCP<MOERTEL::Interface> >::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    if (curr->second->IsComplete() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not Complete()\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
    }
    if (curr->second->IsIntegrated() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not integrated yet\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
    }
  }

  // check whether we have a problemmap_
  if (problemmap_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** No problemmap_ set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check whether we have a constraintsmap_
  if (constraintsmap_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** onstraintsmap is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check for saddlemap_
  if (saddlemap_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** saddlemap_==NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check for inputmatrix
  if (inputmatrix_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** No inputmatrix set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check whether we have M and D matrices
  if (D_==Teuchos::null || M_==Teuchos::null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** Matrix M or D is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // we need a map from lagrange multiplier dofs to primal dofs on the same node
  std::vector<MOERTEL::Node*> nodes(0);
  std::map<int,int> lm_to_dof;
  for (curr=interface_.begin(); curr!=interface_.end(); ++curr)
  {
	Teuchos::RCP<MOERTEL::Interface> inter = curr->second;
    inter->GetNodeView(nodes);
    for (int i=0; i<(int)nodes.size(); ++i)
    {
      if (!nodes[i]->Nlmdof())
        continue;
      const int* dof = nodes[i]->Dof();
      const int* lmdof = nodes[i]->LMDof();
      for (int j=0; j<nodes[i]->Nlmdof(); ++j)
      {
        //cout << "j " << j << " maps lmdof " << lmdof[j] << " to dof " << dof[j] << endl;
        lm_to_dof[lmdof[j]] = dof[j];
      }
    }
  }
  lm_to_dof_ = Teuchos::rcp(new std::map<int,int>(lm_to_dof)); // this is a very useful map for the moertel_ml_preconditioner
  /*
               _              _
              |               |
              |  Arr  Arn  Mr |
         S =  |               |
              |  Anr  Ann  D  |
              |
              |  MrT  D    0  |
              |_          _   |

               _           _
              |            |
              |  Arr  Arn  |
         A =  |            |
              |  Anr  Ann  |
              |_          _|

        1) Ann is square and we need it's Range/DomainMap annmap

               _         _
        WT =  |_ 0 Dinv _|

        2) Build WT (has rowmap/rangemap annmap and domainmap problemmap_)

               _    _
              |     |
              |  Mr |
         B =  |     |
              |  D  |
              |_   _|

        3) Build B (has rowmap/rangemap problemmap_ and domainmap annmap)

        4) Build I, the identity matrix with maps problemmap_,problemmap_);

  */

  int err=0;
  //--------------------------------------------------------------------------
  // 1) create the rangemap of Ann
  std::vector<int> myanngids(problemmap_->NumMyElements());
  int count=0;
  std::map<int,int>::iterator intintcurr;
  for (intintcurr=lm_to_dof.begin(); intintcurr!=lm_to_dof.end(); ++intintcurr)
  {
    if (problemmap_->MyGID(intintcurr->second)==false)
      continue;
    if ((int)myanngids.size()<=count)
      myanngids.resize(myanngids.size()+50);
    myanngids[count] = intintcurr->second;
    ++count;
  }
  myanngids.resize(count);
    int numglobalelements;
  Comm().SumAll(&count,&numglobalelements,1);
  Epetra_Map* annmap = new Epetra_Map(numglobalelements,count,&myanngids[0],0,Comm());
  annmap_ = Teuchos::rcp(annmap);
  myanngids.clear();


  //--------------------------------------------------------------------------
  // 2) create WT
  Epetra_CrsMatrix* WT        = new Epetra_CrsMatrix(Copy,*annmap,1,false);
  for (intintcurr=lm_to_dof.begin(); intintcurr!=lm_to_dof.end(); ++intintcurr)
  {
    int lmdof = intintcurr->first;
    int dof   = intintcurr->second;
    if (D_->MyGRID(lmdof)==false)
      continue;
    int lmlrid = D_->LRID(lmdof);
    int numentries;
    int* indices;
    double* values;
    err = D_->ExtractMyRowView(lmlrid,numentries,values,indices);
    if (err) cout << "D_->ExtractMyRowView returned err=" << err << endl;
    bool foundit = false;
    for (int j=0; j<numentries; ++j)
    {
      int gcid = D_->GCID(indices[j]);
      if (gcid<0) cout << "Cannot find gcid for indices[j]\n";
      //cout << "Proc " << Comm().MyPID() << " lmdof " << lmdof << " dof " << dof << " gcid " << gcid << " val " << values[j] << endl;
      if (gcid==dof)
      {
        double val = 1./values[j];
        err = WT->InsertGlobalValues(dof,1,&val,&dof);
        if (err<0) cout << "WT->InsertGlobalValues returned err=" << err << endl;
        foundit = true;
        break;
      }
    }
    if (!foundit)
    {
      cout << "***ERR*** MOERTEL::Manager::MakeSPDProblem:\n"
           << "***ERR*** Cannot compute inverse of D_\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      cout << "lmdof " << lmdof << " dof " << dof << endl;
      return NULL;
    }
  }
  WT->FillComplete(*problemmap_,*annmap);
  WT_ = Teuchos::rcp(WT);

  //--------------------------------------------------------------------------
  // 3) create B
  // create a temporary matrix with rowmap of the Ann block
  Epetra_CrsMatrix* tmp = new Epetra_CrsMatrix(Copy,*annmap,120);
  std::vector<int> newindices(100);
  for (intintcurr=lm_to_dof.begin(); intintcurr!=lm_to_dof.end(); ++intintcurr)
  {
    int lmdof = intintcurr->first;
    int dof   = intintcurr->second;
    if (D_->MyGRID(lmdof)==false)
      continue;
    int lmlrid = D_->LRID(lmdof);
    if (lmlrid<0) cout << "Cannot find lmlrid for lmdof\n";
    int numentries;
    int* indices;
    double* values;

    // extract and add values from D
    err = D_->ExtractMyRowView(lmlrid,numentries,values,indices);
    if (err) cout << "D_->ExtractMyRowView returned err=" << err << endl;
    if (numentries>(int)newindices.size()) newindices.resize(numentries);
    for (int j=0; j<numentries; ++j)
    {
      newindices[j] = D_->GCID(indices[j]);
      if (newindices[j]<0) cout << "Cannot find gcid for indices[j]\n";
    }
    //cout << "Inserting from D in row " << dof << " cols/val ";
    //for (int j=0; j<numentries; ++j) cout << newindices[j] << "/" << values[j] << " ";
    //cout << endl;
    err = tmp->InsertGlobalValues(dof,numentries,values,&newindices[0]);
    if (err) cout << "tmp->InsertGlobalValues returned err=" << err << endl;

    // extract and add values from M
    err = M_->ExtractMyRowView(lmlrid,numentries,values,indices);
    if (err) cout << "M_->ExtractMyRowView returned err=" << err << endl;
    if (numentries>(int)newindices.size()) newindices.resize(numentries);
    for (int j=0; j<numentries; ++j)
    {
      newindices[j] = M_->GCID(indices[j]);
      if (newindices[j]<0) cout << "Cannot find gcid for indices[j]\n";
    }
    //cout << "Inserting from M in row " << dof << " cols/val ";
    //for (int j=0; j<numentries; ++j) cout << newindices[j] << "/" << values[j] << " ";
    //cout << endl;
    err = tmp->InsertGlobalValues(dof,numentries,values,&newindices[0]);
    if (err) cout << "tmp->InsertGlobalValues returned err=" << err << endl;
  }
  tmp->FillComplete(*(problemmap_.get()),*annmap);

  // B is transposed of tmp
  EpetraExt::RowMatrix_Transpose* trans = new EpetraExt::RowMatrix_Transpose(false);
  Epetra_CrsMatrix* B = &(dynamic_cast<Epetra_CrsMatrix&>(((*trans)(const_cast<Epetra_CrsMatrix&>(*tmp)))));
  delete tmp; tmp = NULL;
  B_ = Teuchos::rcp(new Epetra_CrsMatrix(*B));
  newindices.clear();

  //--------------------------------------------------------------------------
  // 4) create I
  Epetra_CrsMatrix* I = MOERTEL::PaddedMatrix(*problemmap_,1.0,1);
  I->FillComplete(*problemmap_,*problemmap_);
  I_ = Teuchos::rcp(I);

  //--------------------------------------------------------------------------
  // 5) Build BWT = B * WT
  Epetra_CrsMatrix* BWT = MOERTEL::MatMatMult(*B,false,*WT,false,OutLevel());

  //--------------------------------------------------------------------------
  // 6) create CBT = (I-BWT)*A*W*BT
  Epetra_CrsMatrix* spdrhs = new Epetra_CrsMatrix(Copy,*problemmap_,10,false);
  MOERTEL::MatrixMatrixAdd(*BWT,false,-1.0,*spdrhs,0.0);
  MOERTEL::MatrixMatrixAdd(*I,false,1.0,*spdrhs,1.0);
  spdrhs->FillComplete();
  spdrhs->OptimizeStorage();
  spdrhs_ = Teuchos::rcp(spdrhs);

  Epetra_CrsMatrix* tmp1 = MOERTEL::MatMatMult(*inputmatrix_,false,*WT,true,OutLevel());
  Epetra_CrsMatrix* tmp2 = MOERTEL::MatMatMult(*tmp1,false,*B,true,OutLevel());
  delete tmp1;
  Epetra_CrsMatrix* CBT  = MOERTEL::MatMatMult(*spdrhs,false,*tmp2,false,OutLevel());
  delete tmp2;

  //--------------------------------------------------------------------------
  // 7) create spdmatrix_ = A - CBT - (CBT)^T
  spdmatrix_ = Teuchos::rcp(new Epetra_CrsMatrix(Copy,*problemmap_,10,false));
  MOERTEL::MatrixMatrixAdd(*inputmatrix_,false,1.0,*spdmatrix_,0.0);
  MOERTEL::MatrixMatrixAdd(*CBT,false,-1.0,*spdmatrix_,1.0);
  MOERTEL::MatrixMatrixAdd(*CBT,true,-1.0,*spdmatrix_,1.0);
  delete CBT; CBT = NULL;
  spdmatrix_->FillComplete();
  spdmatrix_->OptimizeStorage();

  //--------------------------------------------------------------------------
  // tidy up
  lm_to_dof.clear();
  delete trans; B = NULL;
  // time this process
  double t = time.ElapsedTime();
  if (OutLevel()>5 && Comm().MyPID()==0)
    cout << "MOERTEL (Proc 0): Construct spd system in " << t << " sec\n";

  return spdmatrix_.get();
}
#endif

/*----------------------------------------------------------------------*
 | Return the right hand side matrix of the reduced system   mwgee 01/06|
 *----------------------------------------------------------------------*/
Epetra_CrsMatrix* MOERTEL::Manager::GetRHSMatrix()
{
  if (spdrhs_ != Teuchos::null)
    return spdrhs_.get();
  else
  {
    Epetra_CrsMatrix* tmp = MakeSPDProblem();
    if (!tmp)
    {
      cout << "***ERR*** MOERTEL::Manager::GetRHSMatrix:\n"
           << "***ERR*** Cannot compute reduced system of equations\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
    }
    return spdrhs_.get();
  }
}

/*----------------------------------------------------------------------*
 | Reset the internal solver (public)                        mwgee 12/05|
 *----------------------------------------------------------------------*/
void MOERTEL::Manager::ResetSolver()
{
  solverparams_ = Teuchos::null;
  solver_ = Teuchos::null;
  saddlematrix_ = Teuchos::null;
  spdmatrix_ = Teuchos::null;
  spdrhs_ = Teuchos::null;
  I_ = Teuchos::null;
  WT_ = Teuchos::null;
  B_ = Teuchos::null;
  return;
}

/*----------------------------------------------------------------------*
 | Set solver parameters (public)                            mwgee 12/05|
 *----------------------------------------------------------------------*/
void MOERTEL::Manager::SetSolverParameters(Teuchos::ParameterList& params)
{
  solverparams_ = Teuchos::rcp(&params);
  // MOERTEL does not take ownership of the parameterlist
  solverparams_.release();
  return;
}

/*----------------------------------------------------------------------*
 | solve problem (public)                                    mwgee 12/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::Solve(Teuchos::ParameterList& params, Epetra_Vector& sol, const Epetra_Vector& rhs)
{
  SetSolverParameters(params);
  return Solve(sol,rhs);
}

/*----------------------------------------------------------------------*
 | solve problem (public)                                    mwgee 12/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::Solve(Epetra_Vector& sol, const Epetra_Vector& rhs)
{
  // test for solver parameters
  if (solverparams_==Teuchos::null)
  {
    cout << "***ERR*** MOERTEL::Manager::Solve:\n"
         << "***ERR*** No solver parameters set, use SetSolverParameters(ParameterList& params)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  // test for problemmap_
  if (problemmap_==Teuchos::null)
  {
    cout << "***ERR*** MOERTEL::Manager::Solve:\n"
         << "***ERR*** No problem map set, use SetProblemMap(Epetra_Map* map)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  // test for inputmatrix
  if (inputmatrix_==Teuchos::null)
  {
    cout << "***ERR*** MOERTEL::Manager::Solve:\n"
         << "***ERR*** No inputmatrix set, use SetInputMatrix(Epetra_CrsMatrix* inputmatrix, bool DeepCopy = false)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  // test whether problemmap_ matches RangeMap() of inputmatrix
  if (!problemmap_->PointSameAs(inputmatrix_->RangeMap()))
  {
    cout << "***ERR*** MOERTEL::Manager::Solve:\n"
         << "***ERR*** problem map does not match range map of input matrix\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  // test whether maps of rhs and sol are ok
  if (!problemmap_->PointSameAs(rhs.Map()) ||
      !problemmap_->PointSameAs(sol.Map()) )
  {
    cout << "***ERR*** MOERTEL::Manager::Solve:\n"
         << "***ERR*** problem map does not match map of rhs and/or sol\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  // test whether interfaces are complete and have been integrated
  std::map<int,Teuchos::RCP<MOERTEL::Interface> >::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    if (curr->second->IsComplete() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::Solve:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not IsComplete()\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    if (curr->second->IsIntegrated() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::Solve:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not integrated yet\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }

  // test whether we have M and D matrix
  if (D_==Teuchos::null || M_==Teuchos::null)
  {
    cout << "***ERR*** MOERTEL::Manager::Solve:\n"
         << "***ERR*** Matrix M or D is NULL\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  //---------------------------------------------------------------------------
  // make solution and rhs vector matching the system
  Teuchos::RCP<Epetra_Vector> b = Teuchos::rcp(const_cast<Epetra_Vector*>(&rhs));
  b.release();
  Teuchos::RCP<Epetra_Vector> x = Teuchos::rcp(&sol);
  x.release();

  //---------------------------------------------------------------------------
  // get type of system to be used/generated
  Teuchos::RCP<Epetra_CrsMatrix> matrix = Teuchos::null;
  string system = solverparams_->get("System","None");
  if (system=="None")
  {
    cout << "***WRN*** MOERTEL::Manager::Solve:\n"
         << "***WRN*** parameter 'System' was not chosen, using default\n"
         << "***WRN*** which is 'SaddleSystem'\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    solverparams_->set("System","SaddleSystem");
    system = "SaddleSystem";
  }

  //---------------------------------------------------------------------------
  // build a saddle point system
  if (system=="SaddleSystem"  || system=="saddlesystem"  || system=="SADDLESYSTEM" ||
      system=="Saddle_System" || system=="saddle_system" || system=="SADDLE_SYSTEM")
  {
    if (saddlematrix_==Teuchos::null)
    {
      Epetra_CrsMatrix* tmp = MakeSaddleProblem();
      if (!tmp)
      {
        cout << "***ERR*** MOERTEL::Manager::Solve:\n"
             << "***ERR*** MakeSaddleProblem() returned NULL\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
    matrix = saddlematrix_;
    b = Teuchos::rcp(new Epetra_Vector(*saddlemap_,true));
    b.set_has_ownership();
    x = Teuchos::rcp(new Epetra_Vector(*saddlemap_,false));
    x.set_has_ownership();
    for (int i=0; i<rhs.MyLength(); ++i)
    {
      (*b)[i] = rhs[i];
      (*x)[i] = sol[i];
    }
  }

  //---------------------------------------------------------------------------
  // build a spd system
  else if (system=="SPDSystem"  || system=="spdsystem" || system=="spd_system" ||
           system=="SPD_System" || system=="SPDSYSTEM" || system=="SPD_SYSTEM")
  {
    if (spdmatrix_==Teuchos::null)
    {
      Epetra_CrsMatrix* tmp = MakeSPDProblem();
      if (!tmp)
      {
        cout << "***ERR*** MOERTEL::Manager::Solve:\n"
             << "***ERR*** MakeSPDProblem() returned NULL\n"
             << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
        return false;
      }
    }
    matrix = spdmatrix_;
    // we have to multiply the rhs vector b with spdrhs_ to fit with spdmatrix_
    if (spdrhs_==Teuchos::null)
    {
      cout << "***ERR*** MOERTEL::Manager::Solve:\n"
           << "***ERR*** Cannot build righthandside for spd problem\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    Epetra_Vector *tmp = new Epetra_Vector(b->Map(),false);
    int err = spdrhs_->Multiply(false,*b,*tmp);
    if (err)
    {
      cout << "***ERR*** MOERTEL::Manager::Solve:\n"
           << "***ERR*** spdrhs_->Multiply returned err = " << err << endl
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    b = Teuchos::rcp(tmp);
    b.set_has_ownership();
  }

  //---------------------------------------------------------------------------
  // unknown parameter "System"
  else
  {
    cout << "***ERR*** MOERTEL::Manager::Solve:\n"
         << "***ERR*** Unknown type of parameter 'System': " << system << "\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }

  //---------------------------------------------------------------------------
  // create a mortar solver class instance
  if (solver_==Teuchos::null)
    solver_ = Teuchos::rcp(new MOERTEL::Solver(Comm(),OutLevel()));

  //---------------------------------------------------------------------------
  // solve
  bool ok = solver_->Solve(solverparams_,matrix,x,b,*this);
  if (!ok)
  {
    if (Comm().MyPID()==0)
    cout << "***WRN*** MOERTEL::Manager::Solve:\n"
         << "***WRN*** MOERTEL::Solver::Solve returned an error\n"
         << "***WRN*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
  }

  //---------------------------------------------------------------------------
  // copy solution back to sol if neccessary
  if (x.has_ownership())
  {
    for (int i=0; i<sol.MyLength(); ++i)
      sol[i] = (*x)[i];
  }

  return ok;
}


