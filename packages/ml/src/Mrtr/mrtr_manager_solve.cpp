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
#include "mrtr_manager.H"
#include "EpetraExt_MatrixMatrix.h"  // for adding matrices
#include <EpetraExt_Transpose_RowMatrix.h>
#include "Epetra_Time.h"

/*----------------------------------------------------------------------*
 |  Store the rowmap of the underlying matrix from the application 07/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::SetProblemMap(Epetra_Map* map)
{
  problemmap_ = rcp(new Epetra_Map(*map));
  return true;
}

/*----------------------------------------------------------------------*
 |  Set ptr to the uncoupled matrix on input                       07/05|
 |  Note that MOERTEL::Manager does not make a deep copy here by default   |
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::SetInputMatrix(Epetra_CrsMatrix* inputmatrix, bool DeepCopy)
{
  if (DeepCopy)
  {
    inputmatrix_ = rcp(new Epetra_CrsMatrix(*inputmatrix));
    return true;
  }
  else
  {
    inputmatrix_ = rcp(inputmatrix);
    inputmatrix_.release();
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
int MOERTEL::Manager::MatrixMatrixAdd(const Epetra_CrsMatrix& A, bool transposeA,double scalarA,
                      Epetra_CrsMatrix& B,double scalarB )
{
  //
  //This method forms the matrix-matrix sum B = scalarA * op(A) + scalarB * B, where

  if (!A.Filled())
  {
     cout << "***ERR*** MOERTEL::Manager::MatrixMatrixAdd:\n"
          << "***ERR*** FillComplete was not called on A\n"
          << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
     exit(EXIT_FAILURE);
  }
  if (B.Filled())
  {
     cout << "***ERR*** MOERTEL::Manager::MatrixMatrixAdd:\n"
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
  vector<int> Indices(MaxNumEntries);
  vector<double> Values(MaxNumEntries);

  int NumMyRows = A.NumMyRows();
  int Row, err;

  if( scalarA )
  {
    for( int i = 0; i < NumMyRows; ++i )
    {
      Row = Aprime->GRID(i);
      EPETRA_CHK_ERR(Aprime->ExtractGlobalRowCopy(Row,MaxNumEntries,NumEntries,&Values[0],&Indices[0]));
      if( scalarA != 1.0 )
        for( int j = 0; j < NumEntries; ++j ) Values[j] *= scalarA;
      if( B.Filled() ) {//Sum In Values
        err = B.SumIntoGlobalValues( Row, NumEntries, &Values[0], &Indices[0] );
        assert( err == 0 );
      }
      else {
        err = B.InsertGlobalValues( Row, NumEntries, &Values[0], &Indices[0] );
        assert( err == 0 || err == 1 );
      }
    }
  }

  Indices.clear();
  Values.clear();

  if( Atrans ) delete Atrans;

  return(0);
}

/*----------------------------------------------------------------------*
 | Create the saddle point problem map (private)             mwgee 11/05|
 | Note that this is collective for ALL procs                           |
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::BuildSaddleMap()
{
  // check whether all interfaces are complete and integrated
  map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
  for (curr=interface_.begin(); curr != interface_.end(); ++curr)
  {
    if (curr->second->IsComplete() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::BuildSaddleMap:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not Complete()\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
    if (curr->second->IsIntegrated() == false)
    {
      cout << "***ERR*** MOERTEL::Manager::BuildSaddleMap:\n"
           << "***ERR*** interface " << curr->second->Id() << " is not integrated yet\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
    }
  }
  
  // check whether we have a problemmap_
  if (problemmap_==null)
  {
      cout << "***ERR*** MOERTEL::Manager::BuildSaddleMap:\n"
           << "***ERR*** No problemmap_ set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return false;
  }
  
  // check whether we have a constraintsmap_
  if (constraintsmap_==null)
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
  vector<int> myglobalelements(nummyelements);
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
  saddlemap_ = rcp(new Epetra_Map(numglobalelements,nummyelements,
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
  // check whether all interfaces are complete and integrated
  map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
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
  if (problemmap_==null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** No problemmap_ set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }
  
  // check whether we have a constraintsmap_
  if (constraintsmap_==null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** onstraintsmap is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }
  
  // check for saddlemap_
  if (saddlemap_==null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** saddlemap_==NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check for inputmatrix
  if (inputmatrix_==null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** No inputmatrix set\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }

  // check whether we have M and D matrices
  if (D_==null || M_==null)
  {
      cout << "***ERR*** MOERTEL::Manager::MakeSaddleProblem:\n"
           << "***ERR*** Matrix M or D is NULL\n"
           << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
      return NULL;
  }
  
  // create a matrix for the saddle problem and fill it
  saddlematrix_ = rcp(new Epetra_CrsMatrix(Copy,*saddlemap_,90));

  // add values from inputmatrix
  MatrixMatrixAdd(*inputmatrix_,false,1.0,*saddlematrix_,0.0);
  
  // add values from D_
  MatrixMatrixAdd(*D_,false,1.0,*saddlematrix_,1.0);
  MatrixMatrixAdd(*D_,true,1.0,*saddlematrix_,1.0);  

  // add values from M_
  MatrixMatrixAdd(*M_,false,1.0,*saddlematrix_,1.0);
  MatrixMatrixAdd(*M_,true,1.0,*saddlematrix_,1.0);  

  saddlematrix_->FillComplete();
  saddlematrix_->OptimizeStorage();

  return saddlematrix_.get();
}

/*----------------------------------------------------------------------*
 | Reset the internal solver (public)                        mwgee 12/05|
 *----------------------------------------------------------------------*/
void MOERTEL::Manager::ResetSolver()
{
  solverparams_ = null;
  solver_ = null;
  saddlematrix_ = null;
  return;
}

/*----------------------------------------------------------------------*
 | Set solver parameters (public)                            mwgee 12/05|
 *----------------------------------------------------------------------*/
void MOERTEL::Manager::SetSolverParameters(ParameterList& params)
{
  solverparams_ = rcp(new Teuchos::ParameterList(params));
  return;
}

/*----------------------------------------------------------------------*
 | solve problem (public)                                    mwgee 12/05|
 *----------------------------------------------------------------------*/
bool MOERTEL::Manager::Solve(ParameterList& params, Epetra_Vector& sol, const Epetra_Vector& rhs)
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
  if (solverparams_==null)
  {
    cout << "***ERR*** MOERTEL::Manager::Solve:\n"
         << "***ERR*** No solver parameters set, use SetSolverParameters(ParameterList& params)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  // test for problemmap_
  if (problemmap_==null)
  {
    cout << "***ERR*** MOERTEL::Manager::Solve:\n"
         << "***ERR*** No problem map set, use SetProblemMap(Epetra_Map* map)\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  // test for inputmatrix
  if (inputmatrix_==null)
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
  map<int,RefCountPtr<MOERTEL::Interface> >::iterator curr;
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
  if (D_==null || M_==null)
  {
    cout << "***ERR*** MOERTEL::Manager::Solve:\n"
         << "***ERR*** Matrix M or D is NULL\n"
         << "***ERR*** file/line: " << __FILE__ << "/" << __LINE__ << "\n";
    return false;
  }
  
  //---------------------------------------------------------------------------
  // make solution and rhs vector matching the system
  RefCountPtr<Epetra_Vector> b = rcp(const_cast<Epetra_Vector*>(&rhs)); 
  b.release();
  RefCountPtr<Epetra_Vector> x = rcp(&sol); 
  x.release();

  //---------------------------------------------------------------------------
  // get type of system to be used/generated
  RefCountPtr<Epetra_CrsMatrix> matrix = null;
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
  if (system=="SaddleSystem" || system=="saddlesystem" || system=="SADDLESYSTEM" || 
      system=="Saddle_System" || system=="saddle_system" || system=="SADDLE_SYSTEM")
  {
    if (saddlematrix_==null)
      MakeSaddleProblem();

    matrix = saddlematrix_;
    b = rcp(new Epetra_Vector(*saddlemap_,true)); 
    b.set_has_ownership();
    x = rcp(new Epetra_Vector(*saddlemap_,true)); 
    x.set_has_ownership();
    for (int i=0; i<rhs.MyLength(); ++i)
    {
      (*b)[i] = rhs[i];
      (*x)[i] = sol[i];
    }
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
  if (solver_==null)
    solver_ = rcp(new MOERTEL::Solver(Comm(),OutLevel()));
  
  //---------------------------------------------------------------------------
  // solve
  bool ok = solver_->Solve(solverparams_,matrix,x,b);
  
  //---------------------------------------------------------------------------
  // copy solution back to sol if neccessary
  if (x.has_ownership())
  {
    for (int i=0; i<sol.MyLength(); ++i)
      sol[i] = (*x)[i];
  }
  
  //---------------------------------------------------------------------------
  // output the parameter list to see whether something was unused  
  if (OutLevel()>8 && Comm().MyPID()==0)
  {
    cout << "MOERTEL: Solver Parameters:\n";
    solverparams_->print(cout);
  }

  return true;
}


