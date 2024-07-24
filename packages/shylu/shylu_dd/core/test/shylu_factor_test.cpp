// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file shylu_sfactor_test.cpp

    \brief factor test

    \author Joshua Dennis Booth

    \remark Usage:
    \code mpirun -n np shylu_sfactor.exe

*/
// Will be used to test gFACT as we morph into templated


#include <assert.h>
#include <iostream>
#include <sstream>

#include "Isorropia_config.h" // Just for HAVE_MPI

// Epetra includes
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#endif // HAVE_MPI
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_LinearProblem.h"

// Teuchos includes
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"

// EpetraExt includes
#include "EpetraExt_RowMatrixOut.h"
#include "EpetraExt_MultiVectorOut.h"
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_MultiVectorIn.h"

#ifdef HAVE_SHYLU_DDCORE_TPETRA
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#endif

#include "shylu.h"
#include "shylu_util.h"

#include "EpetraExt_readEpetraLinearSystem.h"

using namespace std;

int main(int argc, char* argv[])
{
#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, 0);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  string pass = "End Result: TEST PASSED";
  string fail = "End Result: TEST FAILED";

  int myPID = Comm.MyPID();
  if(myPID == 0)
    {
      cout << "Starting SFactor Epetra Test" << endl;
    }

  //============================= Get Matrix
  string matrixFileName = "wathenSmall.mtx";
  Epetra_CrsMatrix *A;
  // Epetra_CrsMatrix *AHat;
  // int n = 0;

  int err = EpetraExt::MatrixMarketFileToCrsMatrix(matrixFileName.c_str(), Comm, A);

  if(err!=0 && myPID ==0)
    {
      cout << "Error reading matrix file, info = " << err << endl;
      cout << fail << endl;
      exit(1);
    }

  // not used
  // n = A->NumGlobalRows();

  //=============================Partition/Distribute
  Teuchos::ParameterList isoList;
  //isoList.set("partitioning method", "graph");

  cout << "before partition" << endl;
  Epetra_CrsMatrix *B = balanceAndRedistribute(A,isoList);
  cout << "after partition" << endl;

  shylu_data     slu_data_;
  shylu_config   slu_config_;
  shylu_symbolic slu_sym_;

  slu_config_.sym                 = 1;     //This is
  slu_config_.libName             = "Belos"; //This is
  slu_config_.schurSolver         = " "; //This is
  slu_config_.schurAmesosSolver   = " " ; //This is
  slu_config_.diagonalBlockSolver = "Amesos_Klu"; //This is
  slu_config_.relative_threshold  = 0.0; //This is
  slu_config_.Sdiagfactor         = 0.05; //This is
  slu_config_.schurApproxMethod   = 2;  //1 A22Block  2 Thresholding 3 Guided
  slu_config_.schurPreconditioner = "ILU stand-alone"; //This is
  slu_config_.silent_subiter      = true; //This is
  slu_config_.inner_tolerance     = 1e-5; //This is
  slu_config_.inner_maxiters      = 100; //This is
  slu_config_.overlap             = 0; //This is
  slu_config_.sep_type            = 1; //1 Wide 2 Narrow
  slu_config_.amesosForDiagonal   = true;


  int serr = shylu_symbolic_factor(B, &slu_sym_, &slu_data_, &slu_config_);
  cout << "shylu_symbolic_factor done:" << endl;
  cout << "Return value: " << serr << endl;
  if(serr == 0)
    cout << pass << endl;
  return serr;
}//end main
