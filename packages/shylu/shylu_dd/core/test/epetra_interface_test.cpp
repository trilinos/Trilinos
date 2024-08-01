// @HEADER
// *****************************************************************************
//               ShyLU: Scalable Hybrid LU Preconditioner and Solver
//
// Copyright 2011 NTESS and the ShyLU contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <assert.h>
#include <iostream>
#include <sstream>

//#include "Zoltan2_config.h"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#endif // HAVE_MPI
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
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

// Amesos includes
#include "Amesos.h"
#include "Amesos_BaseSolver.h"

// AztecOO includes
#include "AztecOO.h"

// ML includes
#include "ml_MultiLevelPreconditioner.h"
#include "ml_epetra_utils.h"

#include "Ifpack_ConfigDefs.h"
#include "shylu.h"
#include "shylu_util.h"
#include "Ifpack_ShyLU.h"
#include "Ifpack_ILUT.h"
#include "Ifpack_ILU.h"

#include "BelosConfigDefs.hpp"
#include "BelosLinearProblem.hpp"
#include "BelosEpetraAdapter.hpp"
#include "BelosBlockGmresSolMgr.hpp"

#include "EpetraExt_readEpetraLinearSystem.h"
#include "EpetraExt_CrsMatrixIn.h"

//Our test interfaces
#include "shylu.h"
#include "shylu_partition_interface.hpp"
#include "shylu_directsolver_interface.hpp"

using namespace std;


int main(int argc, char** argv)
{

#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, 0);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif


  bool success = true;
  string pass = "End Result: TEST PASSED";
  string fail = "End Result: TEST FAILED";

  int myPID = Comm.MyPID();
  if(myPID == 0)
    {
      cout << "Starting Epetra interface test" << endl;
    }

  /*----------------Load a test matrix---------------*/
  string matrixFileName = "wathenSmall.mtx";
  Epetra_CrsMatrix *A;
  Epetra_CrsMatrix *AHat;
  Epetra_MultiVector *b;
  Epetra_MultiVector *bHat;
  Epetra_MultiVector *x;
  int n = 0;

  //Get Matrix
  int err = EpetraExt::MatrixMarketFileToCrsMatrix(matrixFileName.c_str(), Comm, A);
  if(err!=0 && myPID ==0)
    {
      cout << "Error reading matrix file, info = " << err << endl;
      cout << fail << endl;
      exit(1);
    }
  n = A->NumGlobalRows();

  //Make b vecotor

  Epetra_Map vecMap(n,0,Comm);
  b = new Epetra_MultiVector(vecMap,1,false);
  b->Random();
  x = new Epetra_MultiVector(vecMap,1,false);
  x->Random();

  cout << "Epetra matrices loaded" << endl;


  /*-----------------have_interface-----------------*/
  /*---The have_interface checks is all the parameter list makes sense---*/
  Teuchos::RCP <Teuchos::ParameterList> pLUList;
  string pListFileName = "ShyLU_epetra_interface.xml";
  pLUList = Teuchos::getParametersFromXmlFile(pListFileName);


  /*----------------partitioning_interface--------------*/
  /*-----------Will use check the epetra matrix on partition_interface------*/


  //Isorropia Test - graph/Parmetis
  pLUList->set("Partitioning Package","Isorropia");
  Teuchos::ParameterList ptemp;
  ptemp = pLUList->sublist("Isorropia Input");

  Teuchos::ParameterList pptemp;
  pptemp = ptemp.sublist("Zoltan");
  pptemp.set("GRAPH_PACKAGE", "Parmetis");
  pptemp.set("DEBUG_LEVEL", "1");

  ptemp.set("partitioning method", "graph");
  ptemp.set("Zoltan", pptemp);
  pLUList->set("Isorropia Input", ptemp);

  cout << " \n\n--------------------BIG BREAK --------------\n\n";
  Teuchos::writeParameterListToXmlOStream(*pLUList, std::cout);


  ShyLU::PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector> partI(A, pLUList.get());
  partI.partition();
  AHat = partI.reorderMatrix();
  bHat = partI.reorderVector(b);

  EpetraExt::RowMatrixToMatlabFile("Epetra_Isorropia_Parmetis.mat", *AHat);


   cout << "Done with graph - parmetis (Zoltan)" << endl;

   /*

   //Isorropia Test - Graph/PT-Scotch
  pLUList->set("Partitioning Package","Isorropia");
  ptemp = pLUList->sublist("Isorropia Input");

  //Teuchos::ParameterList pptemp;
  pptemp = ptemp.sublist("Zoltan");
  pptemp.set("GRAPH_PACKAGE", "scotch");
  pptemp.set("DEBUG_LEVEL", "1");


  ptemp.set("partitioning method", "graph");
  ptemp.set("Zoltan", pptemp);
  pLUList->set("Isorropia Input", ptemp);

  cout << " \n\n--------------------BIG BREAK --------------\n\n";
  Teuchos::writeParameterListToXmlOStream(*pLUList, std::cout);

  PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector> partI2(A, pLUList.get());
  partI2.partition();
  AHat = partI2.reorderMatrix();
  bHat = partI2.reorderVector(b);
  cout << "Done with graph - pt-scotch" << endl;

   */

  //Zoltan2 Test

#if defined(HAVE_SHYLU_DDCORE_ZOLTAN2CORE) && (!defined(HAVE_SHYLU_DDCORE_TPETRA) || defined(HAVE_TPETRA_INST_INT_INT))

   //Isorropia Test - Graph/ParMetis
  pLUList->set("Partitioning Package","Zoltan2");
  ptemp = pLUList->sublist("Zoltan2 Input");
  ptemp.set("algorithm", "parmetis");
  ptemp.set("debug_level", "detailed_status");
  pLUList->set("Zoltan2 Input", ptemp);


  cout << " \n\n--------------------BIG BREAK --------------\n\n";
  Teuchos::writeParameterListToXmlOStream(*pLUList, std::cout);

  ShyLU::PartitionInterface<Epetra_CrsMatrix, Epetra_MultiVector> partI3(A, pLUList.get());
  partI3.partition();
  AHat = partI3.reorderMatrix();
  bHat = partI3.reorderVector(b);
  cout << "Done with graph - parmetis (Zoltan2)" << endl;

  EpetraExt::RowMatrixToMatlabFile("Epetra_Zoltan2_Parmetis.mat", *AHat);

#endif



  /*----------------------Direct Solver Interfaces----------------*/
  //#ifdef HAVE_SHYLU_AMESOS

  //Amesos - klu
  pLUList->set("Direct Solver Package", "Amesos");
  ptemp = pLUList->sublist("Amesos Input");
  pptemp = ptemp.sublist("Amesos_Klu Input");


  pptemp.set("PrintTiming", true);
  pptemp.set("PrintStatus", true);
  ptemp.set("Solver", "Amesos_Klu");
  ptemp.set("Amesos_Klu Input", pptemp);
  pLUList->set("Amesos Input", ptemp);


  cout << " \n\n--------------------BIG BREAK --------------\n\n";
  Teuchos::writeParameterListToXmlOStream(*pLUList, std::cout);

  ShyLU::DirectSolverInterface<Epetra_CrsMatrix, Epetra_MultiVector> directsolver(A, pLUList.get());

  directsolver.factor();
  directsolver.solve(b,x);

  cout << "Done with Amesos-KLU" << endl;

  //#endif

  //Amesos2 -klu2
#ifdef HAVE_SHYLU_AMESOS2

  pLUList->set("Direct Solver Package", "Amesos2");
  ptemp = pLUList->sublist("Amesos2 Input");
  //pptemp = ptemp.sublist("Amesos_Klu Input");


  pptemp.set("PrintTiming", true);
  pptemp.set("PrintStatus", true);
  ptemp.set("Solver", "KLU2");
  //ptemp.set("Amesos_Klu Input", pptemp);
  pLUList->set("Amesos2 Input", ptemp);


  cout << " \n\n--------------------BIG BREAK --------------\n\n";
  Teuchos::writeParameterListToXmlOStream(*pLUList, std::cout);


  

  ShyLU::DirectSolverInterface<Epetra_CrsMatrix, Epetra_MultiVector> directsolver2(A, pLUList.get());

  directsolver2.factor();
  directsolver2.solve(b,x);

  cout << "Done with Amesos-KLU2" << endl;

#endif

  if(success)
    {
      cout << pass << endl;
    }

  (void)bHat;
  return 0;
}

