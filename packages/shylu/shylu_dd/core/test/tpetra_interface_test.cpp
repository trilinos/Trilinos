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

#include "shylu.h"
#include "shylu_partition_interface.hpp"
#include "ShyLU_DDCore_config.h"

//Tperta
#ifdef HAVE_SHYLU_DDCORE_TPETRA
#include <Tpetra_Core.hpp>
#include <Tpetra_Version.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Map.hpp>
#include <MatrixMarket_Tpetra.hpp>
#endif


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
#ifdef HAVE_SHYLU_DDCORE_TPETRA
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#endif

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

  Tpetra::ScopeGuard mpiSession(&argc, &argv);
  Teuchos::RCP <const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int myPID = comm->getRank();

  bool success = true;
  string pass = "End Result: TEST PASSED";
  string fail = "End Result: TEST FAILED";

#if !defined(HAVE_TPETRA_INST_INT_INT)
  if(myPID == 0)
    {
      cout << "Tpetra was not instantiated with INT_INT" << endl;
    }
#else
  typedef double scalar_type;
  typedef int local_o_type;
  typedef int global_o_type;
  //typedef Tpetra::KokkosClassic::DefaultNode::DefaultNodeType node_type;

  typedef Tpetra::Details::DefaultTypes::node_type node_type;

  typedef Tpetra::CrsMatrix<scalar_type, local_o_type, global_o_type, node_type> Matrix_t;
  typedef Tpetra::MultiVector<scalar_type, local_o_type, global_o_type, node_type> Vector_t;

  if(myPID == 0)
    {
      cout << "Starting Tpetra interface test" << endl;
    }

  Teuchos::ParameterList defaultParameters;

  /*----------------Load a test matrix---------------*/
  string matrixFileName = "wathenSmall.mtx";

  //Get Matrix
  Teuchos::RCP<Matrix_t> A = Tpetra::MatrixMarket::Reader<Matrix_t>::readSparseFile(matrixFileName, comm);


  //Note:: Tpetra::MatrixMarket::Reader is providing A->getColMap() wrong and equal A->row
  Teuchos::RCP<Vector_t> x = Teuchos::rcp(new Vector_t(A->getRowMap(), 1));
  Teuchos::RCP<Vector_t> b = Teuchos::rcp(new Vector_t(A->getRowMap(), 1));
  b->randomize();
  x->randomize();

  if(myPID == 0)
    {
      cout << "num_vector     : " << b->getNumVectors() << " "
           << x->getNumVectors() << endl;
      cout << "length         : " << b->getGlobalLength() << " "
           << x->getGlobalLength() << endl;

      cout << "A length       : " << A->getGlobalNumRows() << " " << A->getGlobalNumCols() << endl;
      cout << "A local length : " << A->getLocalNumRows() << " " << A->getLocalNumCols() << endl;
    }

  /*-----------------have_interface-----------------*/
  /*---The have_interface checks is all the parameter list makes sense---*/
  Teuchos::RCP <Teuchos::ParameterList> pLUList;
  string pListFileName = "ShyLU_epetra_interface.xml";
  pLUList = Teuchos::getParametersFromXmlFile(pListFileName);


  /*----------------partitioning_interface--------------*/
  /*-----------Will use check the epetra matrix on partition_interface------*/

  pLUList->set("Partitioning Package","Zoltan2");
  Teuchos::ParameterList ptemp = pLUList->sublist("Zoltan2 Input");
  ptemp.set("algorithm", "parmetis");
  ptemp.set("debug_level", "detailed_status");
  pLUList->set("Zoltan2 Input", ptemp);


  if(myPID == 0)
    {
      cout << " \n\n--------------------BIG BREAK --------------\n\n";
      Teuchos::writeParameterListToXmlOStream(*pLUList, std::cout);
    }


#if defined(HAVE_SHYLU_DDCORE_ZOLTAN2CORE)
  ShyLU::PartitionInterface<Matrix_t, Vector_t> partI3(A.get(), pLUList.get());
  partI3.partition();

  if(myPID == 0)
    {
      cout << "Done with graph - parmetis" << endl;
    }
#else

  success = false;

#endif


#ifdef HAVE_SHYLU_DDCORE_AMESOS2

  pLUList->set("Direct Solver Package", "Amesos2");
  ptemp = pLUList->sublist("Amesos2 Input");
  //pptemp = ptemp.sublist("Amesos_Klu Input");

  //pptemp.set("PrintTiming", true);
  //pptemp.set("PrintStatus", true);
  ptemp.set("Solver", "SuperLU");
  //ptemp.set("Amesos_Klu Input", pptemp);
  pLUList->set("Amesos2 Input", ptemp);


  if(myPID == 0)
    {
      cout << " \n\n--------------------BIG BREAK --------------\n\n";
      Teuchos::writeParameterListToXmlOStream(*pLUList, std::cout);

      cout << "num_vector: " << b->getNumVectors() << " "
           << x->getNumVectors() << endl;
      cout << "length: " << b->getGlobalLength() << " "
           << x->getGlobalLength() << endl;

      cout << "A length" << A->getGlobalNumRows() << " " << A->getGlobalNumCols() << endl;
      cout << "A local length" << A->getLocalNumRows() << " " << A->getLocalNumCols() << endl;
    }


  ShyLU::DirectSolverInterface<Matrix_t, Vector_t> directsolver2(A.get(), pLUList.get());
  directsolver2.factor();
  directsolver2.solve(b.get(),x.get());

//Note: should multiple to set b and x for success

  if(myPID == 0)
    {
      cout << "Done with Amesos2-SuperLU" << endl;
    }
#else

  success = false;

#endif
#endif

  if(myPID == 0)
    {
      if(success)
        cout << pass << endl;
      else
        cout << fail << endl;
    }

  return 0;
}

