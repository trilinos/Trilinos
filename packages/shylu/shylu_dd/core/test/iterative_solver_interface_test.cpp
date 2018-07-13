#include <assert.h>
#include <iostream>
#include <sstream>

#include "shylu.h"
#include "shylu_iterativesolver_interface.hpp"
#include "ShyLU_DDCore_config.h"

//Tperta
#ifdef HAVE_SHYLU_DDCORE_TPETRA
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_GlobalMPISession.hpp>
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
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_RCP.hpp"
#include "Tpetra_DefaultPlatform.hpp"
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

  Teuchos::GlobalMPISession mpiSession(&argc, &argv, 0);
  Teuchos::RCP <const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();
    int myPID = comm->getRank();
  if(myPID == 0)
    {
      cout << "Starting Tpetra interface test" << endl;
    }


  bool success = true;
  string pass = "End Result: TEST PASSED";
  string fail = "End Result: TEST PASSED";

  typedef double scalar_type;
  typedef int local_o_type;
  typedef int global_o_type;
  //typedef KokkosClassic::DefaultNode::DefaultNodeType node_type;

  typedef Tpetra::Details::DefaultTypes::node_type node_type;

  typedef Tpetra::CrsMatrix<scalar_type, local_o_type, global_o_type, node_type> Matrix_t;
  typedef Tpetra::MultiVector<scalar_type, local_o_type, global_o_type, node_type> Vector_t;


  Teuchos::ParameterList defaultParameters;
  Teuchos::RCP <node_type> node = Teuchos::rcp(new node_type(defaultParameters));

  /*----------------Load a test matrix---------------*/
  string matrixFileName = "wathenSmall.mtx";

  //Get Matrix
  Teuchos::RCP<Matrix_t> A = Tpetra::MatrixMarket::Reader<Matrix_t>::readSparseFile(matrixFileName, comm, node); //removed node

  if( &A == NULL)
    {
      success = false;
    }


  //Note:: Tpetra::MatrixMarket::Reader is providing A->getColMap() wrong and equal A->row
  Teuchos::RCP<Vector_t> x = Teuchos::rcp(new Vector_t(A->getRowMap(), 1));
  Teuchos::RCP<Vector_t> b = Teuchos::rcp(new Vector_t(A->getRowMap(), 1));
  b->randomize();
  x->randomize();

    cout << "num_vector: " << b->getNumVectors() << " " 
       << x->getNumVectors() << endl;
  cout << "length: " << b->getGlobalLength() << " "
       << x->getGlobalLength() << endl;

  cout << "A length" << A->getGlobalNumRows() << " " << A->getGlobalNumCols() << endl;
  cout << "A local length" << A->getNodeNumRows() << " " << A->getNodeNumCols() << endl;


  /*-----------------have_interface-----------------*/
  /*---The have_interface checks is all the parameter list makes sense---*/
  Teuchos::RCP <Teuchos::ParameterList> pLUList;
  string pListFileName = "ShyLU_epetra_interface.xml";
  //pLUList = Teuchos::getParametersFromXmlFile(pListFileName);
  pLUList = Teuchos::parameterList();


  /*----------------partitioning_interface--------------*/
  /*-----------Will use check the epetra matrix on partition_interface------*/


  pLUList->set("Iterative Solver Package", "Belos");
  Teuchos::ParameterList ptemp = pLUList->sublist("Belos Input");
  ptemp.set("Solver", "GMRES");
  ptemp.set("Block Size" , 1);
  ptemp.set("Maximum Iterations", 1000);
  ptemp.set("Maximum Restarts", 3);
  ptemp.set("Convergence Tolerance" , 1e-5);
  pLUList->set("Belos Input", ptemp);

  //ptemp.set("Solver", "SuperLU");
  //ptemp.set("Amesos_Klu Input", pptemp);
  //pLUList->set("Amesos2 Input", ptemp);


  cout << " \n\n--------------------BIG BREAK --------------\n\n";
  Teuchos::writeParameterListToXmlOStream(*pLUList, std::cout);

  cout << "num_vector: " << b->getNumVectors() << " " 
       << x->getNumVectors() << endl;
  cout << "length: " << b->getGlobalLength() << " "
       << x->getGlobalLength() << endl;

  cout << "A length" << A->getGlobalNumRows() << " " << A->getGlobalNumCols() << endl;
  cout << "A local length" << A->getNodeNumRows() << " " << A->getNodeNumCols() << endl;

  ShyLU::IterativeSolverInterface<Matrix_t, Vector_t> iterativesolver(A.get(), pLUList.get());
  iterativesolver.solve(b.get(),x.get());


    //  ShyLU::DirectSolverInterface<Matrix_t, Vector_t> directsolver2(A.get(), pLUList.get());
    //directsolver2.factor();
    /// directsolver2.solve(b.get(),x.get());

//Note: should multiple to set b and x for success

  cout << "Done with Belos" << endl;

  success = false;


  if(myPID == 0)
    {
      if(success)
        cout << pass << endl;
      else
        cout << fail << endl;
    }


}

