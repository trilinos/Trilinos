#include <assert.h>
#include <iostream>
#include <sstream>

#include "shylu.h"
#include "shylu_partition_interface.hpp"
#include "ShyLU_config.h"

//Tperta
#ifdef HAVE_SHYLU_TPETRA
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
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_SerialComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Map.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"

// Teuchos includes
#ifdef HAVE_SHYLU_TPETRA
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

    typedef double scalar_type;
  typedef int local_o_type;
  typedef int global_o_type;
  typedef KokkosClassic::DefaultNode::DefaultNodeType node_type;
  typedef Tpetra::CrsMatrix<scalar_type, local_o_type, global_o_type, node_type> Matrix_t;
  typedef Tpetra::MultiVector<scalar_type, local_o_type, global_o_type, node_type> Vector_t;


  
  Teuchos::ParameterList defaultParameters;
  Teuchos::RCP <node_type> node = Teuchos::rcp(new node_type(defaultParameters));
  
  /*----------------Load a test matrix---------------*/
  string matrixFileName = "wathenSmall.mtx";
  
  //Get Matrix
  Teuchos::RCP<Matrix_t> A = Tpetra::MatrixMarket::Reader<Matrix_t>::readSparseFile(matrixFileName, comm, node);


  Teuchos::RCP<Vector_t> x = Teuchos::rcp(new Vector_t(A->getColMap(), 1));
  Teuchos::RCP<Vector_t> b = Teuchos::rcp(new Vector_t(A->getRowMap(), 1));
  b->randomize();
  x->randomize();

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
  
  
  cout << " \n\n--------------------BIG BREAK --------------\n\n";
  Teuchos::writeParameterListToXmlOStream(*pLUList, std::cout);


#ifdef HAVE_SHYLU_ZOLTAN2

  cout << "HSTER";

  ShyLU::PartitionInterface<Matrix_t, Vector_t> partI3(A.get(), pLUList.get());
  partI3.partition();
 
  cout << "Done with graph - parmetis" << endl;

#endif


#ifdef HAVE_SHYLU_AMESOS2
  
  pLUList->set("Direct Solver Package", "Amesos2");
  ptemp = pLUList->sublist("Amesos2 Input");
  //pptemp = ptemp.sublist("Amesos_Klu Input");


  //pptemp.set("PrintTiming", true);
  //pptemp.set("PrintStatus", true);
  ptemp.set("Solver", "SuperLU");
  //ptemp.set("Amesos_Klu Input", pptemp);
  pLUList->set("Amesos2 Input", ptemp);


  cout << " \n\n--------------------BIG BREAK --------------\n\n";
  Teuchos::writeParameterListToXmlOStream(*pLUList, std::cout);
  
  ShyLU::DirectSolverInterface<Matrix_t, Vector_t> directsolver2(A.get(), pLUList.get());
directsolver2.solve(b.get(),x.get());

  cout << "Done with Amesos-KLU2" << endl;
  
#endif

  
}

