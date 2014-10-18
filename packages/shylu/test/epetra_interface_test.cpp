#include <assert.h>
#include <iostream>
#include <sstream>


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
#include "shylu_test_interface.hpp"
#include "have_interface.hpp"

using namespace std;

int main(int argc, char** argv)
{

#ifdef HAVE_MPI
  Teuchos::GlobalMPISession mpiSession(&argc, &argv, 0);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

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
      exit(1);
    }
  n = A->NumGlobalRows();

  //Make b vecotor
  
  Epetra_Map vecMap(n,0,Comm);
  b = new Epetra_MultiVector(vecMap,1,false);
  b->Random();
  x = new Epetra_MultiVector(vecMap,1,false);

  cout << "Epetra matrices loaded" << endl;


  /*-----------------have_interface-----------------*/
  /*---The have_interface checks is all the parameter list makes sense---*/
  Teuchos::RCP <Teuchos::ParameterList> pLUList;
  string pListFileName = "ShyLU_epetra_interface.xml";
  pLUList = Teuchos::getParametersFromXmlFile(pListFileName);

  if(myPID ==0)
    {
      cout << "Starting Have Check on parameter list" << endl;

      int one = test<double>(1.0);
      int interface_test = have_interface <Epetra_CrsMatrix>(*A, *pLUList);
      
    }

  /*----------------partitioning_interface--------------*/
  /*-----------Will use check the epetra matrix on partition_interface------*/
  
   int part_error = shylu_interface::partitioning_interface(A, AHat, b, bHat, *pLUList);
   cout << "partitioned done\n";
}
