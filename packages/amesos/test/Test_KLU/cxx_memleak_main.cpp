#include "Amesos_ConfigDefs.h"

#ifdef HAVE_MPI
#include "mpi.h"
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_Util.h"
#include "Amesos_Klu.h"
#include "Amesos_TestRowMatrix.h"
#include "Teuchos_ParameterList.hpp"
#include "Galeri_Maps.h"
#include "Galeri_CrsMatrices.h"
#include "Galeri_Utils.h"

#include "Teuchos_ParameterList.hpp"

#include <vector>

using namespace Galeri;

//============ //
// main driver //
//============ //


Teuchos::RCP<Epetra_Map> buildMap(Epetra_Comm & Comm)
{
   std::vector<int> vec(4);
   int data[] = {0, 1 , 3 , 4};
  
   for(int i=0;i<4;i++)
      vec[i] = 6*Comm.MyPID()+data[i];

   return Teuchos::rcp(new Epetra_Map(-1,4,&vec[0],0,Comm));
}

Teuchos::RCP<Epetra_CrsMatrix> buildMatrix(Epetra_Map & map)
{
   Teuchos::RCP<Epetra_CrsMatrix> A = Teuchos::rcp(new Epetra_CrsMatrix(Copy,map,3));

   int ind[] = {0,0,0};
   double values[] = {-1,2,-1};

   int * indPtr;
   double * valuesPtr;
   int entries = 3;

   for(int i=0;i<4;i++) {
      int gid = map.GID(i);
      ind[0] = map.GID(i-1);
      ind[1] = map.GID(i);
      ind[2] = map.GID(i+1);

      indPtr = ind;
      valuesPtr = values;
      entries = 3;
      if(i==0) {
         entries = 2;
         indPtr = ind+1;
         valuesPtr = values+1;
      }
      else if(i==3) {
         entries = 2;
      }
      A->InsertGlobalValues(gid,entries,valuesPtr,indPtr);
   }
   A->FillComplete();

   return A;
}

int main(int argc, char *argv[]) 
{
#ifdef HAVE_MPI
  MPI_Init(&argc, &argv);
  Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm Comm;
#endif

  Teuchos::ParameterList GaleriList;

  Teuchos::RCP<Epetra_Map> Map = buildMap(Comm);
  Teuchos::RCP<Epetra_CrsMatrix> Matrix = buildMatrix(*Map);
 
  int NumVectors = 2;
  Epetra_MultiVector x(*Map,NumVectors);
  Epetra_MultiVector x_exact(*Map,NumVectors);
  Epetra_MultiVector b(*Map,NumVectors);
  x_exact.Random();
  Matrix->Apply(x_exact,b);

  // =========== //
  // AMESOS PART //
  // =========== //

  Epetra_LinearProblem Problem(&*Matrix, &x, &b);
  Amesos_Klu Solver(Problem);

  Teuchos::ParameterList List;
  List.set("Reindex", true);

  Solver.SetParameters(List); 

  AMESOS_CHK_ERR(Solver.Solve());
  AMESOS_CHK_ERR(Solver.Solve());

  double norm = ComputeNorm(&*Matrix, &x_exact, &b);
  if (Comm.MyPID() == 0)
    std::cout << "norm = " << norm << std::endl;

  if (norm > 1e-5)
    exit(EXIT_FAILURE);

#ifdef HAVE_MPI
  MPI_Finalize();
#endif

  return(EXIT_SUCCESS);
}
