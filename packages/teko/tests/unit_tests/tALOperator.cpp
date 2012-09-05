/*
 * Author: Zhen Wang
 * Email: wangz@ornl.gov
 *        zhen.wang@alum.emory.edu
 */

/*
 * This tests reads tOpMat.mm (Saddle point matrix),
 * tOpMp.mm (pressure mass matrix) and
 * tOpRhs.mm (right hand side we should get)
 * and test Teko_ALOperator.
 */

#include <iostream>
#include <fstream>
#include <cmath>

// Teuchos
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_UnitTestHarness.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"

// Epetra
#ifdef HAVE_MPI
#include "Epetra_MpiComm.h"
#include "mpi.h"
#else
#include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

// EpetraExt
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// Thyra
#include "Thyra_EpetraLinearOp.hpp"

// Teko
#include "Teko_ALOperator.hpp"

using namespace Teko;
using namespace Teko::Epetra;

// int
// main(int argc, char * argv[])
TEUCHOS_UNIT_TEST(tALOperator, test)
{
   // Build communicator
#ifdef HAVE_MPI
   Epetra_MpiComm Comm(MPI_COMM_WORLD);
#else
   Epetra_SerialComm Comm;
#endif

   // Get process information
   int myPID = Comm.MyPID();
   out << "MPI_PID = " << myPID << ", UNIX_PID = " << getpid() << std::endl;

   // Maps.
   int dim = 2, numVel = 3, numPre = 2, errCode;
   Epetra_Map mapVel(numVel, 0, Comm), mapPre(numPre, 0, Comm);
   Epetra_Map mapAll(numVel * dim + numPre, 0, Comm);

   // Reorder.
   std::vector<int> reorderedVec;
   int numMyLen = mapVel.NumMyElements();
   int * myGlb;
   myGlb = mapVel.MyGlobalElements();
   for(int i = 0; i < dim; i++)
      for(int j = 0; j < numMyLen; j++)
         reorderedVec.push_back(myGlb[j] + numVel * i);
   numMyLen = mapPre.NumMyElements();
   myGlb = mapPre.MyGlobalElements();
   for(int j = 0; j < numMyLen; j++)
      reorderedVec.push_back(myGlb[j] + numVel * dim);

   Teuchos::RCP < Epetra_Map > mapReorder = Teuchos::rcp(
         new Epetra_Map(-1, reorderedVec.size(), &reorderedVec[0], 0, Comm));
   Teuchos::RCP < Epetra_Import > importReorder = Teuchos::rcp(new Epetra_Import(*mapReorder, mapAll));

   std::vector<std::vector<int> > blockedVec;
   numMyLen = mapVel.NumMyElements();
   myGlb = mapVel.MyGlobalElements();
   for(int i = 0; i < dim; i++)
   {
      reorderedVec.clear();
      for(int j = 0; j < numMyLen; j++)
         reorderedVec.push_back(myGlb[j] + numVel * i);
      blockedVec.push_back(reorderedVec);
   }
   numMyLen = mapPre.NumMyElements();
   myGlb = mapPre.MyGlobalElements();
   reorderedVec.clear();
   for(int j = 0; j < numMyLen; j++)
      reorderedVec.push_back(myGlb[j] + numVel * dim);
   blockedVec.push_back(reorderedVec);

   // Read matrices and vector.
   Epetra_CrsMatrix *ptrMat = 0, *ptrMp = 0;
   TEUCHOS_ASSERT(EpetraExt::MatrixMarketFileToCrsMatrix("data/tOpMat.mm", mapAll, ptrMat)==0);
   TEUCHOS_ASSERT(EpetraExt::MatrixMarketFileToCrsMatrix("data/tOpMp.mm", mapPre, ptrMp)==0);
   LinearOp lpMp = Thyra::epetraLinearOp(Teuchos::rcpFromRef(*ptrMp));
   // This vector is computed by Matlab for comparison.
   Epetra_Vector *ptrExact = 0;
   TEUCHOS_ASSERT(EpetraExt::MatrixMarketFileToVector("data/tOpRhs.mm", mapAll, ptrExact)==0);

   // Reorder matrix.
   RCP < Epetra_CrsMatrix > mat = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *mapReorder, ptrMat->GlobalMaxNumEntries()));
   errCode = mat->Import(*ptrMat, *importReorder, Insert);
   errCode = mat->FillComplete();

   // Build augmented Lagrangian-based operator.
   Teko::NS::ALOperator al(blockedVec, mat, lpMp);

   // Initialize vectors.
   Epetra_Vector x(*mapReorder, false), b(*mapReorder, false);
   x.PutScalar(1.0);
   b.PutScalar(0.0);

   // Apply operator.
   al.Apply(x, b);

   // Compare computed vector and exact vector.
   b.Update(-1.0, *ptrExact, 1.0);
   double norm2;
   b.Norm2(&norm2);
   if(norm2 < 1.0e-15)
   {
      out << "Test:ALOperator: Passed." << std::endl;
      errCode = 0;
   }
   else
   {
      out << "Test:ALOperator: Failed." << std::endl;
      errCode = -1;
   }

   delete ptrMat;
   delete ptrMp;
   delete ptrExact;

   TEST_ASSERT(errCode==0);
}
