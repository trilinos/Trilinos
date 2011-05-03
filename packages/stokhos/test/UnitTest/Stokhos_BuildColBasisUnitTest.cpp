#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_RCP.hpp>

#include "Stokhos_InterlacedTestSupport.hpp"

// Stokhos Stochastic Galerkin
#include "Stokhos_BasisInteractionGraph.hpp"
#include "Stokhos_EpetraSparse3Tensor.hpp"
#include "Stokhos_ParallelData.hpp"
#include "Stokhos_AdaptivityManager.hpp"
#include "Stokhos_AdaptivityUtils.hpp"

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
   #include "EpetraExt_MultiMpiComm.h"
   #include "mpi.h"
#else
   #include "Epetra_SerialComm.h"
   #include "EpetraExt_MultiSerialComm.h"
#endif
#include "Epetra_CrsGraph.h"

#include "EpetraExt_RowMatrixOut.h"

Teuchos::RCP<Epetra_CrsGraph> buildTridiagonalGraph(int numUnks,const Epetra_Comm & Comm)
{
   Epetra_Map map(-1,numUnks,0,Comm);
   Teuchos::RCP<Epetra_CrsGraph> graph
      = Teuchos::rcp(new Epetra_CrsGraph(Copy,map,0));
   
   // build tridiagonal graph
   int colCnt = 3;
   int * colPtr = 0;
   int colIndices[3];
   int colOffset[] = {-1, 0, 1};
   for(int myRow=0;myRow<numUnks;myRow++) {
      int row = map.GID(myRow);
      for(int i=0;i<3;i++)
         colIndices[i] = colOffset[i]+row;
      colCnt = 3;
      colPtr = colIndices;

      if(row==0) {
         colCnt = 2;
         colPtr = colIndices+1;
      }
      else if(row==map.NumGlobalElements()-1) 
         colCnt = 2;

      TEUCHOS_ASSERT(graph->InsertGlobalIndices(row,colCnt,colPtr)==0);
   }

   graph->FillComplete();

   return graph;
}

// Test construction of Linear Algebra (LA) data structures
TEUCHOS_UNIT_TEST(tBuildColBasis, test_adapted)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else 
      Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_SerialComm);
   #endif

   int numProc = comm->NumProc();
   int rank = comm->MyPID();

   out << "NumProc = " << numProc << ", Rank = " << rank << std::endl;

   int num_KL = 3;
   int porder = 3;

   Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = buildBasis(num_KL,porder);
   Teuchos::RCP<Epetra_CrsGraph> determGraph = buildTridiagonalGraph(3,*comm);

   std::vector<int> order(3);
   order[0] = 2; order[1] = 3; order[2] = 1;
   std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > sa_BasisPerDRow(3);
   sa_BasisPerDRow[0] = buildBasis(num_KL,1);
   sa_BasisPerDRow[1] = buildBasis(num_KL,1);
   sa_BasisPerDRow[2] = buildBasis(num_KL,order);
   
   std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > sa_BasisPerDCol;
   Stokhos::adapt_utils::buildColBasisFunctions(*determGraph,basis,sa_BasisPerDRow,sa_BasisPerDCol);

   if(numProc==2)
   {   TEST_EQUALITY(sa_BasisPerDCol.size(),4); }
   else 
   {   TEST_EQUALITY(sa_BasisPerDCol.size(),3); }

   for(std::size_t c=0;c<sa_BasisPerDCol.size();c++) {
      int gid = determGraph->ColMap().GID(c);
      Teuchos::RCP<const Stokhos::ProductBasis<int,double> > cBasis = 
         Teuchos::rcp_dynamic_cast<const Stokhos::ProductBasis<int,double> >(sa_BasisPerDCol[c]);

      Teuchos::Array<Teuchos::RCP<const Stokhos::OneDOrthogPolyBasis<int,double> > > cBases
         = cBasis->getCoordinateBases();
     
      for(std::size_t i=0;i<cBases.size();i++) {
         int bOrder = cBases[i]->order();

         if(gid==2)
         {   TEST_EQUALITY(bOrder,order[i]); }
         else if(gid==5)
         {   TEST_EQUALITY(bOrder,order[i]); }
         else
         {   TEST_EQUALITY(bOrder,1); }
      }
   }
}
