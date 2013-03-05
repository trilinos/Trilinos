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

Teuchos::RCP<Epetra_CrsMatrix> buildTridiagonalOp(const Epetra_CrsGraph & graph,double * stencil)
{
   Teuchos::RCP<Epetra_CrsMatrix> result = Teuchos::rcp(new Epetra_CrsMatrix(Copy,graph));
   const Epetra_Map & map = result->RowMap();

   result->PutScalar(0.0);

   // build tridiagonal graph
   int colCnt = 3;
   int * colPtr = 0;
   int colIndices[3];
   int colOffset[] = {-1, 0, 1};
   double * stencilPtr = 0;
   for(int myRow=0;myRow<map.NumMyElements();myRow++) {
      int row = map.GID(myRow);
      for(int i=0;i<3;i++)
         colIndices[i] = colOffset[i]+row;
      colCnt = 3;
      colPtr = colIndices;
      stencilPtr = stencil;

      if(row==0) {
         colCnt = 2;
         colPtr = colIndices+1;
         stencilPtr = stencil+1;
      }
      else if(row==map.NumGlobalElements()-1) 
         colCnt = 2;

      TEUCHOS_ASSERT(result->SumIntoGlobalValues(row,colCnt,stencilPtr,colPtr)==0);
   }

   return result;
}

// some testing infrastructure
template <typename ScalarT,typename OrdinalT>
bool ord_func(std::pair<ScalarT,OrdinalT> const & a, std::pair<ScalarT,OrdinalT> const & b) 
{ return a.first < b.first; }

template <typename ScalarT>
void generate_sorted_order(const std::vector<ScalarT> & values,std::vector<std::size_t> & order)
{
   typedef std::pair<ScalarT,std::size_t> Pair;

   // build vector to sort
   std::vector<Pair> pairValues(values.size());
   for(std::size_t i=0;i<values.size();i++) 
      pairValues[i] = std::make_pair(values[i],i);
  
   // build sorted ordering
   std::sort(pairValues.begin(),pairValues.end(),ord_func<ScalarT,std::size_t>);

   // write out new order
   order.resize(pairValues.size());
   for(std::size_t i=0;i<pairValues.size();i++) 
      order[i] = pairValues[i].second;
}

template <typename ScalarT>
void apply_ordering(std::vector<ScalarT> & values, const std::vector<std::size_t> & order)
{
   typedef std::pair<std::size_t,ScalarT> Pair;

   // build vector to sort
   std::vector<Pair> pairValues(values.size());
   for(std::size_t i=0;i<order.size();i++) 
      pairValues[i] = std::make_pair(order[i],values[i]);
  
   // build sorted ordering
   std::sort(pairValues.begin(),pairValues.end(),ord_func<std::size_t,ScalarT>);

   // write out new values
   for(std::size_t i=0;i<pairValues.size();i++) 
      values[i] = pairValues[i].second;
}

// Test construction of Linear Algebra (LA) data structures
TEUCHOS_UNIT_TEST(tAdaptivityManager, test_interface)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<const Epetra_MpiComm> mpiComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
      Teuchos::RCP<const EpetraExt::MultiComm> multiComm = Teuchos::rcp(new EpetraExt::MultiMpiComm(*mpiComm,-1));
      Teuchos::RCP<const Epetra_Comm> comm = mpiComm;
   #else 
      Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_SerialComm);
      Teuchos::RCP<const EpetraExt::MultiComm> multiComm = Teuchos::rcp(new EpetraExt::MultiSerialComm(-1));
   #endif

   int num_KL = 3;
   int porder = 3;

   Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = buildBasis(num_KL,porder);
   Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk = basis->computeTripleProductTensor();

   std::vector<int> order(3);
   order[0] = 2; order[1] = 3; order[2] = 1;
   std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > sa_BasisPerDRow(3);
   sa_BasisPerDRow[0] = buildBasis(num_KL,1);
   sa_BasisPerDRow[1] = buildBasis(num_KL,1);
   sa_BasisPerDRow[2] = buildBasis(num_KL,order);

   Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);

   Teuchos::RCP<Epetra_CrsGraph> determGraph = buildTridiagonalGraph(3,*comm);
   Stokhos::AdaptivityManager adaptMngr(basis,sa_BasisPerDRow,*determGraph,false,-1);

   TEST_EQUALITY(adaptMngr.getGlobalRowId(0,0),0);
   TEST_EQUALITY(adaptMngr.getGlobalRowId(1,0),sa_BasisPerDRow[0]->size());

   TEST_EQUALITY(adaptMngr.getGlobalColId(0,0),0);
   TEST_EQUALITY(adaptMngr.getGlobalColId(1,0),sa_BasisPerDRow[0]->size());

   TEST_EQUALITY(adaptMngr.getRowStochasticBasisSize(0),sa_BasisPerDRow[0]->size());
   TEST_EQUALITY(adaptMngr.getRowStochasticBasisSize(1),sa_BasisPerDRow[1]->size());
   TEST_EQUALITY(adaptMngr.getRowStochasticBasisSize(2),sa_BasisPerDRow[2]->size());

   TEST_EQUALITY(adaptMngr.getColStochasticBasisSize(0),sa_BasisPerDRow[0]->size());
   TEST_EQUALITY(adaptMngr.getColStochasticBasisSize(1),sa_BasisPerDRow[1]->size());
   TEST_EQUALITY(adaptMngr.getColStochasticBasisSize(2),sa_BasisPerDRow[2]->size());
}

// Test construction of Linear Algebra (LA) data structures
TEUCHOS_UNIT_TEST(tAdaptivityManager, sum_in_op_eq_order)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<const Epetra_MpiComm> mpiComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
      Teuchos::RCP<const EpetraExt::MultiComm> multiComm = Teuchos::rcp(new EpetraExt::MultiMpiComm(*mpiComm,-1));
      Teuchos::RCP<const Epetra_Comm> comm = mpiComm;
   #else 
      Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_SerialComm);
      Teuchos::RCP<const EpetraExt::MultiComm> multiComm = Teuchos::rcp(new EpetraExt::MultiSerialComm(-1));
   #endif

   int num_KL = 2;
   int porder = 2;
   int numDetermUnks = 3;

   double stencil[]  = {-1.0,2.0,-1.0};
   Teuchos::RCP<Epetra_CrsGraph> determGraph = buildTridiagonalGraph(numDetermUnks,*comm);
   Teuchos::RCP<Epetra_CrsMatrix> determOp = buildTridiagonalOp(*determGraph,stencil);

   Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = buildBasis(num_KL,porder);
   Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk = basis->computeTripleProductTensor();
   Teuchos::RCP<Stokhos::EpetraSparse3Tensor> epetraCijk =
         Teuchos::rcp(new Stokhos::EpetraSparse3Tensor(basis,Cijk,multiComm));

   std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > sa_BasisPerDRow(numDetermUnks);
   sa_BasisPerDRow[0] = buildBasis(num_KL,1);
   sa_BasisPerDRow[1] = buildBasis(num_KL,1);
   sa_BasisPerDRow[2] = buildBasis(num_KL,1);

   {
      Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
      params->set("Scale Operator by Inverse Basis Norms",true);
      Stokhos::AdaptivityManager adaptMngr(basis,sa_BasisPerDRow,*determGraph,false,-1,true);

      Teuchos::RCP<Epetra_CrsMatrix> sg_A = adaptMngr.buildMatrixFromGraph();
   
      TEST_EQUALITY(sg_A->NumGlobalRows(),sg_A->NumGlobalCols()); // check for square
      TEST_EQUALITY(sg_A->NumGlobalRows(),3*numDetermUnks);     // check sizing
   
      adaptMngr.sumInOperator(*sg_A,*Cijk,0,*determOp);
   
      int determDof = 1;
      // choose second determrow
      for(int stochBasis=0;stochBasis<sa_BasisPerDRow[determDof]->size();stochBasis++) {
         int numEntries = 0;
         std::vector<std::size_t> order;
         std::vector<int> stor_indices(9), indices;
         std::vector<double> stor_values(9), values;
             
         sg_A->ExtractGlobalRowCopy(adaptMngr.getGlobalRowId(determDof,stochBasis),9,numEntries,&stor_values[0],&stor_indices[0]);
      
         TEST_EQUALITY(numEntries,9);
      
         // perform copy so things are correct size
         for(int i=0;i<numEntries;i++) {
            indices.push_back(stor_indices[i]);
            values.push_back(stor_values[i]);
         }
      
         // order indices and values
         generate_sorted_order(indices,order);     
         apply_ordering(indices,order);     
         apply_ordering(values,order);     
      
         int rowTerm = basis->index(sa_BasisPerDRow[determDof]->term(stochBasis));
         int colTerm0 = basis->index(adaptMngr.getColStochasticBasis(determDof)->term(0));
         int colTerm1 = basis->index(adaptMngr.getColStochasticBasis(determDof)->term(1));
         int colTerm2 = basis->index(adaptMngr.getColStochasticBasis(determDof)->term(2));
         double normValue = basis->norm_squared(rowTerm);
      
         // test middle row values
         TEST_EQUALITY(values[0],stencil[0]*Cijk->getValue(rowTerm,colTerm0,0)/normValue);
         TEST_EQUALITY(values[1],stencil[0]*Cijk->getValue(rowTerm,colTerm1,0)/normValue);
         TEST_EQUALITY(values[2],stencil[0]*Cijk->getValue(rowTerm,colTerm2,0)/normValue);
      
         TEST_EQUALITY(values[3],stencil[1]*Cijk->getValue(rowTerm,colTerm0,0)/normValue);
         TEST_EQUALITY(values[4],stencil[1]*Cijk->getValue(rowTerm,colTerm1,0)/normValue);
         TEST_EQUALITY(values[5],stencil[1]*Cijk->getValue(rowTerm,colTerm2,0)/normValue);
      
         TEST_EQUALITY(values[6],stencil[2]*Cijk->getValue(rowTerm,colTerm0,0)/normValue);
         TEST_EQUALITY(values[7],stencil[2]*Cijk->getValue(rowTerm,colTerm1,0)/normValue);
         TEST_EQUALITY(values[8],stencil[2]*Cijk->getValue(rowTerm,colTerm2,0)/normValue);
      }

      TEST_ASSERT(sg_A->Filled());
   }

   {
      Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
      params->set("Scale Operator by Inverse Basis Norms",false);
      Stokhos::AdaptivityManager adaptMngr(basis,sa_BasisPerDRow,*determGraph,false,-1,false);

      Teuchos::RCP<Epetra_CrsMatrix> sg_A = adaptMngr.buildMatrixFromGraph();
   
      TEST_EQUALITY(sg_A->NumGlobalRows(),sg_A->NumGlobalCols()); // check for square
      TEST_EQUALITY(sg_A->NumGlobalRows(),3*numDetermUnks);     // check sizing
   
      adaptMngr.sumInOperator(*sg_A,*Cijk,0,*determOp);
   
      int determDof = 1;
      // choose second determrow
      for(int stochBasis=0;stochBasis<sa_BasisPerDRow[determDof]->size();stochBasis++) {
         int numEntries = 0;
         std::vector<std::size_t> order;
         std::vector<int> stor_indices(9), indices;
         std::vector<double> stor_values(9), values;
             
         sg_A->ExtractGlobalRowCopy(adaptMngr.getGlobalRowId(determDof,stochBasis),9,numEntries,&stor_values[0],&stor_indices[0]);
      
         TEST_EQUALITY(numEntries,9);
      
         // perform copy so things are correct size
         for(int i=0;i<numEntries;i++) {
            indices.push_back(stor_indices[i]);
            values.push_back(stor_values[i]);
         }
      
         // order indices and values
         generate_sorted_order(indices,order);     
         apply_ordering(indices,order);     
         apply_ordering(values,order);     
      
         int rowTerm = basis->index(sa_BasisPerDRow[determDof]->term(stochBasis));
         int colTerm0 = basis->index(adaptMngr.getColStochasticBasis(determDof)->term(0));
         int colTerm1 = basis->index(adaptMngr.getColStochasticBasis(determDof)->term(1));
         int colTerm2 = basis->index(adaptMngr.getColStochasticBasis(determDof)->term(2));
      
         // test middle row values
         TEST_EQUALITY(values[0],stencil[0]*Cijk->getValue(rowTerm,colTerm0,0));
         TEST_EQUALITY(values[1],stencil[0]*Cijk->getValue(rowTerm,colTerm1,0));
         TEST_EQUALITY(values[2],stencil[0]*Cijk->getValue(rowTerm,colTerm2,0));
      
         TEST_EQUALITY(values[3],stencil[1]*Cijk->getValue(rowTerm,colTerm0,0));
         TEST_EQUALITY(values[4],stencil[1]*Cijk->getValue(rowTerm,colTerm1,0));
         TEST_EQUALITY(values[5],stencil[1]*Cijk->getValue(rowTerm,colTerm2,0));
      
         TEST_EQUALITY(values[6],stencil[2]*Cijk->getValue(rowTerm,colTerm0,0));
         TEST_EQUALITY(values[7],stencil[2]*Cijk->getValue(rowTerm,colTerm1,0));
         TEST_EQUALITY(values[8],stencil[2]*Cijk->getValue(rowTerm,colTerm2,0));
      }

      TEST_ASSERT(sg_A->Filled());
   }
}

TEUCHOS_UNIT_TEST(tAdaptivityManager, sum_in_op_var_order)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<const Epetra_MpiComm> mpiComm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
      Teuchos::RCP<const EpetraExt::MultiComm> multiComm = Teuchos::rcp(new EpetraExt::MultiMpiComm(*mpiComm,-1));
      Teuchos::RCP<const Epetra_Comm> comm = mpiComm;
   #else 
      Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_SerialComm);
      Teuchos::RCP<const EpetraExt::MultiComm> multiComm = Teuchos::rcp(new EpetraExt::MultiSerialComm(-1));
   #endif

   int num_KL = 4;
   int porder = 3;
   int numDetermUnks = 3;

   double stencil[]  = {-1.0,2.0,-1.0};
   Teuchos::RCP<Epetra_CrsGraph> determGraph = buildTridiagonalGraph(numDetermUnks,*comm);
   Teuchos::RCP<Epetra_CrsMatrix> determOp = buildTridiagonalOp(*determGraph,stencil);

   Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = buildBasis(num_KL,porder);
   Teuchos::RCP<Stokhos::Sparse3Tensor<int,double> > Cijk = basis->computeTripleProductTensor();
   Teuchos::RCP<Stokhos::EpetraSparse3Tensor> epetraCijk =
         Teuchos::rcp(new Stokhos::EpetraSparse3Tensor(basis,Cijk,multiComm));

   Stokhos::BasisInteractionGraph big(*basis);

   std::vector<int> vorder(4);
   vorder[0] = 2; vorder[1] = 3; vorder[2] = 2; vorder[3] = 0;
   std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > sa_BasisPerDRow(numDetermUnks);
   sa_BasisPerDRow[0] = buildBasis(num_KL,1);
   sa_BasisPerDRow[1] = buildBasis(num_KL,2);
   sa_BasisPerDRow[2] = buildBasis(num_KL,vorder);

   {
      Teuchos::RCP<Teuchos::ParameterList> params = Teuchos::rcp(new Teuchos::ParameterList);
      params->set("Scale Operator by Inverse Basis Norms",false);
      Stokhos::AdaptivityManager adaptMngr(basis,sa_BasisPerDRow,*determGraph,false,-1,false);

      Teuchos::RCP<Epetra_CrsMatrix> sg_A = adaptMngr.buildMatrixFromGraph();
   
      TEST_EQUALITY(sg_A->NumGlobalRows(),sg_A->NumGlobalCols()); // check for square
   
      adaptMngr.sumInOperator(*sg_A,*Cijk,0,*determOp);

      out << "Summed into" << std::endl;
   
      int determDof = 1;
      // choose second determrow
      for(int stochBasis=0;stochBasis<sa_BasisPerDRow[determDof]->size();stochBasis++) {
         int numEntries = 0;
         std::vector<std::size_t> order;
         std::vector<int> stor_indices(400), indices;
         std::vector<double> stor_values(400), values;
             
         out << "grabbing row " << stochBasis << " values" << std::endl; 
         TEST_ASSERT(sg_A->ExtractGlobalRowCopy(adaptMngr.getGlobalRowId(determDof,stochBasis),400,numEntries,&stor_values[0],&stor_indices[0])==0);
         out << "num entries " << numEntries << std::endl;
      
         TEST_ASSERT(numEntries<400);
      
         // perform copy so things are correct size
         for(int i=0;i<numEntries;i++) {
            indices.push_back(stor_indices[i]);
            values.push_back(stor_values[i]);
         }
      
         // order indices and values
         out << "sort row" << std::endl;
         generate_sorted_order(indices,order);     
         apply_ordering(indices,order);     
         apply_ordering(values,order);     
      
         out << "grabbing row index, and row norm" << std::endl;
         int rowTerm = basis->index(sa_BasisPerDRow[determDof]->term(stochBasis));

         out << "checking matrix" << std::endl;
         // test middle row values
         int offset = 0;
         for(int stochColBasisIndex = 0;stochColBasisIndex<3;stochColBasisIndex++) {
            for(int stochCol=0;stochCol<adaptMngr.getColStochasticBasisSize(stochColBasisIndex);stochCol++) {
               int colTerm = basis->index(adaptMngr.getColStochasticBasis(stochColBasisIndex)->term(stochCol));
   
               if(big(rowTerm,colTerm)) { 
                  TEST_EQUALITY(indices[offset],adaptMngr.getGlobalColId(stochColBasisIndex,stochCol));
                  TEST_EQUALITY(values[offset],stencil[stochColBasisIndex]*Cijk->getValue(rowTerm,colTerm,0));
                  offset++;
               }
            }
            out << "offset = " << offset << std::endl;
         }
      }

      TEST_ASSERT(sg_A->Filled());
   }
}

// Test construction of Linear Algebra (LA) data structures
TEUCHOS_UNIT_TEST(tBuildAdaptLA, test_uniform)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else 
      Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_SerialComm);
   #endif

   int num_KL = 3;
   int porder = 3;

   Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = buildBasis(num_KL,porder);
   Teuchos::RCP<Epetra_CrsGraph> determGraph = buildTridiagonalGraph(10,*comm);

   std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > sa_BasisPerDRow(determGraph->NumMyRows(),basis);

   // build adapted row map for uniform grid
   Teuchos::RCP<Epetra_Map> rowMap;
   std::vector<int> sa_RowGidOffsets;
   {
      int determMyRows = determGraph->RowMap().NumMyElements();
      int determGlobalRows = determGraph->RowMap().NumGlobalElements();
      
      rowMap = Stokhos::adapt_utils::buildAdaptedRowMapAndOffsets(determGraph->Comm(),sa_BasisPerDRow,sa_RowGidOffsets);
   
      // tests some sizes
      TEST_EQUALITY(rowMap->NumMyElements(),determMyRows*basis->size());
      TEST_EQUALITY(rowMap->NumGlobalElements(),determGlobalRows*basis->size());
      TEST_EQUALITY(int(sa_RowGidOffsets.size()),determMyRows);
      TEST_ASSERT(rowMap->LinearMap());
   
      // test some content
      { 
         bool result = true;
         for(std::size_t i=0;i<sa_RowGidOffsets.size();i++)
            result &= (sa_RowGidOffsets[i]==rowMap->GID(i*basis->size()));
         TEST_ASSERT(result);
      }
   }

   // build adapted column map for uniform grid
   std::vector<int> sa_ColGidOffsets;
   {
      Stokhos::adapt_utils::buildAdaptedColOffsets(*determGraph,
                                                   sa_RowGidOffsets,
                                                   sa_ColGidOffsets);

      const Epetra_BlockMap & d_rowMap = determGraph->RowMap();
      const Epetra_BlockMap & d_colMap = determGraph->ColMap();
      
      int determMyCols = d_colMap.NumMyElements();
      
      TEST_EQUALITY(int(sa_ColGidOffsets.size()),determMyCols);

      bool result = true;
      bool checkOne = false;
      for(int localColId=0;localColId<determMyCols;localColId++) { 
         int localRowId = d_rowMap.LID(d_colMap.GID(localColId));

         if(localRowId==-1) 
            continue; // global comlumn not row on this processor

         checkOne = true; 
         result &= (sa_ColGidOffsets[localColId]==sa_RowGidOffsets[localRowId]); 
      }

      TEST_ASSERT(checkOne);
      TEST_ASSERT(result);
   }
}

TEUCHOS_UNIT_TEST(tBuildAdaptLA, test_graph)
{
   #ifdef HAVE_MPI
      Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
   #else 
      Teuchos::RCP<const Epetra_Comm> comm = Teuchos::rcp(new Epetra_SerialComm);
   #endif

   int numProc = comm->NumProc();
   int rank = comm->MyPID();

   out << "NumProc = " << numProc << ", Rank = " << rank << std::endl;

   int numDetermRows = 3;
   int num_KL = 3;
   int porder = 3;

   Teuchos::RCP<const Stokhos::CompletePolynomialBasis<int,double> > basis = buildBasis(num_KL,porder);
   Teuchos::RCP<Epetra_CrsGraph> determGraph = buildTridiagonalGraph(numDetermRows,*comm);

   std::vector<Teuchos::RCP<const Stokhos::ProductBasis<int,double> > > sa_BasisPerDRow(numDetermRows);
   for(int i=0;i<numDetermRows;i+=3) {
      sa_BasisPerDRow[i] = buildBasis(num_KL,1);
      if(i+1<numDetermRows)
         sa_BasisPerDRow[i+1] = buildBasis(num_KL,2);
      if(i+2<numDetermRows)
         sa_BasisPerDRow[i+2] = buildBasis(num_KL,3);
   }

   for(int i=0;i<numDetermRows;i++) {
      out << "Row " << i << ":\n";
      Stokhos::BasisInteractionGraph(*sa_BasisPerDRow[i]).printGraph(out); 
      out << std::endl;
   } 

   for(int i=1;i<numDetermRows;i++) {
      out << "Pair row " << i-1 << ", col " << i << ":\n";
      Stokhos::BasisInteractionGraph(*basis,*sa_BasisPerDRow[i-1],*sa_BasisPerDRow[i]).printGraph(out); 
      out << std::endl;

      out << "Pair row " << i << ", col " << i-1 << ":\n";
      Stokhos::BasisInteractionGraph(*basis,*sa_BasisPerDRow[i],*sa_BasisPerDRow[i-1]).printGraph(out); 
      out << std::endl;

   } 

   Teuchos::RCP<Epetra_CrsGraph> graph = Stokhos::adapt_utils::buildAdaptedGraph(*determGraph,basis,sa_BasisPerDRow);

   TEST_ASSERT(graph!=Teuchos::null);

   // compute expected number of non zeros
   std::size_t nnz = 0;
   for(int i=0;i<numDetermRows;i++) {
      int gid = determGraph->GRID(i);
      int indices[3];
      int numRowEntries = 0;

      determGraph->ExtractGlobalRowCopy(gid,3,numRowEntries,indices); 

      for(int c=0;c<numRowEntries;c++)
         nnz += Stokhos::BasisInteractionGraph(*basis,*sa_BasisPerDRow[i],*sa_BasisPerDRow[indices[c]]).numNonZeros();
   }
 
   TEST_EQUALITY(graph->NumMyNonzeros(),int(nnz));
}
