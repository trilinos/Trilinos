#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <stdexcept>

#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Time.hpp>

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#include "Epetra_SerialComm.h"

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector_def.hpp>
#include <Tpetra_Vector_def.hpp>
#include <Tpetra_RowMatrix_def.hpp>
#include <Tpetra_Experimental_BlockMultiVector.hpp>
#include <Tpetra_Experimental_BlockCrsMatrix.hpp>

#include <Ifpack2_ConfigDefs.hpp>
#include <Ifpack2_Version.hpp>
#include <Ifpack2_UnitTestHelpers.hpp>
#include <Ifpack2_BlockRelaxation.hpp>
#include <Ifpack2_SparseContainer.hpp>
#include <Ifpack2_TriDiContainer.hpp>
#include <Ifpack2_BandedContainer.hpp>
#include <Ifpack2_DenseContainer.hpp>
#include <Ifpack2_OverlappingPartitioner.hpp>
#include <Ifpack2_LinearPartitioner.hpp>
#include <Ifpack2_LinePartitioner.hpp>
#include <Ifpack2_ILUT.hpp>
#include <Ifpack2_Factory.hpp>
#include "Ifpack2_ETIHelperMacros.h"

namespace
{

using std::vector;
using std::string;
using Teuchos::RCP;
using Teuchos::rcp;
typedef unsigned long long u64;
typedef Tpetra::Map<>::node_type NO;
typedef Tpetra::CrsMatrix<double,int,int,NO> TMat;
typedef Tpetra::Map<int,int,NO> TMap;

/**************************************************************
 * Create a Tpetra 2D laplacian matrix for testing relaxation *
 **************************************************************/

RCP<const TMat> create_testmat_t(const RCP<const TMap>& rowmap)
{
  //tile this over each diagonal element:
  // 0 -1  0
  //-1  4 -1
  // 0 -1  0
  RCP<TMat> crsmatrix = rcp(new TMat(rowmap, 3));
  int numRows = rowmap->getGlobalNumElements();
  double tile[3] = {-2, 4, -2};
  int firstCols[] = {0, 1};
  crsmatrix->insertGlobalValues(0, 2, tile + 1, firstCols);
  for(int row = 1; row < numRows - 1; row++)
  {
    int cols[] = {row - 1, row, row + 1};
    crsmatrix->insertGlobalValues(row, 3, tile, cols);
  }
  int lastCols[] = {numRows - 2, numRows - 1};
  crsmatrix->insertGlobalValues(numRows - 1, 2, tile, lastCols);
  crsmatrix->fillComplete();
  return crsmatrix;
}

/***************************************************************
 * Create an Epetra 2D laplacian matrix for testing relaxation *
 ***************************************************************/

RCP<const Epetra_CrsMatrix> create_testmat_e(const Epetra_Map& rowmap)
{
  RCP<Epetra_CrsMatrix> crsmatrix = rcp(new Epetra_CrsMatrix(Epetra_DataAccess::Copy, rowmap, 3));
  int numRows = rowmap.NumGlobalElements();
  double tile[3] = {-2, 4, -2};
  int firstCols[] = {0, 1};
  crsmatrix->InsertGlobalValues(0, 2, tile + 1, firstCols);
  for(int row = 1; row < numRows - 1; row++)
  {
    int cols[] = {row - 1, row, row + 1};
    crsmatrix->InsertGlobalValues(row, 3, tile, cols);
  }
  int lastCols[] = {numRows - 2, numRows - 1};
  crsmatrix->InsertGlobalValues(numRows - 1, 2, tile, lastCols);
  crsmatrix->FillComplete();
  return crsmatrix;
}

/*****************************************
 * Test New Ifpack2 vs ML *
 *****************************************/

TEUCHOS_UNIT_TEST(Ifpack2BlockRelaxation, Performance)
{
  using Teuchos::ArrayRCP;
  //use the same types as epetra for fair comparison
  auto tcomm = Tpetra::getDefaultComm();
  if(tcomm->getSize() == 1)
  {
    //Configure ranges of testing
    const vector<int> blockSizes = {2, 4, 5, 8, 10, 20, 50};
    string container = "TriDi";
    const int localRows = 5000;
    auto tmap = rcp(new TMap(localRows, 0, tcomm));
    auto tmat = create_testmat_t(tmap);
#ifdef EPETRA_MPI
    Epetra_MpiComm ecomm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm ecomm;
#endif
    Epetra_Map emap(localRows, 0, ecomm);
    auto emat = create_testmat_e(emap);
    //set up matrix
    typedef Tpetra::RowMatrix<double,int,int,NO> ROW;
    typedef Tpetra::Vector<double,int,int,NO> TVec;
    TVec lhs(tmap);
    TVec rhs(tmap);
    rhs.putScalar(2.0);
    //set up output file
    FILE* csv = fopen("BlockRelax.csv", "w");
    fputs("lib,blockSize,setup,apply\n", csv);
    //////////////////////
    //-- Ifpack2 test --//
    //////////////////////
    Teuchos::ParameterList ilist;
    ilist.set("partitioner: type", "user");
    ilist.set("relaxation: sweeps", 5);
    ilist.set("relaxation: type", "Gauss-Seidel");
    int* partMap = new int[localRows];
    Teuchos::Time timer("");
    out << "Testing Ifpack2 block G-S\n";
    //target for total running time
    const double maxTotalTime = 6.0;
    //proportion of total time to spend on ifpack2
    //adjust depending on relative speeds to balance trial count
    const double ifpack2TimeProportion = 0.9;
    const double ifpack2TimePerSet = (ifpack2TimeProportion * maxTotalTime) / blockSizes.size();
    const double mlTimePerSet = ((1.0 - ifpack2TimeProportion) * maxTotalTime) / blockSizes.size();
    for(auto blockSize : blockSizes)
    {
      out << "    Testing block size: " << blockSize << '\n';
      for(int i = 0; i < localRows; i++)
        partMap[i] = i / blockSize;
      ArrayRCP<int> PMrcp(partMap, 0, localRows, false);
      ilist.set("partitioner: map", PMrcp);
      ilist.set("partitioner: local parts", partMap[localRows - 1] + 1);
      ilist.set("relaxation: container", container);
      Ifpack2::BlockRelaxation<ROW> relax(tmat);
      relax.setParameters(ilist);
      relax.initialize();
      relax.compute();
      timer.reset();
      timer.start();
      int trials = 0;
      do
      {
        relax.apply(rhs, lhs);
        trials++;
      }
      while(timer.totalElapsedTime(true) < ifpack2TimePerSet);
      timer.stop();
      out << '(' << trials << " trials)\n";
      double applyTime = timer.totalElapsedTime() / trials;
      fprintf(csv, "Ifpack2,%i,%f\n", blockSize, applyTime);
    }
    delete[] partMap;
    /////////////////
    //-- ML test --//
    /////////////////
    out << "Testing ML block gs\n";
    {
      using ML_Epetra::MultiLevelPreconditioner;
      Epetra_Vector X(emap);
      Epetra_Vector Y(emap);
      double* yval;
      Y.ExtractView(&yval);
      for(int i = 0; i < localRows; i++)
        yval[i] = 2;
      double* xval;
      X.ExtractView(&xval);
      Teuchos::ParameterList mllist;
      mllist.set("max levels", 1);
      mllist.set("smoother: type", "block Gauss-Seidel");
      mllist.set("smoother: sweeps", 5);
      mllist.set("smoother: pre or post", "pre");
      mllist.set("smoother: damping factor", 1.0);
      mllist.set("ML output", 0);
      for(auto blockSize : blockSizes)
      {
        //manually set dofs per node
        mllist.set("PDE equations", blockSize);
        MultiLevelPreconditioner prec(*emat, mllist, true);
        prec.ComputePreconditioner();
        auto ml = prec.GetML();
        auto sm = ml->pre_smoother;
        out << "    Testing block size: " << blockSize << '\n';
        timer.reset();
        timer.start();
        int trials = 0;
        do
        {
          ML_Smoother_BlockGS(sm, localRows, xval, localRows, yval);
          trials++;
        }
        while(timer.totalElapsedTime(true) < mlTimePerSet);
        timer.stop();
        out << '(' << trials << " trials)\n";
        double applyTime = timer.totalElapsedTime() / trials;
        fprintf(csv,"ML,%i,%f\n", blockSize, applyTime);
      }
    }
    out << "Done with performance testing: data written to BlockRelax.csv\n";
    fclose(csv);
  }
}

IFPACK2_ETI_MANGLING_TYPEDEFS()

} //anonymous namespace
