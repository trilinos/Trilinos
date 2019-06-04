#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
#include <stdexcept>

#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_Time.hpp>

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_LinearProblem.h"
#if defined(EPETRA_MPI)
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif // defined(EPETRA_MPI)

#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"

#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector_def.hpp>
#include <Tpetra_Vector_def.hpp>
#include <Tpetra_RowMatrix_def.hpp>
#include <Tpetra_Experimental_BlockMultiVector.hpp>
#include <Tpetra_Experimental_BlockCrsMatrix.hpp>

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

namespace
{

using std::vector;
using std::string;
using Teuchos::RCP;
using Teuchos::rcp;
using tpetra_crs_matrix_type = Tpetra::CrsMatrix<double>;
using tpetra_map_type = Tpetra::Map<>;

/**************************************************************
 * Create a Tpetra 2D laplacian matrix for testing relaxation *
 **************************************************************/

RCP<const tpetra_crs_matrix_type> create_testmat_t(const RCP<const tpetra_map_type>& rowmap)
{
  using GO = tpetra_map_type::global_ordinal_type;

  //tile this over each diagonal element:
  // 0 -1  0
  //-1  4 -1
  // 0 -1  0
  RCP<tpetra_crs_matrix_type> crsmatrix = rcp(new tpetra_crs_matrix_type(rowmap, 3));
  GO numRows = rowmap->getGlobalNumElements();
  double tile[3] = {-2, 4, -2};
  GO firstCols[] = {0, 1};
  crsmatrix->insertGlobalValues(0, 2, tile + 1, firstCols);
  for(GO row = 1; row < numRows - GO(1); row++)
  {
    GO cols[] = {row - 1, row, row + 1};
    crsmatrix->insertGlobalValues(row, 3, tile, cols);
  }
  GO lastCols[] = {numRows - GO(2), numRows - GO(1)};
  crsmatrix->insertGlobalValues(numRows - GO(1), 2, tile, lastCols);
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
  constexpr bool debug = true;
  auto myOutPtr = debug ?
    Teuchos::getFancyOStream (Teuchos::rcpFromRef (std::cerr)) :
    Teuchos::rcpFromRef (out);
  auto& myOut = *myOutPtr;

  myOut << "Ifpack2 BlockRelaxation performance comparison: Epetra vs. Tpetra"
        << std::endl;

  using Teuchos::ArrayRCP;
  //use the same types as epetra for fair comparison

  decltype (Tpetra::getDefaultComm()) tcomm;
  try {
    tcomm = Tpetra::getDefaultComm();
  }
  catch (std::exception& e) {
    myOut << "Tpetra::getDefaultComm() threw an exception: " << e.what () << std::endl;
    TEST_ASSERT( false );
    return;
  }
  catch (...) {
    myOut << "Tpetra::getDefaultComm() threw an exception "
      "not a subclass of std::exception" << std::endl;
    TEST_ASSERT( false );
    return;
  }

  if(tcomm->getSize() == 1)
  {
    //Configure ranges of testing
    const vector<int> blockSizes = {2, 4, 5, 8, 10, 20, 50};
    string container = "TriDi";
    using GO = tpetra_map_type::global_ordinal_type;
    const GO localRows = 5000;

    myOut << "Create Tpetra test matrix" << std::endl;

    auto tmap = rcp(new tpetra_map_type(localRows, 0, tcomm));
    auto tmat = create_testmat_t(tmap);

    myOut << "Create Epetra test matrix" << std::endl;
#if defined(EPETRA_MPI)
    Epetra_MpiComm ecomm(MPI_COMM_WORLD);
#else
    Epetra_SerialComm ecomm;
#endif
    Epetra_Map emap(static_cast<int> (localRows), 0, ecomm);
    Teuchos::RCP<const Epetra_CrsMatrix> emat;
    try {
      // mfh 26 Jul 2018: This function call throws.  That means the
      // test is broken; I didn't write that code.  This blocks #3139,
      // so I'm going to just let this test build for now, but not
      // actually run it.
      emat = create_testmat_e(emap);
    }
    catch (std::exception& e) {
      myOut << "create_testmat_e threw an exception: " << e.what () << std::endl;
      TEST_ASSERT( false );
      return;
    }
    catch (...) {
      myOut << "create_testmat_e threw an exception "
        "not a subclass of std::exception" << std::endl;
      TEST_ASSERT( false );
      return;
    }

    myOut << "Create and fill Tpetra Vector(s)" << std::endl;

    //set up matrix
    typedef Tpetra::RowMatrix<double> ROW;
    typedef Tpetra::Vector<double> TVec;
    TVec lhs(tmap);
    TVec rhs(tmap);
    rhs.putScalar(2.0);

    myOut << "Open output file" << std::endl;

    //set up output file
    FILE* csv = fopen("BlockRelax.csv", "w");
    if (csv == nullptr) {
      TEST_ASSERT( false );
      myOut << "Failed to open output file \"BlockRelax.csv\"!" << std::endl;
      return;
    }
    fputs("lib,blockSize,setup,apply\n", csv);

    myOut << "Benchmark Ifpack2" << std::endl;
    //////////////////////
    //-- Ifpack2 test --//
    //////////////////////
    Teuchos::ParameterList ilist;
    ilist.set("partitioner: type", "user");
    ilist.set("relaxation: sweeps", 5);
    ilist.set("relaxation: type", "Gauss-Seidel");
    int* partMap = new int[localRows];
    Teuchos::Time timer("");
    myOut << "Testing Ifpack2 block G-S\n";
    //target for total running time
    const double maxTotalTime = 6.0;
    //proportion of total time to spend on ifpack2
    //adjust depending on relative speeds to balance trial count
    const double ifpack2TimeProportion = 0.9;
    const double ifpack2TimePerSet = (ifpack2TimeProportion * maxTotalTime) / blockSizes.size();
    const double mlTimePerSet = ((1.0 - ifpack2TimeProportion) * maxTotalTime) / blockSizes.size();
    for(auto blockSize : blockSizes)
    {
      myOut << "    Testing block size: " << blockSize << '\n';
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
      myOut << '(' << trials << " trials)\n";
      double applyTime = timer.totalElapsedTime() / trials;
      fprintf(csv, "Ifpack2,%i,%f\n", blockSize, applyTime);
    }
    delete[] partMap;

    myOut << "Benchmark ML" << std::endl;
    /////////////////
    //-- ML test --//
    /////////////////
    myOut << "Testing ML block gs\n";
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
        myOut << "    Testing block size: " << blockSize << '\n';
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
        myOut << '(' << trials << " trials)\n";
        double applyTime = timer.totalElapsedTime() / trials;
        fprintf(csv,"ML,%i,%f\n", blockSize, applyTime);
      }
    }
    myOut << "Done with performance testing: data written to BlockRelax.csv"
          << std::endl;
    fclose(csv);
  }
}

} //anonymous namespace
