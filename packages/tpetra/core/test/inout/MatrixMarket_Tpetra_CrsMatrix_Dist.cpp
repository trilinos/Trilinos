// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This program tests matrix creation and matrix apply using matrices with
// arbitrarily distributed nonzeros (not necessarily row-based distribution).
//
// Create global matrix nonzeros randomly; store all global nonzeros on 
// each proc in a std::map.
// Create distributed vectors with randomized entries using Trilinos' default
// maps 
// For each test (linear row-wise distribution, linear column-wise distribution,
//                random distribution of nonzeros (2D) to processors)
//    distribute matrix nonzeros (each proc selects subset of global nonzeros)
//    create distributed CrsMatrix
//    perform SpMV (nMatvec SpMVs)
//    return result of SpMV
// Compare norms of results of SpMV from all distributions; they should be the
// same.
//
// NOTE:  timings are also reported but should be interpreted carefully.
// This test does not attempt to optimize the distribution of the vectors to
// minimize communication costs.  Moreover, 2D random distribution of nonzeros
// can lead to high communication volume; a better 2D approach would be a 
// block-based approach that better aligns vector entries with matrix entries.

#include "Tpetra_Core.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include "Teuchos_StackedTimer.hpp"
#include "Teuchos_Array.hpp"
#include "MatrixMarket_Tpetra.hpp"

// TestReader:  
// Assume Tpetra's existing matrix-market reader is correct.
//    Read the matrix; do SpMV with resulting matrix; compute norms
//
// Exercise the file reader multiple ways:
//   1D 
//   1D Random
//   1D with 1D partition
//   2D
//   2D Random
//   2D with 1D partition
//   2D with nonzero value = part assignment
// 
// For each, read the matrix; do SpMV; compute norms; compare to baseline

class TestReader {
public:
  using map_t = Tpetra::Map<>;
  using gno_t = map_t::global_ordinal_type;
  using scalar_t = Tpetra::Details::DefaultTypes::scalar_type; 
  using matrix_t = Tpetra::CrsMatrix<scalar_t>;
  using vector_t = Tpetra::Vector<scalar_t>;
  using reader_t = Tpetra::MatrixMarket::Reader<matrix_t>;

  //////////////////////////////
  // Constructor
  TestReader(
    const std::string filename_, 
    const std::string diagonal_,
    const Teuchos::RCP<const Teuchos::Comm<int> > &comm_
  ) : filename(filename_), diagonal(diagonal_), comm(comm_), norm_baseline(3)
  {
    // Compute baseline for testing
    // Use standard MatrixMarket reader from Tpetra to get baseline matrix

    Teuchos::RCP<matrix_t> A_baseline;
    try {
      auto timer = Teuchos::TimeMonitor::getNewTimer("Read:  baseline");
      Teuchos::TimeMonitor tt(*timer);

      A_baseline = reader_t::readSparseFile(filename, comm,
                                            true, false, false);
    }
    catch (std::exception &e) {
      if (comm->getRank() == 0) {
        std::cout << "FAIL:  baseline matrix reading failed "
                  << filename << std::endl;
        std::cout << e.what() << std::endl;
      }
      throw e;
    }

    if (diagonal == "exclude") {
      A_baseline->resumeFill();
      auto rowMap = A_baseline->getRowMap();
      size_t nMyRows = rowMap->getLocalNumElements();
      for (size_t i = 0; i < nMyRows; i++) {
        gno_t gid = rowMap->getGlobalElement(i);
        scalar_t val = Teuchos::ScalarTraits<scalar_t>::zero();
        A_baseline->replaceGlobalValues(gid, 1, &val, &gid);
      }
      A_baseline->fillComplete();
    }

    nRow = A_baseline->getRowMap()->getMaxAllGlobalIndex() 
         + 1;  // Since Trilinos' reader converts one-based to zero-based
    nCol = A_baseline->getColMap()->getMaxAllGlobalIndex() 
         + 1;  // Since Trilinos' reader converts one-based to zero-based
    nNz = A_baseline->getGlobalNumEntries();

    Teuchos::RCP<const map_t> domainMap =
                              Teuchos::rcp(new map_t(nCol, 0, comm));
    x_baseline = Teuchos::rcp(new vector_t(A_baseline->getDomainMap()));
    x_baseline->randomize();

    yin_baseline = Teuchos::rcp(new vector_t(A_baseline->getRangeMap()));
    yin_baseline->randomize();

    // Apply baseline matrix to vectors.
    
    yout_baseline = applyMatrix("baseline", *A_baseline);
    norm_baseline[0] = yout_baseline->normInf();
    norm_baseline[1]=  yout_baseline->norm1();
    norm_baseline[2] = yout_baseline->norm2();
  }


  //////////////////////////////
  // Various combinations to be run for each file
  int run() {

    int ierr = 0;
    int np = comm->getSize();

    // 1D Linear decomposition
    {
      Teuchos::ParameterList params;
      const std::string testname = "1D";
      params.set("diagonal", diagonal);
      params.set("distribution", "1D");
      params.set("randomize", false);
      params.set("useTimers", true);
      ierr += runTest(testname, params);
    }

    // 1D Random
    {
      Teuchos::ParameterList params;
      const std::string testname = "1DRandom";
      params.set("diagonal", diagonal);
      params.set("distribution", "1D");
      params.set("randomize", true);
      params.set("useTimers", true);
      ierr += runTest(testname, params);
    }

    // 2D block decomposition
    {
      Teuchos::ParameterList params;
      const std::string testname = "2D";
      params.set("diagonal", diagonal);
      params.set("distribution", "2D");
      params.set("randomize", false);
      params.set("useTimers", true);
      ierr += runTest(testname, params);
    }

    // 2D Random
    {
      Teuchos::ParameterList params;
      const std::string testname = "2DRandom";
      params.set("diagonal", diagonal);
      params.set("distribution", "2D");
      params.set("randomize", true);
      params.set("useTimers", true);
      ierr += runTest(testname, params);
    }

    // 2D parameter testing:  user-specified npRow
    {
      Teuchos::ParameterList params;
      const std::string testname = "2D_npRow";
      int npRow = 3;
      params.set("diagonal", diagonal);
      params.set("distribution", "2D");
      params.set("randomize", false);
      params.set("nProcessorRows", npRow);
      params.set("useTimers", true);
      try {
        ierr += runTest(testname, params);
      }
      catch (std::exception &e) {
        if (np % npRow) {
          // runTest should fail on when np is not a multiple of npRow; ignore
          // the throw in this case
          if (comm->getRank() == 0) 
            std::cout << "Correctly caught error in " << testname << std::endl;
        }
        else {
          // Test should have passed; this error is real
          throw e;
        }
      }
    }

    // 2D parameter testing:  user-specified npCol
    {
      Teuchos::ParameterList params;
      const std::string testname = "2D_npCol";
      int npCol = 3;
      params.set("diagonal", diagonal);
      params.set("distribution", "2D");
      params.set("randomize", false);
      params.set("nProcessorCols", npCol);
      params.set("useTimers", true);
      try {
        ierr += runTest(testname, params);
      }
      catch (std::exception &e) {
        if (np % npCol) {
          // runTest should fail on when np is not a multiple of npCol; ignore
          // the throw in this case
          if (comm->getRank() == 0) 
            std::cout << "Correctly caught error in " << testname << std::endl;
        }
        else {
          // Test should have passed; this error is real
          throw e;
        }
      }
    }

    // 2D parameter testing:  user-specified npCol and npRow
    {
      Teuchos::ParameterList params;
      const std::string testname = "2D_npCol_npRow";
      int npCol = 3;
      int npRow = 2;
      params.set("diagonal", diagonal);
      params.set("distribution", "2D");
      params.set("randomize", false);
      params.set("nProcessorCols", npCol);
      params.set("nProcessorRows", npRow);
      params.set("useTimers", true);
      try {
        ierr += runTest(testname, params);
      }
      catch (std::exception &e) {
        if (npRow * npCol != np) {
          // runTest should fail on when npRow * npCol != np; ignore
          // the throw in this case
          if (comm->getRank() == 0) 
            std::cout << "Correctly caught error in " << testname << std::endl;
        }
        else {
          // Test should have passed; this error is real
          throw e;
        }
      }
    }

    // LowerTriangularBlock partition
    {
      Teuchos::ParameterList params;
      const std::string testname = "LowerTriangularBlock";
      params.set("diagonal", diagonal);
      params.set("distribution", "LowerTriangularBlock");
      params.set("useTimers", true);
      try {
        ierr += runTestOp(testname, params);
      }
      catch (std::exception &e) {
        int q = int(std::sqrt(float(2 * np)));
        if (q * (q + 1) != 2 * np) {
          // runTest should fail with this processor count; ignore
          // the throw in this case
          if (comm->getRank() == 0) 
            std::cout << "Correctly caught error in " << testname << std::endl;
        }
        else {
          // Test should have passed; this error is real
          throw e;
        }
      }
    }

    // LowerTriangularBlock partition with rows sorted by degree
    {
      Teuchos::ParameterList params;
      const std::string testname = "LowerTriangularBlockSorted";
      params.set("diagonal", diagonal);
      params.set("distribution", "LowerTriangularBlock");
      params.set("sortByDegree", true);
      params.set("useTimers", true);
      try {
        ierr += runTestOp(testname, params);
      }
      catch (std::exception &e) {
        int q = int(std::sqrt(float(2 * np)));
        if (q * (q + 1) != 2 * np) {
          // runTest should fail with this processor count; ignore
          // the throw in this case
          if (comm->getRank() == 0) 
            std::cout << "Correctly caught error in " << testname << std::endl;
        }
        else {
          // Test should have passed; this error is real
          throw e;
        }
      }
    }

    // Todo: add more
    //   1D with 1D partition
    //   2D with 1D partition
    //   2D with nonzero value = part assignment
    return ierr;
  }
  
private:

  //////////////////////////////
  // Read matrix from a MatrixMarket file
  // Distribute the matrix as specified by the parameters
  Teuchos::RCP<matrix_t> readFile(
    const std::string &testname,
    const Teuchos::ParameterList &params,
    Teuchos::RCP<Tpetra::Distribution<gno_t, scalar_t> > &dist
  )
  {
    Teuchos::RCP<matrix_t> Amat;
    try {
      std::string tname = std::string("Read:  ") + testname;
      auto timer = Teuchos::TimeMonitor::getNewTimer(tname);
      Teuchos::TimeMonitor tt(*timer);
      Amat = reader_t::readSparseFile(filename, comm, params, dist);
    }
    catch (std::exception &e) {
      if (comm->getRank() == 0) {
        std::cout << "In test " << testname 
                  << ":  matrix reading failed " << filename 
                  << std::endl;
        std::cout << e.what() << std::endl;
      }
      throw e;
    }
    return Amat;
  }

  //////////////////////////////
  // Apply a matrix to a vector; 
  // compute three norms of product vector
  Teuchos::RCP<vector_t> applyMatrix(
    const std::string testname,
    const Tpetra::Operator<scalar_t> &matrix)
  {
    // Create xvec and yvec using domain and range map of matrix
    vector_t yvec(matrix.getRangeMap());
    vector_t xvec(matrix.getDomainMap());

    // Set xvec to have same values as x_baseline
    Tpetra::Export<> exporterX(x_baseline->getMap(), xvec.getMap());
    xvec.doExport(*x_baseline, exporterX, Tpetra::INSERT);

    // Set yvec to have same values as yin_baseline
    Tpetra::Export<> exporterY(yin_baseline->getMap(), yvec.getMap());
    yvec.doExport(*yin_baseline, exporterY, Tpetra::INSERT);

    // Do matvecs with the matrix; do several so we can time them
    const int nMatvecs = 1;
    std::string tname = std::string("SpMV:  ") + testname;
    auto timer = Teuchos::TimeMonitor::getNewTimer(tname);
    for (int n = 0; n < nMatvecs; n++) {
      Teuchos::TimeMonitor tt(*timer);
      matrix.apply(xvec, yvec, Teuchos::NO_TRANS, 3., 2.);
    }

    // Import result to a vector with default Tpetra map for easy comparisons
    Teuchos::RCP<map_t> defMap = 
             rcp(new map_t(yvec.getGlobalLength(),0,yvec.getMap()->getComm()));
    Teuchos::RCP<vector_t> ydef = 
             rcp(new vector_t(defMap, yvec.getNumVectors()));
    Tpetra::Export<> exportDef(yvec.getMap(), defMap);
    ydef->doExport(yvec, exportDef, Tpetra::INSERT);

    return ydef;
  }
  
  //////////////////////////////
  //  Compare computed norms to the baseline norms
  int compareToBaseline(
    const std::string testname, 
    const Teuchos::RCP<vector_t> &y_test
  )
  {
    const scalar_t epsilon = 10*Teuchos::ScalarTraits<scalar_t>::squareroot(Teuchos::ScalarTraits<scalar_t>::eps());

    int ierr = 0;

    // First compare the norms of the result vector to the baseline
    Teuchos::Array<scalar_t> norm_test(3);
    norm_test[0] = y_test->normInf();
    norm_test[1]=  y_test->norm1();
    norm_test[2] = y_test->norm2();

    for (size_t i = 0; i < static_cast<size_t>(norm_baseline.size()); i++) {
      if (std::abs(norm_baseline[i] - norm_test[i]) > epsilon) {
        ierr++;
        if (comm->getRank() == 0) {
          std::cout << "FAIL in test " << testname << ":  norm " << i << ": "
                    << std::abs(norm_baseline[i] - norm_test[i]) 
                    << " > " << epsilon
                    << std::endl;
          std::cout << "FAIL in test " << testname 
                    << ": baseline " << norm_baseline[i]
                    << ": test " << norm_test[i]
                    << std::endl;
        }
      }
    }

    // If norms match, make sure the vector entries match, too 
    // (Norms are indifferent to errors in permutation)
    if (ierr == 0) {
      //y_test->sync_host();
      //yout_baseline->sync_host();
      auto ytestData = y_test->getLocalViewHost(Tpetra::Access::ReadOnly);
      auto ybaseData = yout_baseline->getLocalViewHost(Tpetra::Access::ReadOnly);
      for (size_t i = 0; i < y_test->getLocalLength(); i++) {
        if (std::abs(ytestData(i,0) - ybaseData(i,0)) > epsilon) ierr++;
      }
      if (ierr > 0) {
        std::cout << "FAIL in test " << testname << ": " 
                  << ierr << " vector entries differ on rank "
                  << comm->getRank() << std::endl;
      }
    }

    int gerr = 0;
    Teuchos::reduceAll<int, int>(*comm, Teuchos::REDUCE_SUM, 1, &ierr, &gerr);
    return gerr;
  }

  //////////////////////////////
  // Each test reads, applies, and compares
  int runTestOp(
    const std::string &testname,
    Teuchos::ParameterList &params
  )
  {
    if (comm->getRank() == 0) 
      std::cout << "\n\nBEGIN " << testname << "\n" << std::endl;

    params.set("verbose", true);
    params.set("callFillComplete", true);

    using dist_t = Tpetra::Distribution<gno_t, scalar_t>;
    Teuchos::RCP<dist_t> dist;

    Teuchos::RCP<matrix_t> Amat = readFile(testname, params, dist);

    using distltb_t = Tpetra::DistributionLowerTriangularBlock<gno_t, scalar_t>;
    Tpetra::LowerTriangularBlockOperator<scalar_t> lto(Amat, 
                               dynamic_cast<distltb_t&>(*dist));

    Teuchos::RCP<vector_t> yvec = applyMatrix(testname, lto);

    return compareToBaseline(testname, yvec);
  }

  //////////////////////////////
  // Each test reads, applies, and compares
  int runTest(
    const std::string &testname,
    Teuchos::ParameterList &params
  )
  {
    if (comm->getRank() == 0) 
      std::cout << "\n\nBEGIN " << testname << "\n" << std::endl;

    params.set("verbose", true);
    params.set("callFillComplete", true);

    Teuchos::RCP<Tpetra::Distribution<gno_t, scalar_t> > dist;  // Not used
    Teuchos::RCP<matrix_t> Amat = readFile(testname, params, dist);

    Teuchos::RCP<vector_t> yvec = applyMatrix(testname, *Amat);

    return compareToBaseline(testname, yvec);
  }
  
private:

  const std::string filename;    // MatrixMarket filename
  const std::string diagonal;    // Special handling of the matrix diagonal:  
                                 // require or exclude?
  const Teuchos::RCP<const Teuchos::Comm<int> > comm;

  size_t nRow, nCol, nNz;                  // dimensions of baseline matrix
  Teuchos::RCP<vector_t> x_baseline;       // SpMV input vector y = by + aAx
  Teuchos::RCP<vector_t> yin_baseline;     // SpMV input vector y = by + aAx
  Teuchos::RCP<vector_t> yout_baseline;    // SpMV output vector y = by + aAx
  Teuchos::Array<scalar_t> norm_baseline;  // Norm of y for baseline matrix
};

////////////////////////////////////////////////////////////////////////////

int main(int narg, char *arg[]) 
{
  Tpetra::ScopeGuard scope(&narg, &arg);
  const Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  Teuchos::RCP<Teuchos::StackedTimer> stackedTimer =
                        rcp(new Teuchos::StackedTimer("main"));
  Teuchos::TimeMonitor::setStackedTimer(stackedTimer);

  int ierr = 0;

  // Get the filename (set up tests with general, pattern, symmetric, etc.)
  std::string filename = "";
  std::string diagonal = "";

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("file", &filename,
                 "Path and filename of the matrix to be read.");
  cmdp.setOption("diagonal", &diagonal,
                 "Options are exclude or require.");
  if (cmdp.parse(narg,arg)!=Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  TestReader test(filename, diagonal, comm);
  ierr += test.run();

  // Output timing info
  stackedTimer->stop("main");

  Teuchos::StackedTimer::OutputOptions options;
  options.output_fraction = true;
  options.output_histogram = true;
  options.output_minmax = true;

  Teuchos::FancyOStream out(Teuchos::rcpFromRef(std::cout));
  out.setOutputToRootOnly(0);

  stackedTimer->report(out, comm, options);
  stackedTimer = Teuchos::null;

  // Report test result
  if (comm->getRank() == 0) {
    if (ierr == 0) 
      std::cout << "End Result: TEST PASSED" << std::endl;
    else  
      std::cout << "End Result: TEST FAILED with " << ierr << " failures" 
                << std::endl;
  }

  return ierr;
}

