/*
// @HEADER
// ***********************************************************************
//
//          Tpetra: Templated Linear Algebra Services Package
//                 Copyright (2008) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// ************************************************************************
// @HEADER
*/

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
    const Teuchos::RCP<const Teuchos::Comm<int> > &comm_
  ) : filename(filename_), comm(comm_), norm_baseline(3)
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
    
    // Set all values in A to one: This is needed because binary readers 
    // do not currently support numeric values
    const scalar_t ONE = Teuchos::ScalarTraits<scalar_t>::one();
    A_baseline->resumeFill();
    A_baseline->setAllToScalar(ONE);
    A_baseline->fillComplete(A_baseline->getDomainMap(), A_baseline->getRangeMap());

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
  int run(const bool perProcess) {

    int ierr = 0;
    int np = comm->getSize();

    // 1D Linear decomposition
    {
      Teuchos::ParameterList params;
      const std::string testname = "1D";
      params.set("distribution", "1D");
      params.set("randomize", false);
      ierr += runTest(testname, params, perProcess);
    }

    // 1D Random
    {
      Teuchos::ParameterList params;
      const std::string testname = "1DRandom";
      params.set("distribution", "1D");
      params.set("randomize", true);
      ierr += runTest(testname, params, perProcess);
    }

    // 2D block decomposition
    {
      Teuchos::ParameterList params;
      const std::string testname = "2D";
      params.set("distribution", "2D");
      params.set("randomize", false);
      ierr += runTest(testname, params, perProcess);
    }

    // 2D Random
    {
      Teuchos::ParameterList params;
      const std::string testname = "2DRandom";
      params.set("distribution", "2D");
      params.set("randomize", true);
      ierr += runTest(testname, params, perProcess);
    }

    // 2D parameter testing:  user-specified npRow
    {
      Teuchos::ParameterList params;
      const std::string testname = "2D_npRow";
      int npRow = 3;
      params.set("distribution", "2D");
      params.set("randomize", false);
      params.set("nProcessorRows", npRow);
      try {
    	ierr += runTest(testname, params, perProcess);
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
      params.set("distribution", "2D");
      params.set("randomize", false);
      params.set("nProcessorCols", npCol);
      try {
    	ierr += runTest(testname, params, perProcess);
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
      params.set("distribution", "2D");
      params.set("randomize", false);
      params.set("nProcessorCols", npCol);
      params.set("nProcessorRows", npRow);
      try {
    	ierr += runTest(testname, params, perProcess);
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
      params.set("distribution", "LowerTriangularBlock");
      try {
    	ierr += runTestLTB(testname, params, perProcess);
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
  // Write one or more binary files depending on the value of perProcess
  void writeBinary(
    const std::string &testname,
    const Teuchos::RCP<matrix_t> &AmatWrite,
    const bool perProcess,
    std::string &binfilename
  )
  {
    // if(perProcess)
    //   writeBinaryPerProcess(testname, AmatWrite);
    // else 
    writeBinaryGlobal(testname, AmatWrite, binfilename);
  }

  //////////////////////////////
  // Write a single binary file 
  // The path for the written files should be unique for each test
  void writeBinaryGlobal(
    const std::string &testname,
    const Teuchos::RCP<matrix_t> &AmatWrite,
    std::string &binfilename
  )
  {
    binfilename = filename + "_" + std::to_string(comm->getSize())
                           + "_" + testname + ".cooBin"; 
    std::ofstream out;

    // Processes will take turns to write the nonzeros they own
    for(int i = 0; i < comm->getSize(); i++){
      if(i == comm->getRank()) {

	// Open the file and write the header if rank 0
	if(i == 0) {
	  out.open(binfilename, std::ios::out | std::ios::binary);
	  unsigned int nRows = static_cast<unsigned int>(AmatWrite->getRowMap()->getMaxAllGlobalIndex()) + 1;
	  unsigned long long  nNzs = static_cast<unsigned long long>(AmatWrite->getGlobalNumEntries());
	  out.write((char *)& nRows, sizeof(unsigned int));
	  out.write((char *)& nNzs, sizeof(unsigned long long));
	}
	else 
	  out.open(binfilename, std::ios::out | std::ios::binary | std::ios::app);
	
	// Get the CrsGraph because we do not need the values
	auto graph = AmatWrite->getCrsGraph();	
	auto rowMap = graph->getRowMap();
	Teuchos::Array<gno_t> gblColInds;
	size_t numEntries = 0;
	
	unsigned int entry[2];
	for(size_t r = 0; r < graph->getNodeNumRows(); r++) {

	  // Get the global index for row r
	  auto gblRow = rowMap->getGlobalElement(static_cast<gno_t>(r));
	  entry[0] = static_cast<unsigned int>(gblRow) + 1;
	  
	  // Get the copy of the row with global column indices
	  numEntries = graph->getNumEntriesInGlobalRow(gblRow);
	  gblColInds.resize(numEntries);
	  graph->getGlobalRowCopy(gblRow, gblColInds(), numEntries);

	  // Write the entries in the row in COO format (i.e., in "rowId colId" pairs)
	  for(size_t c = 0; c < numEntries; c++) {
	    entry[1] = static_cast<unsigned int>(gblColInds[c]) + 1;
	    out.write((char *)entry, sizeof(unsigned int)*2);
	  }
	}
	out.close();
      }
      comm->barrier();
    }
      
  }

  void cleanBinary(
    const std::string &binfilename,
    const bool perProcess
  )
  {
    if(perProcess) {
    }
    else {
      // make sure none of the ranks is currently reading the file
      comm->barrier();
      if(comm->getRank() == 0) {
	try { std::remove(binfilename.c_str()); }
	catch (std::exception &e) {
	  std::cout << "Could not delete file: " << binfilename << std::endl;
	  std::cout << e.what() << std::endl;
	  throw e;
	}	
      }
    }
      
  }

  //////////////////////////////
  // Read matrix from a MatrixMarket file
  // Distribute the matrix as specified by the parameters
  Teuchos::RCP<matrix_t> readFile(
    const std::string &filetoread,
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
      Amat = reader_t::readSparseFile(filetoread, comm, params, dist);
    }
    catch (std::exception &e) {
      if (comm->getRank() == 0) {
        std::cout << "In test " << testname 
                  << ":  matrix reading failed " << filetoread
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
    const scalar_t epsilon = 0.0000001;
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
      y_test->sync_host();
      yout_baseline->sync_host();
      auto ytestData = y_test->getLocalViewHost();
      auto ybaseData = yout_baseline->getLocalViewHost();
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
  //  Compare computed norms of two given vectors 
  int compareTwoVectors(
    const std::string testname, 
    const Teuchos::RCP<vector_t> &y1,
    const Teuchos::RCP<vector_t> &y2
  )
  {
    const scalar_t epsilon = 0.0000001;
    int ierr = 0;

    // First compare the norms of the vectors
    Teuchos::Array<scalar_t> norm1(3);
    norm1[0] = y1->normInf();
    norm1[1]=  y1->norm1(); 
    norm1[2] = y1->norm2();
    Teuchos::Array<scalar_t> norm2(3);
    norm2[0] = y2->normInf();
    norm2[1]=  y2->norm1(); 
    norm2[2] = y2->norm2();

    for (size_t i = 0; i < static_cast<size_t>(norm1.size()); i++) {
      if (std::abs(norm1[i] - norm2[i]) > epsilon) {
        ierr++;
        if (comm->getRank() == 0) {
          std::cout << "FAIL in test " << testname << ":  norm " << i << ": "
                    << std::abs(norm1[i] - norm2[i]) 
                    << " > " << epsilon
                    << std::endl;
          std::cout << "FAIL in test " << testname 
                    << ": vec1 " << norm1[i]
                    << ": vec2 " << norm2[i]
                    << std::endl;
        }
      }
    }

    // If norms match, make sure the vector entries match, too 
    // (Norms are indifferent to errors in permutation)
    if (ierr == 0) {
      y1->sync_host();
      y2->sync_host();
      auto y1Data = y1->getLocalViewHost();
      auto y2Data = y2->getLocalViewHost();
      for (size_t i = 0; i < y1->getLocalLength(); i++) {
        if (std::abs(y1Data(i,0) - y2Data(i,0)) > epsilon) ierr++;
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
  // This test should not compare the current results with baseline's results,
  // because the LTB distribution only keeps the lower triangular entries.
  // Therefore, this test compares the current results with AmatWriter's results.
  int runTestLTB(
    const std::string &testname,
    Teuchos::ParameterList &params,
    const bool perProcess
  )
  {
    if (comm->getRank() == 0) 
      std::cout << "\n\nBEGIN " << testname << "\n" << std::endl;

    params.set("verbose", true);
    params.set("callFillComplete", true);
    params.set("useTimers", true);

    // Create AmatWriter by reading the input mtx file with the required distribution
    using dist_t = Tpetra::Distribution<gno_t, scalar_t>;
    Teuchos::RCP<dist_t> distWriter;
    Teuchos::RCP<matrix_t> AmatWriter = readFile(filename, testname, params, distWriter);

    // Set the values in Amat_writer to one
    const scalar_t ONE = Teuchos::ScalarTraits<scalar_t>::one();
    AmatWriter->resumeFill();
    AmatWriter->setAllToScalar(ONE);
    AmatWriter->fillComplete(AmatWriter->getDomainMap(), AmatWriter->getRangeMap());

    // Get the LTB operator from both matrices
    using distltb_t = Tpetra::DistributionLowerTriangularBlock<gno_t, scalar_t>;
    Tpetra::LowerTriangularBlockOperator<scalar_t> lto_writer(AmatWriter, 
                               dynamic_cast<distltb_t&>(*distWriter));

    // Write the binary file(s) using AmatWriter 
    std::string binfilename;
    writeBinary(testname, AmatWriter, perProcess, binfilename);

    // Read the binary file(s)
    params.set("binary", true);
    Teuchos::RCP<dist_t> distTest;
    Teuchos::RCP<matrix_t> AmatTest = readFile(binfilename, testname, params, distTest);
   
    // Clean-up the binary file(s)
    cleanBinary(binfilename, perProcess);

    Tpetra::LowerTriangularBlockOperator<scalar_t> lto_test(AmatTest, 
                               dynamic_cast<distltb_t&>(*distTest));

    // Apply with AmatWriter
    Teuchos::RCP<vector_t> yout_writer = applyMatrix("AmatWriter", lto_writer);

    // Apply with AmatTest
    Teuchos::RCP<vector_t> yout_test = applyMatrix(testname, lto_test);

    // Compare
    return compareTwoVectors(testname, yout_writer, yout_test);
  }

  //////////////////////////////
  // Each test reads, applies, and compares
  int runTest(
    const std::string &testname,
    Teuchos::ParameterList &params,
    const bool perProcess
  )
  {
    if (comm->getRank() == 0) 
      std::cout << "\n\nBEGIN " << testname << "\n" << std::endl;

    params.set("verbose", true);
    params.set("callFillComplete", true);
    params.set("useTimers", true);

    // Create AmatWriter by reading the input mtx file with the required distribution
    Teuchos::RCP<Tpetra::Distribution<gno_t, scalar_t> > distWriter;  // Not used
    Teuchos::RCP<matrix_t> AmatWriter = readFile(filename, testname, params, distWriter);

    // Write the binary file(s) using AmatWriter
    std::string binfilename;
    writeBinary(testname, AmatWriter, perProcess, binfilename);

    // Read the binary file(s)
    params.set("binary", true);
    Teuchos::RCP<Tpetra::Distribution<gno_t, scalar_t> > distTest;  // Not used
    Teuchos::RCP<matrix_t> AmatTest = readFile(binfilename, testname, params, distTest);

    // Clean-up the binary file(s)
    cleanBinary(binfilename, perProcess);

    // Apply and compare
    Teuchos::RCP<vector_t> yvec = applyMatrix(testname, *AmatTest);
    return compareToBaseline(testname, yvec);
  }
  
private:

  const std::string filename;    // MatrixMarket filename
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

  // Binary formats do not currently support numeric values, so we will set 
  // all entries to one after reading the matrix. This means you can use any 
  // type (pattern, general, etc) of mtx file in this test. 
  std::string filename = "";
  bool perProcess = false;

  Teuchos::CommandLineProcessor cmdp(false,true);
  cmdp.setOption("file", &filename,
                 "Path and filename of the matrix (.mtx) to be read.");
  cmdp.setOption("per-process", "global", &perProcess,
                 "Whether the per-process or global reader should be tested.");
  if (cmdp.parse(narg,arg)!=Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
    return -1;
  }

  TestReader test(filename, comm);
  ierr += test.run(perProcess);

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

