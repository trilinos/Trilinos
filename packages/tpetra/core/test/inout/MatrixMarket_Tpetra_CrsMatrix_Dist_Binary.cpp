// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

// This program tests binary readers with the following distributions: 1D,
// 1DRandom, 2D, 2DRandom, 2D_npRow, 2D_npCol, 2D_npCol_npRow, LowerTriangularBlock. 
// This program is very similar to the one in MatrixMarket_Tpetra_CrsMatrix_Dist.cpp
// The constructor reads the matrix from a matrix market (mtx) file using the 
// standard MatrixMarket reader and creates A_baseline. With random input vector 
//  x_baseline, SpMV is performed on A_baseline and yin_baseline, and the output 
// vector is stored in yout_baseline for future comparisons. If global binary reader 
// is being tested, the constructor also makes rank 0 write a global binary file,
// in which the order of the nonzeros is the same as the mtx file. The destructor 
// at rank 0 removes the global binary file.

// Testing the global reader:
// For each tested distribution, we read the global binary file and create a 
// matrix called AmatTest with the desired distribution. We perform SpMV on AmatTest
// with x_baseline and yin_baseline and compare the output vector with yout_baseline.  

// Testing the per-process reader:
// For each tested distribution, we read the input matrix (mtx) file again 
// and create a matrix called AmatWriter with the desired distribution. We use
// AmatWriter to write the nonzeros in a global binary file or multiple per-process
// binary files. Then we read the binary file(s) and create AmatTest with the same
// distribution as AmatWriter. We perfom SpMV on AmatTest with x_baseline and
// yin_baseline and compare the output vector with yout_baseline. 

// IMPORTANT:
//    Since the current binary formats do not support any numeric values,
//    the binary readers set all numeric values to one. Therefore, we also set 
//    all numeric values to one in this program.


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
  using indices_type = typename matrix_t::nonconst_global_inds_host_view_type;

  //////////////////////////////
  // Constructor
  TestReader(
    const std::string filename_, 
    const bool perProcess_,
    const Teuchos::RCP<const Teuchos::Comm<int> > &comm_
  ) : filename(filename_), perProcess(perProcess_), comm(comm_), norm_baseline(3)
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
    A_baseline->setAllToScalar(ONE);

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

    // Write the global binary file if perProcess==false 
    if(!perProcess)
      writeBinaryGlobal();
  }

  // Destructor deletes the global binary file if it exists
  ~TestReader() {
    if(!perProcess && comm->getRank() == 0) {
      try { std::remove(binfilename.c_str()); }
      catch (std::exception &e) {
    	std::cout << "Could not delete file: " << binfilename << std::endl;
    	std::cout << e.what() << std::endl;
      }	      
    }
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
      params.set("distribution", "1D");
      params.set("randomize", false);
      ierr += runTest(testname, params);
    }

    // 1D Random
    {
      Teuchos::ParameterList params;
      const std::string testname = "1DRandom";
      params.set("distribution", "1D");
      params.set("randomize", true);
      ierr += runTest(testname, params);
    }

    // 2D block decomposition
    {
      Teuchos::ParameterList params;
      const std::string testname = "2D";
      params.set("distribution", "2D");
      params.set("randomize", false);
      ierr += runTest(testname, params);
    }

    // 2D Random
    {
      Teuchos::ParameterList params;
      const std::string testname = "2DRandom";
      params.set("distribution", "2D");
      params.set("randomize", true);
      ierr += runTest(testname, params);
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
      params.set("distribution", "2D");
      params.set("randomize", false);
      params.set("nProcessorCols", npCol);
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
      params.set("distribution", "2D");
      params.set("randomize", false);
      params.set("nProcessorCols", npCol);
      params.set("nProcessorRows", npRow);
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
      params.set("distribution", "LowerTriangularBlock");
      try {
    	ierr += runTestLTB(testname, params);
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
    if(!perProcess)
    {
      Teuchos::ParameterList params;
      const std::string testname = "LowerTriangularBlockSorted";
      params.set("distribution", "LowerTriangularBlock");
      params.set("sortByDegree", true);
      try {
        ierr += runTestLTB(testname, params);
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
  // Write a single global binary file, which will be used for all tests (1D, 2D, etc) 
  void writeBinaryGlobal()
  {
    binfilename = filename + "_" + std::to_string(comm->getSize()) + ".cooBin"; 
   
    int err = 0;

    // Only rank 0 will read the mtx file and write the corresponding binary file
    // Other ranks will wait for rank 0 to finish
    if(comm->getRank() == 0) {

      unsigned int nRows, nCols, entry[2];
      unsigned long long nNzs;
      std::string line;
      
      // Open the input file and skip the header
      // NOTE: symmetric keyword is ignored too.
      std::ifstream in(filename, std::ios::in);
      do
	std::getline (in, line);
      while(line[0] == '%');
      
      std::stringstream sstream(line);
      sstream >> nRows >> nCols >> nNzs;
      if(nRows != nCols)
	err = 1;

      // Proceed with the binary file only if the input mtx is square
      if(err == 0) {
	
	// Open the output file and write the header
	std::ofstream out(binfilename, std::ios::out | std::ios::binary);
	out.write((char *)&nRows, sizeof(unsigned int));
	out.write((char *)&nNzs, sizeof(unsigned long long));
	
	// Write nonzeros
	while(std::getline(in, line)) {
	  std::stringstream sstream2(line);
	  sstream2 >> entry[0] >> entry[1];
	  out.write((char *)entry, sizeof(unsigned int)*2);	
	}
	out.close();
      }

      in.close();
    }

    Teuchos::broadcast(*comm, 0, 1, &err);

    if(err == 1) {
      throw std::runtime_error( "Input matrix " + filename +  " is not square.");
    }
      
  }

  //////////////////////////////
  // Write per-process binary files: each rank writes its own nonzeros to a separate file 
  // The path for the written files should be unique for each test
  void writeBinaryPerProcess(
    const std::string &testname,
    const Teuchos::RCP<matrix_t> &AmatWrite
  )
  {
    // Open the file
    binfilename = filename + "_" + std::to_string(comm->getSize()) + "_" + testname; 
    std::string binrankfilename = binfilename + "." + std::to_string(comm->getRank()) + ".cooBin";
    std::ofstream out(binrankfilename, std::ios::out | std::ios::binary);

    // Write the header
    unsigned int nRows = static_cast<unsigned int>(AmatWrite->getRowMap()->getMaxAllGlobalIndex()) + 1;
    unsigned long long  nNzs = static_cast<unsigned long long>(AmatWrite->getLocalNumEntries());
    out.write((char *)& nRows, sizeof(unsigned int));
    out.write((char *)& nRows, sizeof(unsigned int));
    out.write((char *)& nNzs, sizeof(unsigned long long));

    // Get the CrsGraph because we do not need the values
    auto graph = AmatWrite->getCrsGraph();	
    auto rowMap = graph->getRowMap();
    indices_type gblColInds;
    size_t numEntries = 0;

    // Write the nonzeros
    unsigned int entry[2];
    for(size_t r = 0; r < graph->getLocalNumRows(); r++) {

      // Get the global index for row r
      auto gblRow = rowMap->getGlobalElement(static_cast<gno_t>(r));
      entry[0] = static_cast<unsigned int>(gblRow) + 1;
      
      // Get the copy of the row with global column indices
      numEntries = graph->getNumEntriesInGlobalRow(gblRow);
      Kokkos::resize(gblColInds,numEntries);
      graph->getGlobalRowCopy(gblRow, gblColInds, numEntries);
      
      // Write the entries in the row in COO format (i.e., in "rowId colId" pairs)
      for(size_t c = 0; c < numEntries; c++) {
	entry[1] = static_cast<unsigned int>(gblColInds[c]) + 1;
	out.write((char *)entry, sizeof(unsigned int)*2);
      }
    }

    out.close();
    
  }

  //////////////////////////////
  // Remove the binary file(s) with the given path
  void cleanBinaryPerProcess(
  )
  {
    // There exists a binary file for each process, so each process removes its own file.  
    std::string binrankfilename = binfilename + "." + std::to_string(comm->getRank()) + ".cooBin";
    try { std::remove(binrankfilename.c_str()); }
    catch (std::exception &e) {
      std::cout << "Could not delete file: " << binrankfilename << std::endl;
      std::cout << e.what() << std::endl;
      throw e;
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
  // Run the per-process reader test or the global reader test 
  int runTest(
    const std::string &testname,
    Teuchos::ParameterList &params
  )
  {
    params.set("verbose", true);
    params.set("callFillComplete", true);
    params.set("useTimers", true);

    if(perProcess)
      return runTestPerProcess(testname, params);
    else
      return runTestGlobal(testname, params);
  }

  //////////////////////////////
  // Each test reads the binary file, applies, and compares
  int runTestGlobal(
    const std::string &testname,
    Teuchos::ParameterList &params
  )
  {
    if (comm->getRank() == 0) 
      std::cout << "\n\nBEGIN GLOBAL TEST " << testname << "\n" << std::endl;

    // Read the binary file
    params.set("binary", true);  // This invokes the binary reader
    Teuchos::RCP<Tpetra::Distribution<gno_t, scalar_t> > distTest;  // Not used
    Teuchos::RCP<matrix_t> AmatTest = readFile(binfilename, testname, params, distTest);

    // Apply and compare
    Teuchos::RCP<vector_t> yvec = applyMatrix(testname, *AmatTest);
    return compareToBaseline(testname, yvec);
  }

  //////////////////////////////
  // Each test reads the mtx file, creates a distributed matrix AmatWriter,
  //  writes the per-process binary files using AmatWriter, reads the per-process
  //  files, applies, and compares
  int runTestPerProcess(
    const std::string &testname,
    Teuchos::ParameterList &params
  )
  {
    if (comm->getRank() == 0) 
      std::cout << "\n\nBEGIN PER-PROCESS TEST " << testname << "\n" << std::endl;

    // Create AmatWriter by reading the input mtx file with the required distribution
    Teuchos::RCP<Tpetra::Distribution<gno_t, scalar_t> > distWriter;  // Not used
    Teuchos::RCP<matrix_t> AmatWriter = readFile(filename, testname, params, distWriter);

    // Write the binary files using AmatWriter
    writeBinaryPerProcess(testname, AmatWriter);

    // Read the binary files
    params.set("binary", true);  // This invokes binary readers.
    params.set("readPerProcess", true);
    Teuchos::RCP<Tpetra::Distribution<gno_t, scalar_t> > distTest;  // Not used
    Teuchos::RCP<matrix_t> AmatTest = readFile(binfilename, testname, params, distTest);

    // Clean-up the per-process binary files
    cleanBinaryPerProcess();

    // Apply and compare
    Teuchos::RCP<vector_t> yvec = applyMatrix(testname, *AmatTest);
    return compareToBaseline(testname, yvec);
  }

  //////////////////////////////
  // Each test reads, writes binary file(s), reads the binary file(s), applies, and compares
  // This test should not compare the current results with baseline's results,
  // because the LTB distribution only keeps the lower triangular entries.
  // Therefore, this test compares the current results with AmatWriter's results.
  int runTestLTB(
    const std::string &testname,
    Teuchos::ParameterList &params
  )
  {
    params.set("verbose", true);
    params.set("callFillComplete", true);
    params.set("useTimers", true);

    if(perProcess)
      return runTestLTBPerProcess(testname, params);
    else
      return runTestLTBGlobal(testname, params);
  }

  //////////////////////////////
  // Each test reads the binary file, applies, and compares
  int runTestLTBGlobal(
    const std::string &testname,
    Teuchos::ParameterList &params
  )
  {
    if (comm->getRank() == 0) 
      std::cout << "\n\nBEGIN GLOBAL TEST " << testname << "\n" << std::endl;

    // Read the binary file
    params.set("binary", true);   // This invokes the binary reader 
    Teuchos::RCP<Tpetra::Distribution<gno_t, scalar_t> > distTest;  // Not used
    Teuchos::RCP<matrix_t> AmatTest = readFile(binfilename, testname, params, distTest);

    // Get the LTB operator from the test matrix
    using distltb_t = Tpetra::DistributionLowerTriangularBlock<gno_t, scalar_t>;
    Tpetra::LowerTriangularBlockOperator<scalar_t> lto_test(AmatTest, 
                               dynamic_cast<distltb_t&>(*distTest));

    // Apply and compare
    Teuchos::RCP<vector_t> yvec = applyMatrix(testname, lto_test);
    return compareToBaseline(testname, yvec);
  }

  //////////////////////////////
  // Each test reads the mtx file, creates a distributed matrix AmatWriter,
  //  writes the per-process binary files using AmatWriter, reads the per-process
  //  files, applies, and compares
  int runTestLTBPerProcess(
    const std::string &testname,
    Teuchos::ParameterList &params
  )
  {
    if (comm->getRank() == 0) 
      std::cout << "\n\nBEGIN PER-PROCESS TEST " << testname << "\n" << std::endl;

    // Create AmatWriter by reading the input mtx file with the required distribution
    Teuchos::RCP<Tpetra::Distribution<gno_t, scalar_t> > distWriter;  // Not used
    Teuchos::RCP<matrix_t> AmatWriter = readFile(filename, testname, params, distWriter);

    // Write the binary files using AmatWriter
    writeBinaryPerProcess(testname, AmatWriter);

    // Read the binary files
    params.set("binary", true);  // This invokes binary readers.
    params.set("readPerProcess", true);
    Teuchos::RCP<Tpetra::Distribution<gno_t, scalar_t> > distTest;  // Not used
    Teuchos::RCP<matrix_t> AmatTest = readFile(binfilename, testname, params, distTest);

    // Clean-up the per-process binary files
    cleanBinaryPerProcess();

    // Get the LTB operator from the test matrix
    using distltb_t = Tpetra::DistributionLowerTriangularBlock<gno_t, scalar_t>;
    Tpetra::LowerTriangularBlockOperator<scalar_t> lto_test(AmatTest, 
                               dynamic_cast<distltb_t&>(*distTest));

    // Apply and compare
    Teuchos::RCP<vector_t> yvec = applyMatrix(testname, lto_test);
    return compareToBaseline(testname, yvec);
  }

  
private:

  const std::string filename;    // MatrixMarket filename
  const bool perProcess;         // Test the per-process reader instead of global
  std::string binfilename;       // Binary (global) filename
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

  TestReader test(filename, perProcess, comm);
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

