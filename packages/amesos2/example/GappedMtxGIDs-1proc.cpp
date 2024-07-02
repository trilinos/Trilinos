// @HEADER
// *****************************************************************************
//           Amesos2: Templated Direct Sparse Solver Package
//
// Copyright 2011 NTESS and the Amesos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/**
   \file   SimpleSolveNonContigMap.cpp
   \author Nathan Ellingwood <ndellin@sandia.gov>
   \author Mark Hoemmen <mhoemme@sandia.gov>
   \date   Wed Feb 27 15:12:39 2019

   \brief  Example of Amesos2 usage reading a matrix and  map with non-contiguous GIDs


   This example uses a custom matrix market file reader added by Mark Hoemmen to
   read a matrix and map with gapped or non-contiguous GIDs

   This example then solves a sparse system of linear equations using the
   Amesos2 interface

   The example is hard-coded for 1 MPI processes
*/

#include <iostream>
#include <cstdio>
#include <iostream>

#include <Teuchos_RCP.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_MultiVector.hpp>
#include <Tpetra_CrsMatrix.hpp>

#include <MatrixMarket_Tpetra.hpp>
#include <Tpetra_Import.hpp>

#include "Amesos2.hpp"
#include "Amesos2_Version.hpp"


using map_type = Tpetra::Map<>;
using LO = map_type::local_ordinal_type;
using GO = map_type::global_ordinal_type;
using Scalar = double;
using MAT = Tpetra::CrsMatrix<Scalar>;
using MV = Tpetra::MultiVector<Scalar>;
using reader_type = Tpetra::MatrixMarket::Reader<MAT>;


// Fwd decl - codes from Mark Hoemmen
bool
readEntryFromFile (GO& gblRowInd, GO& gblColInd, Scalar& val, const std::string& s);


// NOTE: This is not a general routine, despite occasional efforts at
// being so.  It expects a particular input file and only works with a
// particular test, namely the test for Issue #3647.
Teuchos::RCP<MAT>
readCrsMatrixFromFile (const std::string& matrixFilename,
                       Teuchos::RCP<Teuchos::FancyOStream> & fos,
                       const Teuchos::RCP<const map_type>& rowMap,
                       const Teuchos::RCP<const map_type>& domainMap,
                       const Teuchos::RCP<const map_type>& rangeMap,
                       const bool convert_to_zero_base = true,
                       const int header_size=2);


// ------------------------------------------------

int main(int argc, char *argv[]) {
  using Tpetra::global_size_t;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::tuple;
  using std::endl;

  typedef double Scalar;
  typedef Tpetra::Map<>::local_ordinal_type LO;
  typedef Tpetra::Map<>::global_ordinal_type GO;

  typedef Tpetra::CrsMatrix<Scalar,LO,GO> MAT;
  typedef Tpetra::MultiVector<Scalar,LO,GO> MV;

  Tpetra::ScopeGuard tpetraScope(&argc,&argv);
  {

    auto comm = Tpetra::getDefaultComm();
    const int myRank = comm->getRank();

    Teuchos::oblackholestream blackhole;

    // command-line processing
    bool allprint = true;
    bool verbose = false;

#if 1
    std::string mtx_name = "gap-ids-1procA.mm";
    std::string map_name = "gap-ids-1proc-rowmap.mm";
    std::string rhs_name = "gap-ids-1procRhs.mm";
#else
    std::string mtx_name = "badMatrix.mm";
    std::string map_name = "rowMap.mm";
    std::string rhs_name = "badRhs.mm";
#endif

    std::string solverName = "KLU2";

    int num_header_lines = 2;
    bool convert_mtx_to_zero_base = true;
    bool make_contiguous = true;

    Teuchos::CommandLineProcessor cmdp(false,true);
    cmdp.setOption("allPrint","rankZeroPrint",&allprint,"Print output for all ranks or rank zero.");
    cmdp.setOption("verbose","quiet",&verbose,"Print messages and results.");
    cmdp.setOption("mapFilename", &map_name, "Name of Map file");
    cmdp.setOption("matrixFilename", &mtx_name, "Name of sparse matrix "
                   "file in Matrix Market format, to load as the linear system "
                   "matrix A");
    cmdp.setOption("rhsFilename", &rhs_name, "Name of dense (multi)vector "
                   "file in Matrix Market format, to load as the linear system "
                   "right-hand side(s) B");
    cmdp.setOption("solverName", &solverName, "Amesos2 solver name");
    cmdp.setOption("numMtxHeaderLines", &num_header_lines, "Header lines to remove from matrix market file");
    cmdp.setOption("convertMtxZeroBase","keepMtxOneBase",&convert_mtx_to_zero_base, "Set this option to leave matrix with index base 1");
    cmdp.setOption("makeContiguous","isContiguous",&make_contiguous, "Set this option to makeContiguous if matrix has gapped row ids");

    if (cmdp.parse(argc,argv) != Teuchos::CommandLineProcessor::PARSE_SUCCESSFUL) {
      return -1;
    }

    std::ostream& out = ( (allprint || (myRank == 0)) ? std::cout : blackhole );
    RCP<Teuchos::FancyOStream> fos = Teuchos::fancyOStream(Teuchos::rcpFromRef(out));

    // Say hello
    *fos << myRank << " : " << " Hello - Running GappedMtxGIDs test" << endl;

    // Step 1: Read in data
    const global_size_t numProcs = comm->getSize();
    *fos << "Running on numProcs = " << numProcs << endl;

    RCP< const map_type > rowMap = Tpetra::MatrixMarket::Reader< MAT >::readMapFile(map_name, comm);

    // Scope-guarded check of rowMap
    {
      std::ostringstream outMap;
      Tpetra::MatrixMarket::Writer<MAT>::writeMap (outMap, *rowMap);
      std::istringstream inMap (outMap.str ());
      auto rowMap2 = reader_type::readMap (inMap, comm);

      const bool same = rowMap->isSameAs (*rowMap2);
      TEUCHOS_TEST_FOR_EXCEPTION
        (! same, std::logic_error, "map -> readMapFile -> writeMap -> readMap "
         "does not result in the same Map.");
    }

    if ( myRank == 0 && verbose ) {
      *fos << "\nrowMap->describe output:" << endl;
      rowMap->describe(*fos, Teuchos::VERB_EXTREME);
    }

    RCP<MAT> A;
    RCP<const map_type> domainMap;
    RCP<const map_type> rangeMap;
    domainMap = rowMap;
    rangeMap = rowMap;

    if (mtx_name.find ("badMatrix") != std::string::npos) {
      *fos << "badMatrix.mm input matrix - special case" << endl;
      num_header_lines = 4;
      convert_mtx_to_zero_base=false;
      *fos << "  num_header_lines to remove: " << num_header_lines << endl;
      *fos << "  convert mtx base hard-coded to false in this case. checking value: " << convert_mtx_to_zero_base << endl;
      A = readCrsMatrixFromFile (mtx_name, fos, rowMap, domainMap, rangeMap, convert_mtx_to_zero_base, num_header_lines);
    }
    else {
      *fos << "Preparing to read in matrix file with gapped rowids" << endl;
      *fos << "  num_header_lines to remove: " << num_header_lines << endl;
      *fos << "  convert mtx base: " << convert_mtx_to_zero_base << endl;
      A = readCrsMatrixFromFile (mtx_name, fos, rowMap, domainMap, rangeMap, convert_mtx_to_zero_base, num_header_lines);
    }

    if ( myRank == 0 && verbose ) {
      *fos << "A->describe" << endl;
      A->describe(*fos, Teuchos::VERB_EXTREME);
    }


    RCP<MV> RHS;
    RHS = Tpetra::MatrixMarket::Reader<MAT>::readDenseFile (rhs_name, comm, rangeMap);
    if ( myRank == 0 && verbose ) {
      *fos << "RHS->describe" << endl;
      RHS->describe(*fos, Teuchos::VERB_EXTREME);
    }


    // Amesos2 compute solution
    const int numVectors = 1;
    RCP<MV> Xhat = rcp(new MV(rowMap,numVectors));


    if( !Amesos2::query(solverName) ){
      std::cerr << solverName << " not enabled.  Exiting..." << endl;
      return EXIT_SUCCESS;        // Otherwise CTest will pick it up as
                                  // failure, which it isn't really
    }

    RCP<Amesos2::Solver<MAT,MV> > solver = Amesos2::create<MAT,MV>(solverName, A, Xhat, RHS);

    if( myRank == 0 ){
      *fos << "\nAmesos2: Solver description:  " << solver->description() << endl;
    }

    // Create a Teuchos::ParameterList to hold solver parameters
    if( myRank == 0 ){
      *fos << "Amesos2: set Params" << endl;
    }
    Teuchos::ParameterList amesos2_params("Amesos2");
    if ( make_contiguous ) {
      if( myRank == 0 ) { *fos << "  set IsContigous==false in solver parameter list" << endl; }
      amesos2_params.sublist(solverName).set("IsContiguous", false, "Are GIDs Contiguous");
    }
    solver->setParameters( Teuchos::rcpFromRef(amesos2_params) );

    if( myRank == 0 ){
      *fos << "Amesos2: call symbolic" << endl;
    }
    solver->symbolicFactorization();

    if( myRank == 0 ){
      *fos << "Amesos2: call numeric" << endl;
    }
    solver->numericFactorization();

    if( myRank == 0 ){
      *fos << "Amesos2: call solve" << endl;
    }
    solver->solve();


    // Check result
    {
      // Compute: RHShat = A*Xhat
      RCP<MV> RHShat;
      RHShat = rcp(new MV(rowMap,numVectors));
      A->apply(*Xhat, *RHShat);

      // norm2(RHShat - RHS)
      typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
      Teuchos::Array<Magnitude> xhatnorms(numVectors);
      RHS->update(-1.0, *RHShat, 1.0);
      RHS->norm2(xhatnorms());

      if( myRank == 0 ){
        *fos << "\nsolve finished" << endl;
        *fos << "\nnorm2(Ax - b) = " << xhatnorms << endl;
      }
    }

    // Store modified RHS and re-solve, reusing the symbolic and numeric results
    RCP<MV> newRHS;
    newRHS = Tpetra::MatrixMarket::Reader<MAT>::readDenseFile (rhs_name, comm, rangeMap);
    newRHS->scale(2.0);
    if( myRank == 0 ){
      *fos << "\n\nRetest: created new RHS to test" << endl;
    }

    solver->setB( newRHS );
    if( myRank == 0 ){
      *fos << "Amesos2: call solve with newRHS" << endl;
    }
    solver->solve();

    // Check result
    {
      // Compute: RHShat = A*Xhat
      RCP<MV> RHShat;
      RHShat = rcp(new MV(rowMap,numVectors));
      A->apply(*Xhat, *RHShat);

      // norm2(RHShat - newRHS)
      typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
      Teuchos::Array<Magnitude> xhatnorms(numVectors);
      newRHS->update(-1.0, *RHShat, 1.0);
      newRHS->norm2(xhatnorms());

      if( myRank == 0 ){
        *fos << "\nsecond solve finished" << endl;
        *fos << "\nnorm2(Ax - b) = " << xhatnorms << endl;
      }
    }


  } // end ScopeGuard

  return 0;
} // end main



bool
readEntryFromFile (GO& gblRowInd, GO& gblColInd, Scalar& val, const std::string& s)
{
  if (s.size () == 0 || s.find ("%") != std::string::npos) {
    return false; // empty line or comment line
  }
  std::istringstream in (s);

  if (! in) {
    return false;
  }
  in >> gblRowInd;
  if (! in) {
    return false;
  }
  in >> gblColInd;
  if (! in) {
    return false;
  }
  in >> val;
  return true;
}

// NOTE: This is not a general routine, despite occasional efforts at
// being so.  It expects a particular input file and only works with a
// particular test, namely the test for Issue #3647.
Teuchos::RCP<MAT>
readCrsMatrixFromFile (const std::string& matrixFilename,
                       Teuchos::RCP<Teuchos::FancyOStream> & fos,
                       const Teuchos::RCP<const map_type>& rowMap,
                       const Teuchos::RCP<const map_type>& domainMap,
                       const Teuchos::RCP<const map_type>& rangeMap,
                       const bool convert_to_zero_base,
                       const int header_size)
{
  auto comm = rowMap->getComm ();
  const int myRank = comm->getRank ();

  std::ifstream inFile;
  int opened = 0;
  if (myRank == 0) 
  {
    try {
      inFile.open (matrixFilename);
      if (inFile) {
        opened = 1;
      }
    }
    catch (...) {
      opened = 0;
    }
  }
  Teuchos::broadcast<int, int> (*comm, 0, Teuchos::outArg (opened));
  TEUCHOS_TEST_FOR_EXCEPTION
    (opened == 0, std::runtime_error, "readCrsMatrixFromFile: "
     "Failed to open file \"" << matrixFilename << "\" on Process 0.");

  using Teuchos::RCP;
  RCP<MAT> A;

  if (myRank == 0) 
  {
    std::string line;

    // Skip the first N lines.  This is a hack, specific to the input file in question.

    //*fos << "  Reading matrix market file. Skip " << header_size << "  header lines" << std::endl;
    for ( int i = 0; i < header_size; ++i ) {
      std::getline (inFile, line);
    } 

    std::map<GO, size_t> counts;
    Teuchos::Array<Scalar> vals;
    Teuchos::Array<GO> gblRowInds;
    Teuchos::Array<GO> gblColInds;
    while (inFile) {
      std::getline (inFile, line);
      GO gblRowInd {};
      GO gblColInd {};
      Scalar val {};
      const bool gotLine = readEntryFromFile (gblRowInd, gblColInd, val, line);
      if (gotLine) {
        //*fos << "  read mtx rank: " << myRank << " |  gblRowInd = " << gblRowInd << "  gblColInd = " << gblColInd << std::endl;
        if ( convert_to_zero_base ) {
          gblRowInd -= 1 ;
          gblColInd -= 1 ;
        }
        counts[gblRowInd]++;
        vals.push_back(val);
        gblRowInds.push_back(gblRowInd);
        gblColInds.push_back(gblColInd);
      }
    }

    // Max number of entries in any row
    using pair_type = decltype(counts)::value_type;
    auto pr = std::max_element(
        std::begin(counts),
        std::end(counts),
        [] (pair_type const& p1, pair_type const& p2){ return p1.second < p2.second; }
    );
    size_t maxCount = (counts.empty()) ? size_t(0) : pr->second;
    A = Teuchos::rcp(new MAT(rowMap, maxCount));
    for (typename Teuchos::Array<GO>::size_type i=0; i<gblRowInds.size(); i++) {
      A->insertGlobalValues (gblRowInds[i], gblColInds(i,1), vals(i,1));
    }
  }

  A->fillComplete (domainMap, rangeMap);
  return A;
}

