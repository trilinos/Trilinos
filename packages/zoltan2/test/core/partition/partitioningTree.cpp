// @HEADER
// *****************************************************************************
//   Zoltan2: A package of combinatorial algorithms for scientific computing
//
// Copyright 2012 NTESS and the Zoltan2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Zoltan2_PartitioningProblem.hpp>
#include <Zoltan2_XpetraCrsMatrixAdapter.hpp>
#include <Zoltan2_XpetraMultiVectorAdapter.hpp>
#include <Zoltan2_TestHelpers.hpp>
#include <iostream>
#include <limits>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_FancyOStream.hpp>
#include <Teuchos_CommandLineProcessor.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_Vector.hpp>
#include <MatrixMarket_Tpetra.hpp>

using Teuchos::RCP;

typedef zlno_t z2TestLO;
typedef zgno_t z2TestGO;
typedef zscalar_t z2TestScalar;

typedef Tpetra::CrsMatrix<z2TestScalar, z2TestLO, z2TestGO> SparseMatrix_t;
typedef Tpetra::Vector<z2TestScalar, z2TestLO, z2TestGO> Vector;
typedef Vector::node_type Node;

typedef Tpetra::MultiVector<z2TestScalar, z2TestLO, z2TestGO,znode_t> tMVector_t;


typedef Zoltan2::XpetraCrsMatrixAdapter<SparseMatrix_t,tMVector_t> SparseMatrixAdapter_t;

typedef SparseMatrixAdapter_t::part_t part_t;

typedef Zoltan2::XpetraMultiVectorAdapter<tMVector_t> MultiVectorAdapter_t;



int testForRCB(SparseMatrixAdapter_t &matAdapter, int myrank, part_t numparts,
  RCP<tMVector_t> coords, RCP<const Teuchos::Comm<int> > comm);
int testForPHG(SparseMatrixAdapter_t &matAdapter, int myrank, part_t numparts,
  RCP<tMVector_t> coords, RCP<const Teuchos::Comm<int> >);
int testForMJ(SparseMatrixAdapter_t &matAdapter, int myrank, part_t numparts,
  RCP<tMVector_t> coords, RCP<const Teuchos::Comm<int> >);


////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
int main(int narg, char** arg)
{
  Tpetra::ScopeGuard tscope(&narg, &arg);
  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();
  int me = comm->getRank();

  //////////////////////////////////////////////////////////////////////

  std::string inputFile = "";        // Matrix Market or Zoltan file to read
  std::string inputPath = testDataFilePath;  // Directory with input file
  bool distributeInput = true;
  int success = 0;
  part_t numParts = 8;


  //////////////////////////////////////////////////////////////////////
  // Read run-time options.
  //////////////////////////////////////////////////////////////////////
  Teuchos::CommandLineProcessor cmdp (false, false);
  cmdp.setOption("inputPath", &inputPath,
                 "Path to the MatrixMarket or Zoltan file to be read; "
                 "if not specified, a default path will be used.");
  cmdp.setOption("inputFile", &inputFile,
                 "Name of the Matrix Market or Zoltan file to read; "
                 "");
  cmdp.setOption("distribute", "no-distribute", &distributeInput,
                "for Zoltan input files only, "
                "indicate whether or not to distribute "
                "input across the communicator");
  cmdp.setOption("numParts", &numParts,
                 "Global number of parts;");

  Teuchos::CommandLineProcessor::EParseCommandLineReturn
    parseReturn= cmdp.parse( narg, arg );

  if( parseReturn == Teuchos::CommandLineProcessor::PARSE_HELP_PRINTED )
  {
    return 0;
  }
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Construct matrix from file
  //////////////////////////////////////////////////////////////////////
  RCP<UserInputForTests> uinput;

  if (inputFile != "")   // Input file specified; read a matrix
  {
    uinput = rcp(new UserInputForTests(inputPath, inputFile, comm,
                                       true, distributeInput));
  }
  else
  {
    std::cout << "Input file must be specified." << std::endl;
  }

  RCP<SparseMatrix_t> origMatrix = uinput->getUITpetraCrsMatrix();


  if (me == 0)
  {
    std::cout << "NumRows     = " << origMatrix->getGlobalNumRows() << std::endl
         << "NumNonzeros = " << origMatrix->getGlobalNumEntries() << std::endl
         << "NumProcs = " << comm->getSize() << std::endl
         << "NumParts = " << numParts << std::endl;
  }

  if (origMatrix->getGlobalNumRows() < 40)
  {
    Teuchos::FancyOStream out(Teuchos::rcp(&std::cout,false));
    origMatrix->describe(out, Teuchos::VERB_EXTREME);
  }
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  // Get coordinates, corresponding to graph vertices
  //////////////////////////////////////////////////////////////////////
  RCP<tMVector_t> coords;
  try
  {
    coords = uinput->getUICoordinates();
  }
  catch(...)
  {
    if (me == 0)
      std::cout << "FAIL: get coordinates" << std::endl;
    return 1;
  }
  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  ////// Create an input adapter for the Tpetra matrix.
  //////////////////////////////////////////////////////////////////////
  SparseMatrixAdapter_t matAdapter(origMatrix);

  MultiVectorAdapter_t *ca = NULL;

  try
  {
    ca = new MultiVectorAdapter_t(coords);
  }
  catch(...)
  {
    if (me == 0)
      std::cout << "FAIL: vector adapter" << std::endl;
    return 1;
  }

  matAdapter.setCoordinateInput(ca);
  //////////////////////////////////////////////////////////////////////

  // MDM - MJ disabled - not implemented yet
  // int errMJ  = testForMJ(matAdapter, me, numParts, coords, comm); -check err
  int errRCB = testForRCB(matAdapter, me, numParts, coords, comm);
  int errPHG = testForPHG(matAdapter, me, numParts, coords, comm);

  // Currently just for internal development - sweep from 1 - 25 numParts
  // and validate each to check for any bad results.
  // #define BRUTE_FORCE_SEARCH
  #ifdef BRUTE_FORCE_SEARCH
  for(int setParts = 1; setParts <= 25; ++setParts) {
    // numParts is a parameter so should not be set like for for production
    // code - this is just a 2nd development check to run try a bunch of tests
    numParts = setParts;
    std::cout << "Brute force testing num parts: " << numParts << std::endl;
    // MJ not implemented yet
    // if(testForMJ(matAdapter, me, numParts, coords, comm) != 0) { return 1; }
    if(testForRCB(matAdapter, me, numParts, coords, comm) != 0) { return 1; }
    if(testForPHG(matAdapter, me, numParts, coords, comm) != 0) { return 1; }
  }
  #endif

  delete ca;

  // if(errMJ==0)
  if(errRCB==0)
  if(errPHG==0)
  {
    std::cout << "PASS" << std::endl;
    return success;
  }
  return 1;
}
////////////////////////////////////////////////////////////////////////////////


////////////////////////////////////////////////////////////////////////////////
// Runs validation on the results to make sure the arrays are sensible
////////////////////////////////////////////////////////////////////////////////
bool validate(part_t numTreeVerts,
  std::vector<part_t> permPartNums,
  std::vector<part_t> splitRangeBeg,
  std::vector<part_t> splitRangeEnd,
  std::vector<part_t> treeVertParents
)
{
  // by design numTreeVerts does not include the root
  if(numTreeVerts != static_cast<part_t>(splitRangeBeg.size()) - 1) {
    return false;
  }
  if(numTreeVerts != static_cast<part_t>(splitRangeEnd.size()) - 1) {
    return false;
  }
  if(numTreeVerts != static_cast<part_t>(treeVertParents.size()) - 1) {
    return false;
  }
  // now search every node and validate it's properties
  for(part_t n = 0; n <= numTreeVerts; ++n) {
    if(n < static_cast<part_t>(permPartNums.size())) {
      // terminal so children should just be range 1
      if(splitRangeEnd[n] != splitRangeBeg[n] + 1) {
        std::cout << "Invalid terminal - range should be 1" << std::endl;
        return false;
      }
      if(splitRangeBeg[n] != n) {
        std::cout << "Invalid terminal - not pointing to myself!" << std::endl;
        return false;
      }
    }
    else {
      part_t beg = splitRangeBeg[n];
      part_t end = splitRangeEnd[n];
      part_t count = end - beg;
      std::vector<bool> findChildren(count, false);
      for(part_t n2 = 0; n2 <= numTreeVerts; ++n2) {
        if(treeVertParents[n2] == n) {
          part_t beg2 = splitRangeBeg[n2];
          part_t end2 = splitRangeEnd[n2];
          for(part_t q = beg2; q < end2; ++q) {
            part_t setIndex = q - beg; // rel to parent so beg, not beg2
            if(setIndex < 0 || setIndex >= count) {
              std::cout << "Found bad child index - invalid results!" << std::endl;
              return false;
            }
            if(findChildren[setIndex]) {
              std::cout << "Found child twice - invalid results!" << std::endl;
              return false;
            }
            findChildren[setIndex] = true;
          }
        }
      }
      for(part_t q = 0; q < count; ++q) {
        if(!findChildren[q]) {
          std::cout << "Did not find all children. Invalid results!" << std::endl;
          return false;
        }
      }
    }
    if(n == numTreeVerts) {
      // this is the root
      if(splitRangeBeg[n] != 0) {
        std::cout << "Root must start at 0!" << std::endl;
        return false;
      }
      if(splitRangeEnd[n] != static_cast<part_t>(permPartNums.size())) {
        std::cout << "Root must contain all parts!" << std::endl;
        return false;
      }
    }
  }
  return true;
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// General test function to analyze solution and print out some information
////////////////////////////////////////////////////////////////////////////////
int analyze(Zoltan2::PartitioningProblem<SparseMatrixAdapter_t> &problem,
  RCP<const Teuchos::Comm<int> > comm)
{
  part_t numTreeVerts = 0;
  std::vector<part_t> permPartNums; // peritab in Scotch
  std::vector<part_t> splitRangeBeg;
  std::vector<part_t> splitRangeEnd;
  std::vector<part_t> treeVertParents;

  Zoltan2::PartitioningSolution<SparseMatrixAdapter_t> solution =
    problem.getSolution();

  solution.getPartitionTree(numTreeVerts,permPartNums,splitRangeBeg,
                            splitRangeEnd,treeVertParents);

  comm->barrier(); // for tidy output...
  if(comm->getRank() == 0) { // for now just plot for rank 0

    // print the acquired information about the tree

    // Header
    std::cout << std::endl << "Printing partition tree info..." << std::endl << std::endl;

    // numTreeVerts
    std::cout << "  numTreeVerts: " << numTreeVerts << std::endl << std::endl;

    // array index values 0 1 2 3 ...
    std::cout << "  part array index:";
    for(part_t n = 0; n < static_cast<part_t>(permPartNums.size()); ++n) {
      std::cout << std::setw(4) << n << " ";
    }
    std::cout << std::endl;

    // permParNums
    std::cout << "  permPartNums:    ";
    for(part_t n = 0; n < static_cast<part_t>(permPartNums.size()); ++n) {
      std::cout << std::setw(4) << permPartNums[n] << " ";
    }
    std::cout << std::endl << std::endl;

    // node index values 0 1 2 3 ...
    std::cout << "  node index:      ";
    for(part_t n = 0; n < static_cast<part_t>(splitRangeBeg.size()); ++n) {
      std::cout << std::setw(4) << n << " ";
    }
    std::cout << std::endl;

    // splitRangeBeg
    std::cout << "  splitRangeBeg:   ";
    for(part_t n = 0; n < static_cast<part_t>(splitRangeBeg.size()); ++n) {
      std::cout << std::setw(4) << splitRangeBeg[n] << " ";
    }
    std::cout << std::endl;

    // splitRangeEnd
    std::cout << "  splitRangeEnd:   ";
    for(part_t n = 0; n < static_cast<part_t>(splitRangeEnd.size()); ++n) {
      std::cout << std::setw(4) << splitRangeEnd[n] << " ";
    }
    std::cout << std::endl;

    // treeVertParents
    std::cout << "  treeVertParents: ";
    for(part_t n = 0; n < static_cast<part_t>(treeVertParents.size()); ++n) {
      std::cout << std::setw(4) << treeVertParents[n] << " ";
    }
    std::cout << std::endl << std::endl;
  }
  comm->barrier(); // for tidy output...

  if(!validate(numTreeVerts, permPartNums, splitRangeBeg, splitRangeEnd,
    treeVertParents)) {
    return 1;
  }

  return 0;
}

////////////////////////////////////////////////////////////////////////////////
// Test partitioning tree access for RCB partitioning.  This partitioning tree
// should be binary..
////////////////////////////////////////////////////////////////////////////////
int testForRCB(SparseMatrixAdapter_t &matAdapter, int me, part_t numParts,
  RCP<tMVector_t> coords, RCP<const Teuchos::Comm<int> > comm)
{

  //////////////////////////////////////////////////////////////////////
  ////// Specify problem parameters
  //////////////////////////////////////////////////////////////////////
  Teuchos::ParameterList params;

  params.set("num_global_parts", numParts);
  params.set("partitioning_approach", "partition");
  params.set("algorithm", "rcb");
  params.set("keep_partition_tree", true);

  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  ////// Create and solve partitioning problem
  //////////////////////////////////////////////////////////////////////
  Zoltan2::PartitioningProblem<SparseMatrixAdapter_t> problem(&matAdapter, &params);

  try
  {
    if (me == 0) std::cout << "Calling solve() " << std::endl;
    problem.solve();
    if (me == 0) std::cout << "Done solve() " << std::endl;
  }
  catch (std::runtime_error &e)
  {
    std::cout << "Runtime exception returned from solve(): " << e.what();
    if (!strncmp(e.what(), "BUILD ERROR", 11)) {
      // Catching build errors as exceptions is OK in the tests
      std::cout << " PASS" << std::endl;
      return 0;
    }
    else {
      // All other runtime_errors are failures
      std::cout << " FAIL" << std::endl;
      return -1;
    }
  }
  catch (std::logic_error &e)
  {
    std::cout << "Logic exception returned from solve(): " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  catch (std::bad_alloc &e)
  {
    std::cout << "Bad_alloc exception returned from solve(): " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  catch (std::exception &e)
  {
    std::cout << "Unknown exception returned from solve(). " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  //////////////////////////////////////////////////////////////////////


  Zoltan2::PartitioningSolution<SparseMatrixAdapter_t> solution = problem.getSolution();

  bool binary = solution.isPartitioningTreeBinary();

  if (binary == false)
  {
    std::cout << "RCB should produce a binary partitioning tree. FAIL" << std::endl;
    return -1;
  }

  return analyze(problem, comm);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Test partitioning tree access for PHG partitioning.  This partitioning tree
// should be balanced.
////////////////////////////////////////////////////////////////////////////////
int testForPHG(SparseMatrixAdapter_t &matAdapter, int me, part_t numParts,
  RCP<tMVector_t> coords, RCP<const Teuchos::Comm<int> > comm)
{

  //////////////////////////////////////////////////////////////////////
  ////// Specify problem parameters
  //////////////////////////////////////////////////////////////////////
  Teuchos::ParameterList params;

  params.set("num_global_parts", numParts);
  params.set("partitioning_approach", "partition");
  params.set("algorithm", "zoltan");
  params.set("keep_partition_tree", true);

  Teuchos::ParameterList &zparams = params.sublist("zoltan_parameters",false);
  zparams.set("LB_METHOD","phg");
  zparams.set("FINAL_OUTPUT", "1");

  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  ////// Create and solve partitioning problem
  //////////////////////////////////////////////////////////////////////
  Zoltan2::PartitioningProblem<SparseMatrixAdapter_t> problem(&matAdapter, &params);

  try
  {
    if (me == 0) std::cout << "Calling solve() " << std::endl;
    problem.solve();
    if (me == 0) std::cout << "Done solve() " << std::endl;
  }
  catch (std::runtime_error &e)
  {
    std::cout << "Runtime exception returned from solve(): " << e.what();
    if (!strncmp(e.what(), "BUILD ERROR", 11)) {
      // Catching build errors as exceptions is OK in the tests
      std::cout << " PASS" << std::endl;
      return 0;
    }
    else {
      // All other runtime_errors are failures
      std::cout << " FAIL" << std::endl;
      return -1;
    }
  }
  catch (std::logic_error &e)
  {
    std::cout << "Logic exception returned from solve(): " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  catch (std::bad_alloc &e)
  {
    std::cout << "Bad_alloc exception returned from solve(): " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  catch (std::exception &e)
  {
    std::cout << "Unknown exception returned from solve(). " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  //////////////////////////////////////////////////////////////////////

  Zoltan2::PartitioningSolution<SparseMatrixAdapter_t> solution = problem.getSolution();

  bool binary = solution.isPartitioningTreeBinary();

  if (binary == false)
  {
    std::cout << "PHG should produce a binary partitioning tree. FAIL" << std::endl;
    return -1;
  }

  return analyze(problem, comm);
}
////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
// Test partitioning tree access for MJ partitioning.  This partitioning tree
// should be balanced.
////////////////////////////////////////////////////////////////////////////////
int testForMJ(SparseMatrixAdapter_t &matAdapter, int me, part_t numParts,
  RCP<tMVector_t> coords, RCP<const Teuchos::Comm<int> > comm)
{

  //////////////////////////////////////////////////////////////////////
  ////// Specify problem parameters
  //////////////////////////////////////////////////////////////////////
  Teuchos::ParameterList params;

  params.set("num_global_parts", numParts);
  params.set("algorithm", "multijagged");
  params.set("rectilinear", true);
  // params.set("mj_keep_part_boxes", true); // allows getPartBoxesView on solution
  params.set("keep_partition_tree", true);

  //////////////////////////////////////////////////////////////////////

  //////////////////////////////////////////////////////////////////////
  ////// Create and solve partitioning problem
  //////////////////////////////////////////////////////////////////////
  Zoltan2::PartitioningProblem<SparseMatrixAdapter_t> problem(&matAdapter, &params);

  try
  {
    if (me == 0) std::cout << "Calling solve() " << std::endl;
    problem.solve();
    if (me == 0) std::cout << "Done solve() " << std::endl;
  }
  catch (std::runtime_error &e)
  {
    std::cout << "Runtime exception returned from solve(): " << e.what();
    if (!strncmp(e.what(), "BUILD ERROR", 11)) {
      // Catching build errors as exceptions is OK in the tests
      std::cout << " PASS" << std::endl;
      return 0;
    }
    else {
      // All other runtime_errors are failures
      std::cout << " FAIL" << std::endl;
      return -1;
    }
  }
  catch (std::logic_error &e)
  {
    std::cout << "Logic exception returned from solve(): " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  catch (std::bad_alloc &e)
  {
    std::cout << "Bad_alloc exception returned from solve(): " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  catch (std::exception &e)
  {
    std::cout << "Unknown exception returned from solve(). " << e.what()
         << " FAIL" << std::endl;
    return -1;
  }
  //////////////////////////////////////////////////////////////////////

  Zoltan2::PartitioningSolution<SparseMatrixAdapter_t> solution = problem.getSolution();

  // copied from the MultiJaggedTest.cpp to inspect boxes
  // To activate set mj_keep_part_boxes true above
  /*
  std::vector<Zoltan2::coordinateModelPartBox<scalar_t, part_t> >
    & pBoxes = solution.getPartBoxesView();
  int coordDim = coords->getNumVectors();
  std::cout << std::endl;
  std::cout << "Plot final boxes..." << std::endl;
  for (size_t i = 0; i < pBoxes.size(); i++) {
    zscalar_t *lmin = pBoxes[i].getlmins();
    zscalar_t *lmax = pBoxes[i].getlmaxs();;
    int dim = pBoxes[i].getDim();
    std::set<int> * pNeighbors = pBoxes[i].getNeighbors();
    std::vector<int> * pGridIndices = pBoxes[i].getGridIndices();

    std::cout << me << " pBox " << i << " pid " << pBoxes[i].getpId()
              << " dim: " << dim
              << " (" << lmin[0] << "," << lmin[1] << ") "
              << "x"
              << " (" << lmax[0] << "," << lmax[1] << ")" << std::endl;
  }
  */

  bool binary = solution.isPartitioningTreeBinary();

  if (binary == true)
  {
    std::cout << "MJ should not produce a binary partitioning tree for this problem. FAIL" << std::endl;
    return -1;
  }

  return analyze(problem, comm);
}
////////////////////////////////////////////////////////////////////////////////




