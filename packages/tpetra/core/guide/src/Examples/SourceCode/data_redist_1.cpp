// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Tpetra_Core.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <iostream>
// Timer for use in example().
Teuchos::RCP<Teuchos::Time> exportTimer;

// Create and return a simple example CrsMatrix, with row distribution
// over the given Map.
//
// CrsMatrixType: The type of the Tpetra::CrsMatrix specialization to use.
template<class CrsMatrixType>
Teuchos::RCP<const CrsMatrixType>
createMatrix(const Teuchos::RCP<const typename CrsMatrixType::map_type>& map)
{
  using Teuchos::arcp;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using Teuchos::tuple;
  typedef Tpetra::global_size_t GST;
  // Fetch typedefs from the Tpetra::CrsMatrix.
  using SC = typename CrsMatrixType::scalar_type;
  using GO = typename CrsMatrixType::global_ordinal_type;
  using LO = typename CrsMatrixType::local_ordinal_type;

  // Create a timer for sparse matrix creation.
  RCP<Time> timer = TimeMonitor::getNewCounter("Sparse matrix creation");

  // Time the whole scope of this routine, not counting timer lookup.
  TimeMonitor monitor(*timer);

  // Create a Tpetra::Matrix using the Map, with dynamic allocation.
  RCP<CrsMatrixType> A(new CrsMatrixType(map, 3));

  // Add rows one at a time.  Off diagonal values will always be -1.
  const SC two    = static_cast<SC>( 2.0);
  const SC negOne = static_cast<SC>(-1.0);
  const GST numGlobalIndices = map->getGlobalNumElements();

  // const size_t numMyElements = map->getLocalNumElements();
  // The list of global elements owned by this MPI process.
  ArrayView<const GO> myGlobalElements = map->getLocalElementList();
  typedef typename ArrayView<const GO>::const_iterator iter_type;
  for(iter_type it = myGlobalElements.begin(); it != myGlobalElements.end(); ++it) {
    const LO i_local = *it;
    const GO i_global = map->getGlobalElement(i_local);
    // Can't insert local indices without a column map, so we insert
    // global indices here.
    if(i_global == 0) {
      A->insertGlobalValues(i_global,
                            tuple(i_global, i_global+1),
                            tuple(two, negOne));
    } else if(static_cast<GST>(i_global) == numGlobalIndices - 1) {
      A->insertGlobalValues(i_global,
                            tuple(i_global-1, i_global),
                            tuple(negOne, two));
    } else {
      A->insertGlobalValues(i_global,
                            tuple(i_global-1, i_global, i_global+1),
                            tuple(negOne, two, negOne));
    }
  }
  // Finish up the matrix.
  A->fillComplete();
  return A;
}

void
example(const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
        std::ostream& out, std::ostream& err)
{
  using std::endl;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  using GST = Tpetra::global_size_t;

  // Set up Tpetra typedefs.
  using SC = Tpetra::CrsMatrix<>::scalar_type;
  using crs_matrix_type = Tpetra::CrsMatrix<SC>;
  using map_type = Tpetra::Map<>;
  using GO = Tpetra::Map<>::global_ordinal_type;

  // The global number of rows in the matrix A to create.
  const GST numGlobalIndices = 10 * comm->getSize();
  const GO indexBase = 0;

  // Construct a Map that is global (not locally replicated), but puts
  // all the equations on MPI Proc 0.
  RCP<const map_type> procZeroMap;
  {
    const int myRank = comm->getRank();
    const size_t numLocalIndices = (myRank == 0) ? numGlobalIndices : 0;
    procZeroMap = rcp(new map_type(numGlobalIndices, numLocalIndices, indexBase, comm));
  }

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  RCP<const map_type> globalMap =
    rcp(new map_type(numGlobalIndices, indexBase, comm, Tpetra::GloballyDistributed));

  // Create a sparse matrix using procZeroMap.
  RCP<const crs_matrix_type> A = createMatrix<crs_matrix_type>(procZeroMap);

  // We've created a sparse matrix that lives entirely on Process 0.
  // Now we want to distribute it over all the processes.
  //
  // Redistribute the matrix.
  RCP<crs_matrix_type> B;
  {
    // We created exportTimer in main().
    TimeMonitor monitor(*exportTimer); // Time the redistribution

    // Make an export object with procZeroMap as the source Map, and
    // globalMap as the target Map.
    typedef Tpetra::Export<> export_type;
    export_type exporter(procZeroMap, globalMap);

    // Make a new sparse matrix whose row map is the global Map.
    B = rcp(new crs_matrix_type(globalMap, 0));

    // Redistribute the data, NOT in place, from matrix A (which lives
    // entirely on Proc 0) to matrix B (which is distributed evenly over
    // the processes).
    B->doExport(*A, exporter, Tpetra::INSERT);
  }

  // We time redistribution of B separately from fillComplete().
  B->fillComplete();
}

int
main(int argc, char *argv[])
{
  using Teuchos::RCP;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  Teuchos::oblackholestream blackHole;

  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  {
    // Get a communicator corresponding to MPI_COMM_WORLD
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

    const int myRank = comm->getRank();

    std::ostream& out = (myRank == 0) ? std::cout : blackHole;
    std::ostream& err = (myRank == 0) ? std::cerr : blackHole;

    // Make global timer for sparse matrix redistribution.
    // We will use (start and stop) this timer in example().
    exportTimer = TimeMonitor::getNewCounter("Sparse matrix redistribution");
    example(comm, out, err); // Run the whole example.

    // Summarize global performance timing results, for all timers
    // created using TimeMonitor::getNewCounter().
    TimeMonitor::summarize(out);

    // Make sure that the timer goes away before main() exits.
    exportTimer = Teuchos::null;

    // This tells the Trilinos test framework that the test passed.
    if(myRank == 0) {
      std::cout << "End Result: TEST PASSED" << std::endl;
    }
    // ScopeGuard's destructor calls MPI_Finalize, if its constructor
    // called MPI_Init.  Likewise, it calls Kokkos::finalize, if its
    // constructor called Kokkos::initialize.
  }
  return 0;
}
