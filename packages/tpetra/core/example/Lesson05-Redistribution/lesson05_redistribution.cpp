// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*!
\example lesson05_redistribution.cpp
\brief Tpetra parallel data redistribution, using Map and Export.

\ref Tpetra_Lesson05 explains this example in detail.
*/

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
createMatrix (const Teuchos::RCP<const typename CrsMatrixType::map_type>& map)
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
  typedef typename CrsMatrixType::scalar_type scalar_type;
  typedef typename CrsMatrixType::local_ordinal_type LO;
  typedef typename CrsMatrixType::global_ordinal_type GO;

  // Create a timer for sparse matrix creation.
  RCP<Time> timer = TimeMonitor::getNewCounter ("Sparse matrix creation");

  // Time the whole scope of this routine, not counting timer lookup.
  TimeMonitor monitor (*timer);

  // Create a Tpetra::Matrix using the Map, with dynamic allocation.
  RCP<CrsMatrixType> A (new CrsMatrixType (map, 3));

  // Add rows one at a time.  Off diagonal values will always be -1.
  const scalar_type two    = static_cast<scalar_type> ( 2.0);
  const scalar_type negOne = static_cast<scalar_type> (-1.0);

  const GST numGlobalIndices = map->getGlobalNumElements ();
  // const size_t numMyElements = map->getLocalNumElements ();

  // The list of global elements owned by this MPI process.
  ArrayView<const GO> myGlobalElements = map->getLocalElementList ();

  typedef typename ArrayView<const GO>::const_iterator iter_type;
  for (iter_type it = myGlobalElements.begin(); it != myGlobalElements.end(); ++it) {
    const LO i_local = *it;
    const GO i_global = map->getGlobalElement (i_local);

    // Can't insert local indices without a column map, so we insert
    // global indices here.
    if (i_global == 0) {
      A->insertGlobalValues (i_global,
                             tuple (i_global, i_global+1),
                             tuple (two, negOne));
    } else if (static_cast<GST> (i_global) == numGlobalIndices - 1) {
      A->insertGlobalValues (i_global,
                             tuple (i_global-1, i_global),
                             tuple (negOne, two));
    } else {
      A->insertGlobalValues (i_global,
                             tuple (i_global-1, i_global, i_global+1),
                             tuple (negOne, two, negOne));
    }
  }

  // Finish up the matrix.
  A->fillComplete ();
  return A;
}


void
example (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
         std::ostream& out)
{
  using std::endl;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  typedef Tpetra::global_size_t GST;

  // Set up Tpetra typedefs.
  typedef Tpetra::CrsMatrix<> crs_matrix_type;
  typedef Tpetra::Map<> map_type;
  typedef Tpetra::Map<>::global_ordinal_type global_ordinal_type;

  const int myRank = comm->getRank ();

  // The global number of rows in the matrix A to create.  We scale
  // this relative to the number of (MPI) processes, so that no matter
  // how many MPI processes you run, every process will have 10 rows.
  const GST numGlobalIndices = 10 * comm->getSize ();
  const global_ordinal_type indexBase = 0;

  // Construct a Map that is global (not locally replicated), but puts
  // all the equations on MPI Proc 0.
  if (myRank == 0) {
    out << "Construct Process 0 Map" << endl;
  }
  RCP<const map_type> procZeroMap;
  {
    const size_t numLocalIndices = (myRank == 0) ? numGlobalIndices : 0;
    procZeroMap = rcp (new map_type (numGlobalIndices, numLocalIndices,
                                     indexBase, comm));
  }

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  if (myRank == 0) {
    out << "Construct global Map" << endl;
  }
  RCP<const map_type> globalMap =
    rcp (new map_type (numGlobalIndices, indexBase, comm,
                       Tpetra::GloballyDistributed));

  // Create a sparse matrix using procZeroMap.
  if (myRank == 0) {
    out << "Create sparse matrix using Process 0 Map" << endl;
  }
  RCP<const crs_matrix_type> A = createMatrix<crs_matrix_type> (procZeroMap);

  //
  // We've created a sparse matrix that lives entirely on Process 0.
  // Now we want to distribute it over all the processes.
  //

  // Redistribute the matrix.  Since both the source and target Maps
  // are one-to-one, we could use either an Import or an Export.  If
  // only the source Map were one-to-one, we would have to use an
  // Import; if only the target Map were one-to-one, we would have to
  // use an Export.  We do not allow redistribution using Import or
  // Export if neither source nor target Map is one-to-one.
  RCP<crs_matrix_type> B;
  if (myRank == 0) {
    out << "Redistribute sparse matrix" << endl;
  }
  {
    // We created exportTimer in main().  It's a global timer.
    // Actually starting and stopping the timer is local, but
    // computing timer statistics (e.g., in TimeMonitor::summarize(),
    // called in main()) is global.  There are ways to restrict the
    // latter to any given MPI communicator; the default is
    // MPI_COMM_WORLD.
    TimeMonitor monitor (*exportTimer); // Time the redistribution

    // Make an export object with procZeroMap as the source Map, and
    // globalMap as the target Map.  The Export type has the same
    // template parameters as a Map.  Note that Export does not depend
    // on the Scalar template parameter of the objects it
    // redistributes.  You can reuse the same Export for different
    // Tpetra object types, or for Tpetra objects of the same type but
    // different Scalar template parameters (e.g., Scalar=float or
    // Scalar=double).
    typedef Tpetra::Export<> export_type;
    export_type exporter (procZeroMap, globalMap);

    // Make a new sparse matrix whose row map is the global Map.
    B = rcp (new crs_matrix_type (globalMap, 0));

    // Redistribute the data, NOT in place, from matrix A (which lives
    // entirely on Proc 0) to matrix B (which is distributed evenly over
    // the processes).
    B->doExport (*A, exporter, Tpetra::INSERT);
  }

  // We time redistribution of B separately from fillComplete().
  B->fillComplete ();
}


int
main (int argc, char *argv[])
{
  using Teuchos::TimeMonitor;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    auto comm = Tpetra::getDefaultComm ();
    const int myRank = comm->getRank ();
    // const int numProcs = comm->getSize ();

    // Make global timer for sparse matrix redistribution.
    // We will use (start and stop) this timer in example().
    exportTimer =
      TimeMonitor::getNewCounter ("Sparse matrix redistribution");
    example (comm, std::cout); // Run the whole example.

    // Summarize global performance timing results, for all timers
    // created using TimeMonitor::getNewCounter().
    TimeMonitor::summarize (std::cout);

    // Make sure that the timer goes away before main() exits.
    exportTimer = Teuchos::null;

    // This tells the Trilinos test framework that the test passed.
    if (myRank == 0) {
      std::cout << "End Result: TEST PASSED" << std::endl;
    }
  }
  return 0;
}
