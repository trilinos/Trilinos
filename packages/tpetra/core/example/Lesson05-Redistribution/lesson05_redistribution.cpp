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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
// @HEADER

/*!
\example lesson05_redistribution.cpp
\brief Tpetra parallel data redistribution, using Map and Export.

\ref Tpetra_Lesson05 explains this example in detail.
*/

#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>
#include <Teuchos_TimeMonitor.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Version.hpp>
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
  // const size_t numMyElements = map->getNodeNumElements ();

  // The list of global elements owned by this MPI process.
  ArrayView<const GO> myGlobalElements = map->getNodeElementList ();

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
         std::ostream& out,
         std::ostream& err)
{
  using std::endl;
  using Teuchos::ParameterList;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;
  typedef Tpetra::global_size_t GST;

  // Set up Tpetra typedefs.
  typedef double scalar_type;
  typedef Tpetra::CrsMatrix<scalar_type> crs_matrix_type;
  typedef Tpetra::Map<> map_type;
  typedef Tpetra::Map<>::global_ordinal_type global_ordinal_type;

  // Print out the Tpetra software version information.
  out << Tpetra::version () << endl << endl;

  // The global number of rows in the matrix A to create.  We scale
  // this relative to the number of (MPI) processes, so that no matter
  // how many MPI processes you run, every process will have 10 rows.
  const GST numGlobalIndices = 10 * comm->getSize ();
  const global_ordinal_type indexBase = 0;

  // Construct a Map that is global (not locally replicated), but puts
  // all the equations on MPI Proc 0.
  RCP<const map_type> procZeroMap;
  {
    const int myRank = comm->getRank ();
    const size_t numLocalIndices = (myRank == 0) ? numGlobalIndices : 0;
    procZeroMap = rcp (new map_type (numGlobalIndices, numLocalIndices,
                                     indexBase, comm));
  }

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  RCP<const map_type> globalMap =
    rcp (new map_type (numGlobalIndices, indexBase, comm,
                       Tpetra::GloballyDistributed));

  // Create a sparse matrix using procZeroMap.
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
  using Teuchos::RCP;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

  const int myRank = comm->getRank ();
  // const int numProcs = comm->getSize ();
  std::ostream& out = (myRank == 0) ? std::cout : blackHole;
  std::ostream& err = (myRank == 0) ? std::cerr : blackHole;

  // Make global timer for sparse matrix redistribution.
  // We will use (start and stop) this timer in example().
  exportTimer = TimeMonitor::getNewCounter ("Sparse matrix redistribution");

  example (comm, out, err); // Run the whole example.

  // Summarize global performance timing results, for all timers
  // created using TimeMonitor::getNewCounter().
  TimeMonitor::summarize (out);

  // Make sure that the timer goes away before main() exits.
  exportTimer = Teuchos::null;

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }

  return 0;
}
