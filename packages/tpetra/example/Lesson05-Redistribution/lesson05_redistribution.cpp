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


template<class TpetraMatrixType>
class MatrixFactory {
public:
  // Fetch typedefs from the Tpetra::CrsMatrix.
  typedef typename TpetraMatrixType::scalar_type scalar_type;
  typedef typename TpetraMatrixType::local_ordinal_type local_ordinal_type;
  typedef typename TpetraMatrixType::global_ordinal_type global_ordinal_type;
  typedef typename TpetraMatrixType::node_type node_type;

  // The type of Map objects to use.
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

  // Create and return a simple example CrsMatrix, with row
  // distribution over the given Map.
  Teuchos::RCP<const TpetraMatrixType>
  create (const Teuchos::RCP<const map_type>& map) const
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

    // Create a timer for sparse matrix creation.
    RCP<Time> timer = TimeMonitor::getNewCounter ("Sparse matrix creation");

    // Time the whole scope of this routine, not counting timer lookup.
    TimeMonitor monitor (*timer);

    // Create a Tpetra::Matrix using the Map, with dynamic allocation.
    RCP<TpetraMatrixType> A = rcp (new TpetraMatrixType (map, 3));

    // Add rows one at a time.  Off diagonal values will always be -1.
    const scalar_type two    = static_cast<scalar_type>( 2.0);
    const scalar_type negOne = static_cast<scalar_type>(-1.0);

    const GST numGlobalElements = map->getGlobalNumElements ();
    //    const size_t numMyElements = map->getNodeNumElements ();

    // The list of global elements owned by this MPI process.
    ArrayView<const global_ordinal_type> myGlobalElements =
      map->getNodeElementList ();

    typedef typename ArrayView<const global_ordinal_type>::const_iterator iter_type;
    for (iter_type it = myGlobalElements.begin(); it != myGlobalElements.end(); ++it) {
      const local_ordinal_type i_local = *it;
      const global_ordinal_type i_global = map->getGlobalElement (i_local);

      // Can't insert local indices without a column map, so we insert
      // global indices here.
      if (i_global == 0) {
        A->insertGlobalValues (i_global,
                              tuple (i_global, i_global+1),
                              tuple (two, negOne));
      } else if (static_cast<GST> (i_global) == numGlobalElements - 1) {
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
};


template<class NodeType>
void
example (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
         const Teuchos::RCP<NodeType>& node,
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

  // Print out the Tpetra software version information.
  out << Tpetra::version() << endl << endl;

  // Set up Tpetra typedefs.
  typedef double scalar_type;
  typedef int local_ordinal_type;
  typedef long global_ordinal_type;
  typedef NodeType node_type;
  typedef Tpetra::CrsMatrix<scalar_type, local_ordinal_type,
                            global_ordinal_type, node_type> matrix_type;

  // The type of the Tpetra::Map that describes how the matrix is distributed.
  typedef Tpetra::Map<local_ordinal_type, global_ordinal_type, node_type> map_type;

  // The global number of rows in the matrix A to create.  We scale
  // this relative to the number of (MPI) processes, so that no matter
  // how many MPI processes you run, every process will have 10 rows.
  const GST numGlobalElements = 10 * comm->getSize ();

  const global_ordinal_type indexBase = 0;

  // Construct a Map that is global (not locally replicated), but puts
  // all the equations on MPI Proc 0.
  RCP<const map_type> procZeroMap;
  {
    const size_t numLocalElements =
      (comm->getRank() == 0) ? numGlobalElements : 0;
    procZeroMap = rcp (new map_type (numGlobalElements, numLocalElements,
                                     indexBase, comm, node));
  }

  // Construct a Map that puts approximately the same number of
  // equations on each processor.
  RCP<const map_type> globalMap =
    rcp (new map_type (numGlobalElements, indexBase, comm,
                       Tpetra::GloballyDistributed, node));

  // Create a sparse matrix using procZeroMap.
  RCP<const matrix_type> A = MatrixFactory<matrix_type> ().create (procZeroMap);

  //
  // We've created a sparse matrix that lives entirely on MPI Proc 0.
  // Now we want to distribute it over all the processes.
  //

  // Make a timer for sparse matrix redistribution.
  RCP<Time> exportTimer =
    TimeMonitor::getNewCounter ("Sparse matrix redistribution");

  // Redistribute the matrix.  Since both the source and target Maps
  // are one-to-one, we could use either an Import or an Export.  If
  // only the source Map were one-to-one, we would have to use an
  // Import; if only the target Map were one-to-one, we would have to
  // use an Export.  We do not allow redistribution using Import or
  // Export if neither source nor target Map is one-to-one.
  RCP<matrix_type> B;
  {
    TimeMonitor monitor (*exportTimer); // Time the redistribution

    // Make an export object with procZeroMap as the source Map, and
    // globalMap as the target Map.  The Export type has the same
    // template parameters as a Map.  Note that Export does not depend
    // on the Scalar template parameter of the objects it
    // redistributes.  You can reuse the same Export for different
    // Tpetra object types, or for Tpetra objects of the same type but
    // different Scalar template parameters (e.g., Scalar=float or
    // Scalar=double).
    typedef Tpetra::Export<local_ordinal_type, global_ordinal_type,
                           node_type> export_type;
    export_type exporter (procZeroMap, globalMap);

    // Make a new sparse matrix whose row map is the global Map.
    B = rcp (new matrix_type (globalMap, 0));

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
  using std::endl;
  using Teuchos::RCP;
  using Teuchos::Time;
  using Teuchos::TimeMonitor;

  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  typedef Kokkos::DefaultNode::DefaultNodeType node_type;
  RCP<node_type> node = Kokkos::DefaultNode::getDefaultNode ();

  const int myRank = comm->getRank();
  //const int numProcs = comm->getSize();
  std::ostream& out = (myRank == 0) ? std::cout : blackHole;
  std::ostream& err = (myRank == 0) ? std::cerr : blackHole;

  // Make a timer for sparse matrix redistribution.
  //
  // If you are using Trilinos 10.6 instead of the development branch
  // (10.7), just delete this line of code, and make the other change
  // mentioned above.
  RCP<Time >exportTimer =
    TimeMonitor::getNewCounter ("Sparse matrix redistribution");

  // Run the whole example: create the sparse matrix, and compute the preconditioner.
  example (comm, node, out, err);

  // Summarize global performance timing results, for all timers
  // created using TimeMonitor::getNewCounter().
  TimeMonitor::summarize (out);

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }

  return 0;
}
