// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*!
\example lesson02_read_modify_vec.cpp
\brief Read and modify the entries of a vector (Tpetra::Vector),
  using local indices.

\ref Tpetra_Lesson02 explains this example in detail.
*/

#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_CommHelpers.hpp>

void
exampleRoutine (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                std::ostream& out)
{
  using std::endl;
  using Teuchos::Array;
  using Teuchos::ArrayRCP;
  using Teuchos::ArrayView;
  using Teuchos::outArg;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using Teuchos::REDUCE_SUM;
  using Teuchos::reduceAll;

  const int myRank = comm->getRank ();

  // Print out the Tpetra software version information.
  if (myRank == 0) {
    out << Tpetra::version () << endl << endl;
  }

  // Type of the Tpetra::Map specialization to use.
  using map_type = Tpetra::Map<>;

  // The type of the Tpetra::Vector specialization to use.
  using vector_type = Tpetra::Vector<>;

  // The "Scalar" type is the type of the values stored in the Tpetra::Vector.
  using scalar_type = Tpetra::Vector<>::scalar_type;

  // The "LocalOrdinal" (LO) type is the type of "local" indices.
  // The typedef is commented out to avoid "unused typedef" warnings.
  //
  //using local_ordinal_type = Tpetra::Vector<>::local_ordinal_type;

  // The "GlobalOrdinal" (GO) type is the type of "global" indices.
  using global_ordinal_type = Tpetra::Vector<>::global_ordinal_type;

  //////////////////////////////////////////////////////////////////////
  // Create a Tpetra Map
  //////////////////////////////////////////////////////////////////////

  // The total (global, i.e., over all MPI processes) number of
  // entries in the Map.
  //
  // For this example, we scale the global number of entries in the
  // Map with the number of MPI processes.  That way, you can run this
  // example with any number of MPI processes and every process will
  // still have a positive number of entries.
  const Tpetra::global_size_t numGlobalEntries = comm->getSize () * 5;

  // Index base of the Map.  We choose zero-based (C-style) indexing.
  const global_ordinal_type indexBase = 0;

  // Construct a Map that puts the same number of equations on each
  // MPI process.
  RCP<const map_type> contigMap =
    rcp (new map_type (numGlobalEntries, indexBase, comm,
                       Tpetra::GloballyDistributed));

  //////////////////////////////////////////////////////////////////////
  // Create a Tpetra Vector
  //////////////////////////////////////////////////////////////////////

  // Create a Vector with the Map we created above.
  // This version of the constructor will fill in the vector with zeros.
  vector_type x (contigMap);

  //////////////////////////////////////////////////////////////////////
  // Fill the Vector with a single number, or with random numbers
  //////////////////////////////////////////////////////////////////////

  // Set all entries of x to 42.0.
  x.putScalar (42.0);

  // norm2() is a collective, so we need to call it on all processes
  // in the Vector's communicator.
  auto x_norm2 = x.norm2 ();
  if (myRank == 0) {
    out << "Norm of x (all entries are 42.0): " << x_norm2 << endl;
  }

  // Set the entries of x to (pseudo)random numbers.  Please don't
  // consider this a good parallel pseudorandom number generator.
  x.randomize ();

  // Print the norm of x.
  x_norm2 = x.norm2 ();
  if (myRank == 0) {
    out << "Norm of x (random numbers): " << x_norm2 << endl;
  }

  //////////////////////////////////////////////////////////////////////
  // Read the entries of the Vector
  //////////////////////////////////////////////////////////////////////

  {
    // Get a const persisting view of the entries in the Vector.
    // "Const" means that you cannot modify the entries.  "Persisting"
    // means that the view persists beyond the lifetime of the Vector.
    // Even after the Vector's destructor is called, the view won't go
    // away.  If the Vector's entries change during the lifetime of a
    // persisting view, the entries of the persisting view are
    // undefined.
    //
    // An ArrayRCP acts like an array or std::vector, but is
    // reference-counted like std::shared_ptr or Teuchos::RCP.  You
    // may decrement the reference count manually by assigning
    // Teuchos::null to it.  We put this code in an inner scope (in an
    // extra pair of {}) so that the ArrayRCP will fall out of scope
    // (and therefore be deallocated) before the next example, which
    // modifies the entries of the Vector.
    ArrayRCP<const scalar_type> x_data = x.get1dView ();

    // x_data.size () may be longer than the number of local rows in
    // the Vector, so be sure to ask the Vector for its dimensions,
    // rather than the ArrayRCP.
    const size_t localLength = x.getLocalLength ();

    // Count the local number of entries less than 0.5.
    // Use local indices to access the entries of x_data.
    size_t localCount = 0;
    for (size_t k = 0; k < localLength; ++k) {
      if (x_data[k] < 0.5) {
        ++localCount;
      }
    }

    // "reduceAll" is a type-safe templated version of MPI_Allreduce.
    // "outArg" is like taking the address using &, but makes it more
    // clear that its argument is an output argument of a function.
    size_t globalCount = 0;
    reduceAll<int, size_t> (*comm, REDUCE_SUM, localCount,
                            outArg (globalCount));

    // Find the total number of entries less than 0.5, over all
    // processes in the Vector's communicator.  Note the trick for
    // pluralizing the word "entry" conditionally on globalCount.
    if (myRank == 0) {
      out << "x has " << globalCount << " entr"
          << (globalCount != 1 ? "ies" : "y")
          << " less than 0.5." << endl;
    }
  }

  //////////////////////////////////////////////////////////////////////
  // Modify the entries of the Vector
  //////////////////////////////////////////////////////////////////////

  {
    // Get a nonconst persisting view of the entries in the Vector.
    // "Nonconst" means that you may modify the entries.
    // "Persisting" means that the view persists beyond the lifetime of the Vector.
    // Even after the Vector's destructor is called, the view won't go away.
    // If you create two nonconst persisting views of the same Vector,
    // and modify the entries of one view during the lifetime of the other view,
    // the entries of the other view are undefined.
    ArrayRCP<scalar_type> x_data = x.get1dViewNonConst ();

    // Use local indices to access the entries of x_data.
    // x_data.size () may be longer than the number of local rows in the Vector,
    // so be sure to ask the Vector for its dimensions, rather than the ArrayRCP.
    const size_t localLength = x.getLocalLength ();
    for (size_t k = 0; k < localLength; ++k) {
      // Add k (the local index) to every entry of x.  Treat
      // scalar_type as a function to convert k to scalar_type.
      x_data[k] += scalar_type (k);
    }

    // Print the norm of x.
    x_norm2 = x.norm2 ();
    if (myRank == 0) {
      out << "Norm of x (modified random numbers): " << x_norm2 << endl;
    }
  }
}

//
// The same main() driver routine as in the first Tpetra lesson.
//
int
main (int argc, char *argv[])
{
  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    auto comm = Tpetra::getDefaultComm ();
    exampleRoutine (comm, std::cout);
    // This tells the Trilinos test framework that the test passed.
    if (comm->getRank () == 0) {
      std::cout << "End Result: TEST PASSED" << std::endl;
    }
  }
  return 0;
}
