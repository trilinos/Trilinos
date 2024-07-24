// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/*!
\example lesson02_init_map_vec.cpp
\brief Create data distributions (Tpetra::Map), and create vectors
  (Tpetra::Vector) distributed according to those distributions.

\ref Tpetra_Lesson02 explains this example in detail.
*/

#include <Tpetra_Core.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_Comm.hpp>
#include <Teuchos_OrdinalTraits.hpp>

void
exampleRoutine (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                std::ostream& out)
{
  using std::endl;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  const int myRank = comm->getRank ();
  if (myRank == 0) {
    // Print out the Tpetra software version information.
    out << Tpetra::version () << endl << endl;
  }

  //
  // The first thing you might notice that makes Tpetra objects
  // different than their Epetra counterparts, is that Tpetra objects
  // take several template parameters.  These template parameters give
  // Tpetra its features of being able to solve very large problems
  // (of more than 2 billion unknowns) and to exploit intranode
  // parallelism.
  //
  // Most Tpetra objects come with default values of those template
  // parameters.  In many cases, you might not have to specify _any_
  // of those values explicitly!  You can also control some of their
  // default values, like that of the Node type, when configuring
  // Trilinos.
  //
  // It's common to begin a Tpetra application with some typedefs to
  // make the code more concise and readable.  They also make the code
  // more maintainable, since you can change the typedefs without
  // changing the rest of the program.  You may also want to
  // abbreviate the template parameters.  Standard abbreviations are
  // "LO" for local_ordinal_type and "GO" for global_ordinal_type.
  //

  // The "Scalar" type is the type of the values stored in the Tpetra
  // objects.  Valid Scalar types include real or complex
  // (std::complex<T>) floating-point types, or more exotic objects
  // with similar behavior.  We use the default type here, which we
  // get from Vector.  "Vector<>" means that we let all template
  // parameters' values revert to their defaults.  The default type is
  // 'double', a 64-bit double-precision binary floating-point value.
  typedef Tpetra::Vector<>::scalar_type scalar_type;

  // The "LocalOrdinal" (LO) type is the type of "local" indices.
  // Both Epetra and Tpetra index local entries differently than
  // global entries.  Tpetra exploits this so that you can use a
  // shorter integer type for local indices.  This saves bandwidth
  // when computing sparse matrix-vector products.  We use the default
  // LO type here, which we get from Tpetra::Vector.  We could also
  // get it from Tpetra::Map.

  // This line is commented out because we don't actually use this
  // type in the code below.  Leaving the typedef in that case will
  // make the compiler emit "unused typedef" warnings.
  //
  //typedef Tpetra::Vector<>::local_ordinal_type local_ordinal_type;

  // The "GlobalOrdinal" (GO) type is the type of "global" indices.
  // We use the default GO type here, which we get from
  // Tpetra::Vector.  We could also get it from Tpetra::Map.
  typedef Tpetra::Vector<>::global_ordinal_type global_ordinal_type;

  // The Kokkos "Node" type describes the type of shared-memory
  // parallelism that Tpetra will use _within_ an MPI process.  The
  // available Node types depend on Trilinos' build options and the
  // availability of certain third-party libraries.  In almost all
  // cases, the default setting will do.  You may set the default Node
  // type when configuring Trilinos.  In this case, we access the
  // default Node type using the typedef in Tpetra::Vector.  Almost
  // all Tpetra classes have default template parameter values.

  // This line is commented out because we don't actually use this
  // type in the code below.  Leaving the typedef in that case will
  // make the compiler emit "unused typedef" warnings.
  //typedef Tpetra::Vector<>::node_type node_type;

  // Maps know how to convert between local and global indices, so of
  // course they are templated on the local and global Ordinal types.
  // They are also templated on the Kokkos Node type, because Tpetra
  // objects that use Tpetra::Map are.  It's important not to mix up
  // Maps for different Kokkos Node types.  In this case, we use all
  // default template parameters, which are the same as the
  // corresponding template parameters of Vector.
  typedef Tpetra::Map<> map_type;

  //////////////////////////////////////////////////////////////////////
  // Create some Tpetra Map objects
  //////////////////////////////////////////////////////////////////////

  //
  // Like Epetra, Tpetra has local and global Maps.  Local maps
  // describe objects that are replicated over all participating MPI
  // processes.  Global maps describe distributed objects.  You can do
  // imports and exports between local and global maps; this is how
  // you would turn locally replicated objects into distributed
  // objects and vice versa.
  //

  // numLocalEntries: The local (on the calling MPI process) number of
  // entries (indices) in the first Map that we create.  Tpetra
  // expects a size_t for this value.
  const size_t numLocalEntries = 5;

  // numGlobalEntries: The total (global, i.e., over all MPI
  // processes) number of entries (indices) in the Map.  Tpetra
  // expects Tpetra::global_size_t for this value.  This type is at
  // least 64 bits long on 64-bit machines.
  //
  // For this example, we scale the global number of entries in the
  // Map with the number of MPI processes.  That way, you can run this
  // example with any number of MPI processes and every process will
  // still have a positive number of entries.
  const Tpetra::global_size_t numGlobalEntries =
    comm->getSize () * numLocalEntries;

  // Tpetra can index the entries of a Map starting with 0 (C style),
  // 1 (Fortran style), or any base you want.  1-based indexing is
  // handy when interfacing with Fortran.  We choose 0-based indexing
  // here.
  const global_ordinal_type indexBase = 0;

  //
  // Create some Maps.  All Map constructors must be called as a
  // collective over the input communicator.  Not all Map constructors
  // necessarily require communication, but some do, so it's best to
  // treat them all as collectives.
  //

  // Construct a Map that puts the same number of equations on each
  // processor.  The resulting Map is "contiguous and uniform."
  //
  // Maps should be considered immutable objects.  This is why we
  // create it as a "const map_type".  If you want a new data
  // distribution, create a new Map.
  RCP<const map_type> contigMap =
    rcp (new map_type (numGlobalEntries, indexBase, comm));

  // contigMap is contiguous by construction.  Test this at run time.
  // Lesson 01 introduced the TEUCHOS_TEST_FOR_EXCEPTION macro, which
  // throws an exception of the given type (second argument) with the
  // given message (third argument), if the first argument is true.
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! contigMap->isContiguous (), std::logic_error,
    "The supposedly contiguous Map isn't contiguous.");

  // contigMap2: Create a Map which is the same as contigMap, but uses
  // a different Map constructor.  This one asks for the number of
  // entries on each MPI process.  The resulting Map is "contiguous"
  // but not necessarily uniform, since the numbers of entries on
  // different MPI processes may differ.  In this case, the number of
  // entries on each MPI process is the same, but that doesn't always
  // have to be the case.
  RCP<const map_type> contigMap2 =
    rcp (new map_type (numGlobalEntries, numLocalEntries, indexBase, comm));

  // Since contigMap and contigMap2 have the same communicators, and
  // the same number of entries on all MPI processes in their
  // communicators, they are "the same."
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! contigMap->isSameAs (*contigMap2), std::logic_error,
    "contigMap should be the same as contigMap2, but it's not.");

  // contigMap3: Use the same Map constructor as contigMap3, but don't
  // specify the global number of entries.  This is helpful if you
  // only know how many entries each MPI process has, but don't know
  // the global number.  Instead of numGlobalEntries, we use the
  // equivalent of Epetra's -1 for Tpetra::global_size_t (which might
  // be unsigned, so don't use -1!!!), which we call "INVALID" (an
  // "invalid value" used as a flag).
  const Tpetra::global_size_t INVALID =
    Teuchos::OrdinalTraits<Tpetra::global_size_t>::invalid ();
  RCP<const map_type> contigMap3 =
    rcp (new map_type (INVALID, numLocalEntries, indexBase, comm));

  // Even though we made contigMap3 without specifying the global
  // number of entries, it should still be the same as contigMap2.
  TEUCHOS_TEST_FOR_EXCEPTION(
    ! contigMap2->isSameAs (*contigMap3), std::logic_error,
    "contigMap2 should be the same as contigMap3, but it's not.");

  // Create a Map which has the same number of global entries per
  // process as contigMap, but distributes them differently, in
  // round-robin (1-D cyclic) fashion instead of contiguously.
  RCP<const map_type> cyclicMap;
  {
    // We'll use the version of the Map constructor that takes, on
    // each MPI process, a list of the global entries in the Map
    // belonging to that process.  You can use this constructor to
    // construct an overlapping (also called "not 1-to-1") Map, in
    // which one or more entries are owned by multiple processes.  We
    // don't do that here; we make a nonoverlapping (also called
    // "1-to-1") Map.
    Array<global_ordinal_type>::size_type numEltsPerProc = 5;
    Array<global_ordinal_type> elementList (numEltsPerProc);

    const int numProcs = comm->getSize ();
    for (Array<global_ordinal_type>::size_type k = 0; k < numEltsPerProc; ++k) {
      elementList[k] = myRank + k*numProcs;
    }

    cyclicMap = rcp (new map_type (numGlobalEntries, elementList, indexBase, comm));
  }

  // If there's more than one MPI process in the communicator,
  // then cyclicMap is definitely NOT contiguous.
  TEUCHOS_TEST_FOR_EXCEPTION(
    comm->getSize () > 1 && cyclicMap->isContiguous (),
    std::logic_error,
    "The cyclic Map claims to be contiguous.");

  // contigMap and cyclicMap should always be compatible.  However, if
  // the communicator contains more than 1 process, then contigMap and
  // cyclicMap are NOT the same.
  TEUCHOS_TEST_FOR_EXCEPTION(! contigMap->isCompatible (*cyclicMap),
    std::logic_error,
    "contigMap should be compatible with cyclicMap, but it's not.");
  TEUCHOS_TEST_FOR_EXCEPTION(comm->getSize() > 1 && contigMap->isSameAs (*cyclicMap),
    std::logic_error,
    "contigMap should be compatible with cyclicMap, but it's not.");

  //////////////////////////////////////////////////////////////////////
  // We have maps now, so we can create vectors.
  //////////////////////////////////////////////////////////////////////

  // Since Tpetra::Vector takes four template parameters, its type is
  // long.  Here, we use all default values for the template
  // parameters, which helps with the length.  However, I still prefer
  // to use a typedef for encapsulation, so that we only have to
  // change one line of code if we decide to change the template
  // parameters of Vector.
  typedef Tpetra::Vector<> vector_type;

  // Create a Vector with the contiguous Map.  This version of the
  // constructor will fill in the vector with zeros.
  vector_type x (contigMap);

  // The two-argument copy constructor with second argument
  // Teuchos::Copy performs a deep copy.  x and y have the same Map.
  // The one-argument copy constructor does a _shallow_ copy.
  vector_type y (x, Teuchos::Copy);

  // Create a Vector with the 1-D cyclic Map.  Calling the constructor
  // with false for the second argument leaves the data uninitialized,
  // so that you can fill it later without paying the cost of
  // initially filling it with zeros.
  vector_type z (cyclicMap, false);

  // Set the entries of z to (pseudo)random numbers.  Please don't
  // consider this a good parallel pseudorandom number generator.
  z.randomize ();

  // Set the entries of x to all ones.
  //
  // The code below works because scalar_type=double.  In general, you
  // may use the commented-out line of code, if the conversion from
  // float to scalar_type is not defined for your scalar type.
  x.putScalar (1.0);
  //x.putScalar (Teuchos::ScalarTraits<scalar_type>::one());

  // See comment above about type conversions to scalar_type.
  const scalar_type alpha = 3.14159;
  const scalar_type beta = 2.71828;
  const scalar_type gamma = -10.0;

  // x = beta*x + alpha*z
  //
  // This is a legal operation!  Even though the Maps of x and z are
  // not the same, their Maps are compatible.  Whether it makes sense
  // or not depends on your application.
  x.update (alpha, z, beta);

  // See comment above about type conversions from float to scalar_type.
  y.putScalar (42.0);
  // y = gamma*y + alpha*x + beta*z
  y.update (alpha, x, beta, z, gamma);

  // Compute the 2-norm of y.
  //
  // The norm may have a different type than scalar_type.  For
  // example, if scalar_type is complex, then the norm is real.
  // Tpetra::MultiVector and Tpetra::Vector give us the type of the
  // norm.
  //
  // If you are using an older version of Tpetra, this code might not
  // work.  Try the commented-out line instead in that case.
  typedef Tpetra::Vector<>::mag_type mag_type;
  //typedef Teuchos::ScalarTraits<scalar_type>::magnitudeType mag_type;
  const mag_type theNorm = y.norm2 ();

  // Print the norm of y on Proc 0.
  if (myRank == 0) {
    out << "Norm of y: " << theNorm << endl;
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
    // Never allow Tpetra objects to persist past ScopeGuard's
    // destructor.
    auto comm = Tpetra::getDefaultComm ();
    exampleRoutine (comm, std::cout);

    // This tells the Trilinos test framework that the test passed.
    if (comm->getRank () == 0) {
      std::cout << "End Result: TEST PASSED" << std::endl;
    }
  }
  return 0;
}
