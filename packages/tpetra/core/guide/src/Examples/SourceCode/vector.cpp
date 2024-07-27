// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_Vector.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_Core.hpp>

int
main (int argc, char *argv[])
{

  using std::endl;
  using Teuchos::rcp;
  using Teuchos::RCP;
  using Teuchos::Array;

  // Since Tpetra::Vector takes four template parameters, its type is
  // long.  Here, we use all default values for the template
  // parameters, which helps with the length.
  using vector_type = Tpetra::Vector<>;
  using scalar_type = typename vector_type::scalar_type;
  using global_ordinal_type = typename vector_type::global_ordinal_type;
  using map_type = Tpetra::Map<>;

  Teuchos::oblackholestream blackHole;
  Tpetra::ScopeGuard tpetraScope(&argc, &argv);
  {
    // Get a communicator corresponding to MPI_COMM_WORLD
    Teuchos::RCP<const Teuchos::Comm<int>> comm = Tpetra::getDefaultComm();

    const int myRank = comm->getRank ();
    std::ostream& out = (myRank == 0) ? std::cout : blackHole;

    // Create some Tpetra Map objects
    const size_t numLocalEntries = 5;
    const Tpetra::global_size_t numGlobalEntries =
      comm->getSize () * numLocalEntries;
    const global_ordinal_type indexBase = 0;

    RCP<const map_type> contigMap =
      rcp (new map_type (numGlobalEntries, indexBase, comm));
    TEUCHOS_TEST_FOR_EXCEPTION(
      ! contigMap->isContiguous (), std::logic_error,
      "The supposedly contiguous Map isn't contiguous.");

    // Create a Map which has the same number of global entries per
    // process as contigMap, but distributes them differently, in
    // round-robin (1-D cyclic) fashion instead of contiguously.
    RCP<const map_type> cyclicMap;
    {
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

    // Create some vectors
    //
    // Create a Vector with the contiguous Map.  This version of the
    // constructor will fill in the vector with zeros.
    vector_type x(contigMap);

    // The two-argument copy constructor with second argument
    // Teuchos::Copy performs a deep copy.  x and y have the same Map.
    // The one-argument copy constructor does a _shallow_ copy.
    vector_type y(x, Teuchos::Copy);

    // Create a Vector with the 1-D cyclic Map.  Calling the constructor
    // with false for the second argument leaves the data uninitialized,
    // so that you can fill it later without paying the cost of
    // initially filling it with zeros.
    vector_type z(cyclicMap, false);

    // Set the entries of z to (pseudo)random numbers.  Please don't
    // consider this a good parallel pseudorandom number generator.
    z.randomize();

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
    typedef typename vector_type::mag_type mag_type;
    //typedef Teuchos::ScalarTraits<scalar_type>::magnitudeType mag_type;
    const mag_type theNorm = y.norm2 ();

    // Print the norm of y on Proc 0.
    out << "Norm of y: " << theNorm << endl;

    // This tells the Trilinos test framework that the test passed.
    if (myRank == 0) {
        std::cout << "End Result: TEST PASSED" << std::endl;
    }
    // ScopeGuard's destructor calls MPI_Finalize, if its constructor
    // called MPI_Init.  Likewise, it calls Kokkos::finalize, if its
    // constructor called Kokkos::initialize.
  }
  return 0;
}
