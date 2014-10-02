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
\example lesson02_init_map_vec.cpp
\brief Create data distributions (Tpetra::Map), and create vectors
  (Tpetra::Vector) distributed according to those distributions.

\ref Tpetra_Lesson02 explains this example in detail.
*/

#include <Tpetra_DefaultPlatform.hpp>
#include <Tpetra_Vector.hpp>
#include <Tpetra_Version.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

void
exampleRoutine (const Teuchos::RCP<const Teuchos::Comm<int> >& comm,
                std::ostream& out)
{
  using std::endl;
  using Teuchos::Array;
  using Teuchos::ArrayView;
  using Teuchos::RCP;
  using Teuchos::rcp;

  // Print out the Tpetra software version information.
  out << Tpetra::version () << endl << endl;

  //
  // The first thing you might notice that makes Tpetra objects
  // different than their Epetra counterparts, is that Tpetra objects
  // take several template parameters.  These template parameters give
  // Tpetra its features of being able to solve very large problems
  // (of more than 2 billion unknowns) and to exploit intranode
  // parallelism.  Most Tpetra objects come with default values of
  // those template parameters.  In many cases, you might not have to
  // specify _any_ of those values explicitly!  You can also control
  // some of their default values, like that of the Node type, when
  // configuring Trilinos.
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
  // parameters' values revert to their defaults.
  typedef Tpetra::Vector<>::scalar_type scalar_type;

  // The "LocalOrdinal" (LO) type is the type of "local" indices.
  // Both Epetra and Tpetra index local entries differently than
  // global entries.  Tpetra exploits this so that you can use a
  // shorter integer type for local indices.  This saves bandwidth
  // when computing sparse matrix-vector products.  We use the default
  // LO type here, which we get from Tpetra::Vector.  We could also
  // get it from Tpetra::Map.
  typedef Tpetra::Vector<>::local_ordinal_type local_ordinal_type;

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
  typedef Tpetra::Vector<>::node_type node_type;

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

  // The total (global, i.e., over all MPI processes) number of
  // entries in the Map.  Tpetra's global_size_t type is an unsigned
  // type and is at least 64 bits long on 64-bit machines.
  //
  // For this example, we scale the global number of entries in the
  // Map with the number of MPI processes.  That way, you can run this
  // example with any number of MPI processes and every process will
  // still have a positive number of entries.
  const Tpetra::global_size_t numGlobalEntries = comm->getSize() * 5;

  // Tpetra can index the entries of a Map starting with 0 (C style),
  // 1 (Fortran style), or any base you want.  1-based indexing is
  // handy when interfacing with Fortran.  We choose 0-based indexing
  // here.
  const global_ordinal_type indexBase = 0;

  // Construct a Map that puts the same number of equations on each
  // processor.  It's typical to create a const Map.  Maps should be
  // considered immutable objects.  If you want a new data
  // distribution, create a new Map.
  RCP<const map_type> contigMap =
    rcp (new map_type (numGlobalEntries, indexBase, comm,
                       Tpetra::GloballyDistributed));

  // contigMap is contiguous by construction.
  TEUCHOS_TEST_FOR_EXCEPTION(! contigMap->isContiguous(), std::logic_error,
    "The supposedly contiguous Map isn't contiguous.");

  // Let's create a second Map.  It will have the same number of
  // global entries per process, but will distribute them
  // differently, in round-robin (1-D cyclic) fashion instead of
  // contiguously.
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

    const int numProcs = comm->getSize();
    const int myRank = comm->getRank();
    for (Array<global_ordinal_type>::size_type k = 0; k < numEltsPerProc; ++k)
      elementList[k] = myRank + k*numProcs;

    cyclicMap = rcp (new map_type (numGlobalEntries, elementList, indexBase, comm));
  }

  // If there's more than one MPI process in the communicator,
  // then cyclicMap is definitely NOT contiguous.
  TEUCHOS_TEST_FOR_EXCEPTION(
    comm->getSize() > 1 && cyclicMap->isContiguous(),
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
  RCP<vector_type> x = rcp (new vector_type (contigMap));

  // The copy constructor performs a deep copy.
  // x and y have the same Map.
  RCP<vector_type> y = rcp (new vector_type (*x));

  // Create a Vector with the 1-D cyclic Map.  Calling the constructor
  // with false for the second argument leaves the data uninitialized,
  // so that you can fill it later without paying the cost of
  // initially filling it with zeros.
  RCP<vector_type> z = rcp (new vector_type (contigMap, false));

  // Set the entries of z to (pseudo)random numbers.  Please don't
  // consider this a good parallel pseudorandom number generator.
  z->randomize ();

  // Set the entries of x to all ones.
  //
  // The code below works because scalar_type=double, and C++
  // defines the implicit conversion from the integer 1 to double.
  // In general, you may use the commented-out line of code, if
  // the conversion from int to scalar_type is not defined for your
  // scalar type.
  x->putScalar (1);
  //x->putScalar (Teuchos::ScalarTraits<scalar_type>::one());

  // See comment above about type conversions to scalar_type.
  const scalar_type alpha = 3.14159;
  const scalar_type beta = 2.71828;
  const scalar_type gamma = -10;

  // x = beta*x + alpha*z
  //
  // This is a legal operation!  Even though the Maps of x and z are
  // not the same, their Maps are compatible.  Whether it makes sense
  // or not depends on your application.
  x->update (alpha, *z, beta);

  // See comment above about type conversions to scalar_type.
  y->putScalar (42);
  // y = gamma*y + alpha*x + beta*z
  y->update (alpha, *x, beta, *z, gamma);

  // Compute the 2-norm of y.
  //
  // The norm may have a different type than scalar_type.
  // For example, if scalar_type is complex, then the norm is real.
  // The ScalarTraits "traits class" gives us the type of the norm.
  typedef Teuchos::ScalarTraits<scalar_type>::magnitudeType magnitude_type;
  const magnitude_type theNorm = y->norm2 ();

  // Print the norm of y on Proc 0.
  out << "Norm of y: " << theNorm << endl;
}

//
// The same main() driver routine as in the first Tpetra lesson.
//
int
main (int argc, char *argv[])
{
  using std::endl;
  using Teuchos::RCP;

  Teuchos::oblackholestream blackHole;
  Teuchos::GlobalMPISession mpiSession (&argc, &argv, &blackHole);
  RCP<const Teuchos::Comm<int> > comm =
    Tpetra::DefaultPlatform::getDefaultPlatform ().getComm ();

  const int myRank = comm->getRank ();
  std::ostream& out = (myRank == 0) ? std::cout : blackHole;

  // We have a communicator and an output stream.
  // Let's do something with them!
  exampleRoutine (comm, out);

  // This tells the Trilinos test framework that the test passed.
  if (myRank == 0) {
    std::cout << "End Result: TEST PASSED" << std::endl;
  }

  return 0;
}
