//@HEADER
// ************************************************************************
// 
//               Tpetra: Linear Algebra Services Package 
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ************************************************************************
//@HEADER

#include <Teuchos_Array.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"


// prototype
template <class Scalar, class Ordinal>
Scalar power_method(const Teuchos::RCP<const Tpetra::Operator<Scalar,Ordinal> > &A, size_t niters, typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance, bool verbose);


int main(int argc, char *argv[]) {
  Teuchos::oblackholestream blackhole;
  Teuchos::GlobalMPISession mpiSession(&argc,&argv,&blackhole);
  typedef double Scalar;
  typedef Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  typedef int Ordinal;
  using Tpetra::global_size_t;

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::DefaultPlatform::getDefaultPlatform().getComm();

  size_t myRank = comm->getRank();
  size_t numProc = comm->getSize();
  bool verbose = (myRank==0);

  if (verbose) {
    std::cout << Tpetra::version() << std::endl << std::endl;
  }
  std::cout << *comm;

  // Get the number of local equations from the command line
  if (argc != 2) {
    if (verbose) {
      std::cout << "Usage: " << argv[0] << " number_of_equations" << std::endl;
    }
    std::exit(1);
  }
  const global_size_t numGlobalElements = std::atoi(argv[1]);

  if (numGlobalElements < numProc) {
    if (verbose) {
      std::cout << "numGlobalBlocks = " << numGlobalElements 
                << " cannot be less than the number of processors = " << numProc << std::endl;
    }
    std::exit(1);
  }

  // Construct a Map that puts approximately the same number of 
  // equations on each processor.

  const Ordinal indexBase = 0;
  Teuchos::RCP<const Tpetra::Map<Ordinal> > map; 
  map = Teuchos::rcp( new Tpetra::Map<Ordinal>(numGlobalElements, indexBase, comm) );

  // Get update list and number of local equations from newly created map.

  const size_t numMyElements = map->getNodeNumElements();

  Teuchos::ArrayView<const Ordinal> myGlobalElements = map->getNodeElementList();

  // Create an OTeger vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor

  Teuchos::ArrayRCP<size_t> NumNz = Teuchos::arcp<size_t>(numMyElements);

  // We are building a tridiagonal matrix where each row has (-1 2 -1)
  // So we need 2 off-diagonal terms (except for the first and last equation)

  for (size_t i=0; i < numMyElements; ++i) {
    if (myGlobalElements[i] == 0 || static_cast<global_size_t>(myGlobalElements[i]) == numGlobalElements-1) {
      // boundary
      NumNz[i] = 2;
    }
    else {
      NumNz[i] = 3;
    }
  }

  // Create a Tpetra::Matrix using the Map, with a static allocation dictated by NumNz
  Teuchos::RCP< Tpetra::CrsMatrix<Scalar,Ordinal> > A;
  A = Teuchos::rcp( new Tpetra::CrsMatrix<Scalar,Ordinal>(map, NumNz, Tpetra::StaticProfile) );
  
  // We are done with NumNZ
  NumNz = Teuchos::null;

  // Add  rows one-at-a-time
  // Off diagonal values will always be -1
  const Scalar two    = static_cast<Scalar>( 2.0);
  const Scalar negOne = static_cast<Scalar>(-1.0);
  for (size_t i=0; i<numMyElements; i++) {
    if (myGlobalElements[i] == 0) {
      A->insertGlobalValues( myGlobalElements[i],
                             Teuchos::tuple<Ordinal>( myGlobalElements[i], myGlobalElements[i]+1 ),
                             Teuchos::tuple<Scalar> ( two, negOne ) );
    }
    else if (static_cast<global_size_t>(myGlobalElements[i]) == numGlobalElements-1) {
      A->insertGlobalValues( myGlobalElements[i],
                             Teuchos::tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i] ),
                             Teuchos::tuple<Scalar> ( negOne, two ) );
    }
    else {
      A->insertGlobalValues( myGlobalElements[i],
                             Teuchos::tuple<Ordinal>( myGlobalElements[i]-1, myGlobalElements[i], myGlobalElements[i]+1 ),
                             Teuchos::tuple<Scalar> ( negOne, two, negOne ) );
    }
  }

  // Finish up
  A->fillComplete(Tpetra::DoOptimizeStorage);
  if (verbose) std::cout << std::endl << A->description() << std::endl << std::endl;

  // Create vectors for Power method

  // variable needed for iteration
  Scalar lambda;
  const size_t niters = static_cast<size_t>(numGlobalElements*10);
  const Scalar tolerance = 1.0e-2;

  // Iterate
  lambda = power_method<Scalar,Ordinal>(A, niters, tolerance, verbose);

  // Increase diagonal dominance
  if (verbose) {
    std::cout << "\nIncreasing magnitude of first diagonal term, solving again\n"
              << std::endl;
  }

  if (A->getRowMap()->isNodeGlobalElement(0)) {
    // get a copy of the row with with global index 0
    // modify the diagonal entry of that row
    // submit the modified values to the matrix
    const Ordinal ID = 0;
    size_t numVals = A->getNumEntriesInGlobalRow(ID);
    Teuchos::Array<Scalar>  rowvals(numVals);
    Teuchos::Array<Ordinal> rowinds(numVals);
    A->getGlobalRowCopy(ID, rowinds, rowvals, numVals);       // Get A(0,:)
    for (size_t i=0; i<numVals; i++) {
      if (rowinds[i] == ID) {
        // we have found the diagonal; modify it and break the loop
        rowvals[i] *= 10.0;
        break;
      }
    }
    A->replaceGlobalValues(ID, rowinds(), rowvals());
  }

  // Iterate (again)
  lambda = power_method<Scalar,Ordinal>(A, niters, tolerance, verbose);  

  /* end main
   */
  if (verbose) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  return 0;
}


template <class Scalar, class Ordinal>
Scalar power_method(const Teuchos::RCP<const Tpetra::Operator<Scalar,Ordinal> > &A, size_t niters, typename Teuchos::ScalarTraits<Scalar>::magnitudeType tolerance, bool verbose) {
  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
  // create three vectors; do not bother initializing them to zero
  Tpetra::Vector<Scalar,Ordinal> q(A->getRangeMap(), false);
  Tpetra::Vector<Scalar,Ordinal> z(A->getRangeMap(), false);
  Tpetra::Vector<Scalar,Ordinal> resid(A->getRangeMap(), false);

  // Fill z with random numbers
  z.randomize();

  // Variables needed for iteration
  Scalar lambda = static_cast<Scalar>(0.0);
  Magnitude normz, residual = static_cast<Magnitude>(0.0);

  const Scalar one  = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();

  for (size_t iter = 0; iter < niters; ++iter) {
    normz = z.norm2();                            // Compute 2-norm of z
    q.scale(one/normz, z);                        // Set q = z / normz
    A->apply(q, z);                               // Compute z = A*q
    lambda = q.dot(z);                            // Approximate maximum eigenvalue: lamba = dot(q,z)
    if ( iter % 100 == 0 || iter + 1 == niters ) {
      resid.update(one, z, -lambda, q, zero);     // Compute A*q - lambda*q
      residual = resid.norm2();
      if (verbose) {
        std::cout << "Iter = " << iter << "  Lambda = " << lambda 
                  << "  Residual of A*q - lambda*q = " 
                  << residual << std::endl;
      }
    } 
    if (residual < tolerance) {
      break;
    }
  }
  return lambda;
}
