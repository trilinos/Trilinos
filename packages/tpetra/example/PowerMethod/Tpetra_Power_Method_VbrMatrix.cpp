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
#include "Tpetra_BlockMap.hpp"
#include "Tpetra_BlockMultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_VbrMatrix.hpp"


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

  // Get the number of global equations from the command line
  if (argc != 2) {
    if (verbose) {
      std::cout << "Usage: " << argv[0] << " number_of_equations" << std::endl;
    }
    return -1;
  }
  const int blockSize = 3;
  global_size_t numGlobalElements = std::atoi(argv[1]);

  //make sure numGlobalElements is an integer multiple of blockSize.
  numGlobalElements += blockSize - numGlobalElements%blockSize;
  global_size_t numGlobalBlocks = numGlobalElements/blockSize;

  if (numGlobalBlocks < numProc) {
    if (verbose) {
      std::cout << "numGlobalBlocks = " << numGlobalBlocks 
                << " cannot be less than the number of processors = " << numProc << std::endl;
    }
    return -1;
  }

  if (verbose) {
    std::cout << "numGlobalBlocks = " << numGlobalBlocks << ", numGlobalElements (point-rows): " << numGlobalElements << std::endl;
  }

  // Construct a Map that puts approximately the same number of equations on each processor.

  Ordinal indexBase = 0;
  Teuchos::RCP<const Tpetra::BlockMap<Ordinal> > blkmap = Teuchos::rcp(new Tpetra::BlockMap<Ordinal>(numGlobalBlocks, blockSize, indexBase, comm));

  // Get update list and number of local equations from newly created map.

  const size_t numMyBlocks = blkmap->getNodeNumBlocks();

  Teuchos::ArrayView<const Ordinal> myGlobalBlocks = blkmap->getNodeBlockIDs();

  // Create an OTeger vector NumNz that is used to build the Petra Matrix.
  // NumNz[i] is the Number of OFF-DIAGONAL term for the ith global equation 
  // on this processor

  // Create a Tpetra::VbrMatrix using the BlockMap.
  Teuchos::RCP< Tpetra::VbrMatrix<Scalar,Ordinal> > A;
  A = Teuchos::rcp( new Tpetra::VbrMatrix<Scalar,Ordinal>(blkmap, 3) );
  
  // Add  rows one-at-a-time.
  // We will build a block-tridiagonal matrix where the diagonal block has 2
  // on the diagonal and the off-diagonal blocks have -1 on the diagonal.
  Teuchos::Array<Scalar> two(blockSize*blockSize, 0);
  two[0] = 2.0; two[blockSize+1] = 2.0; two[blockSize*2+2] = 2.0;
  Teuchos::Array<Scalar> negOne(blockSize*blockSize, 0);
  negOne[0] = -1.0; negOne[blockSize+1] = -1.0; negOne[blockSize*2+2] = -1.0;

  for (size_t i=0; i<numMyBlocks; i++) {
    if (myGlobalBlocks[i] != 0) {
      A->setGlobalBlockEntry( myGlobalBlocks[i], myGlobalBlocks[i]-1,
                              blockSize, blockSize, blockSize, negOne() );
    }

    //always set the diagonal block:
    A->setGlobalBlockEntry( myGlobalBlocks[i], myGlobalBlocks[i],
                            blockSize, blockSize, blockSize, two() );

    if (static_cast<global_size_t>(myGlobalBlocks[i]) != numGlobalBlocks-1) {
      A->setGlobalBlockEntry( myGlobalBlocks[i], myGlobalBlocks[i]+1,
                              blockSize, blockSize, blockSize, negOne() );
    }
  }

  // Finish up matrix fill
  A->fillComplete();
  if (verbose) std::cout << std::endl << A->description() << std::endl << std::endl;

  // variable needed for iteration
  Scalar lambda; (void)lambda;
  const size_t niters = static_cast<size_t>(numGlobalElements*10);
  const Scalar tolerance = 1.0e-2;

  // Iterate
  lambda = power_method<Scalar,Ordinal>(A, niters, tolerance, verbose);

  // Increase diagonal dominance
  if (verbose) {
    std::cout << "\nIncreasing magnitude of first diagonal term, solving again\n"
              << std::endl;
  }

  if (A->getBlockRowMap()->getPointMap()->isNodeGlobalElement(0)) {
    // modify the diagonal entry of the row with global-id 0
    const Ordinal ID = 0;
    Ordinal numPtRows, numPtCols;
    Teuchos::ArrayRCP<Scalar> blockEntry;
    A->getGlobalBlockEntryViewNonConst(ID,ID, numPtRows,numPtCols, blockEntry);
    blockEntry[0] *= 10.0;
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
  const bool NO_INITIALIZE_TO_ZERO = false;
  // create three vectors; do not bother initializing q to zero, as we will fill it with random below
  Tpetra::Vector<Scalar,Ordinal> z(A->getRangeMap(), NO_INITIALIZE_TO_ZERO),
                                 q(A->getRangeMap(), NO_INITIALIZE_TO_ZERO),
                                 r(A->getRangeMap(), NO_INITIALIZE_TO_ZERO);
  // Fill z with random numbers
  z.randomize();
  // Variables needed for iteration
  const Scalar ONE  = Teuchos::ScalarTraits<Scalar>::one();
  const Scalar ZERO = Teuchos::ScalarTraits<Scalar>::zero();
  Scalar lambda = static_cast<Scalar>(0.0);
  Magnitude normz, residual = static_cast<Magnitude>(0.0);
  // power iteration
  for (size_t iter = 0; iter < niters; ++iter) {
    normz = z.norm2();          // Compute 2-norm of z
    q.scale(ONE/normz, z);      // Set q = z / normz
    A->apply(q, z);             // Compute z = A*q
    lambda = q.dot(z);          // Approximate maximum eigenvalue: lamba = dot(q,z)
    if ( iter % 100 == 0 || iter + 1 == niters ) {
      r.update(ONE, z, -lambda, q, ZERO);     // Compute A*q - lambda*q
      residual = Teuchos::ScalarTraits<Scalar>::magnitude(r.norm2() / lambda);
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
