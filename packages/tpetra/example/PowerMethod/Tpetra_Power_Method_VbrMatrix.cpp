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
  typedef Scalar scalar_type;
  typedef Teuchos::ScalarTraits<scalar_type> STS;
  typedef typename STS::magnitudeType magnitude_type;
  typedef Teuchos::ScalarTraits<magnitude_type> STM;

  const bool NO_INITIALIZE_TO_ZERO = false;
  // create three vectors; do not bother initializing q to zero, as we will fill it with random below
  Tpetra::Vector<scalar_type, Ordinal> z (A->getRangeMap (), NO_INITIALIZE_TO_ZERO);

  // mfh 26 Nov 2012: For some reason, not initializing r to zero
  // causes the norm of r to be NaN in the loop below.  This fixes
  // a test failure on my workstation.
  Tpetra::Vector<scalar_type, Ordinal> q (A->getRangeMap ()); //, NO_INITIALIZE_TO_ZERO);
  Tpetra::Vector<scalar_type, Ordinal> r (A->getRangeMap ()); //, NO_INITIALIZE_TO_ZERO);

  // Fill z with random numbers
  z.randomize ();
  // Variables needed for iteration
  scalar_type lambda = STS::zero ();
  magnitude_type normz, residual = STM::zero ();

  // power iteration
  for (size_t iter = 0; iter < niters; ++iter) {
    normz = z.norm2 ();               // Compute 2-norm of z
    q.scale (STS::one () / normz, z); // Compute q = z / normz
    A->apply (q, z);                  // Compute z = A*q
    lambda = q.dot (z);               // Approximate maximum eigenvalue: lambda = dot(q, z)
    if ( iter % 100 == 0 || iter + 1 == niters ) {
      r.update (STS::one (), z, -lambda, q, STS::zero ());     // Compute A*q - lambda*q

      magnitude_type r_norm = r.norm2 ();
      residual = STS::magnitude (r_norm / lambda);
      if (verbose) {
        std::cout << "Iter = " << iter << ": lambda = " << lambda
                  << ", ||A*q - lambda*q||_2 = " << r_norm
                  << ", ||A*q - lambda*q||_2 / lambda = " << residual
                  << std::endl;
      }
    }
    if (residual < tolerance) {
      break;
    }
  }
  return lambda;
}
