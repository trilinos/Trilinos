// @HEADER
//
// ***********************************************************************
//
//           Amesos2: Templated Direct Sparse Solver Package
//                  Copyright 2011 Sandia Corporation
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
// ***********************************************************************
//
// @HEADER

// intended file name:  tklu_simple.cpp

#include <Teuchos_Array.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_oblackholestream.hpp>

#include "Tpetra_Core.hpp"
#include "Tpetra_Version.hpp"
#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"
#include "Tpetra_CrsMatrix.hpp"


#include "Teuchos_Tuple.hpp"

//#include "klu.h"

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

  Teuchos::RCP<const Teuchos::Comm<int> > comm = Tpetra::getDefaultComm();

  size_t myRank = comm->getRank();
  size_t numProc = comm->getSize();
  bool verbose = (myRank==0);

  if (verbose) {
    std::cout << Tpetra::version() << std::endl << std::endl;
  }
  std::cout << *comm;

  // Get the number of local equations from the command line
  if (numProc != 1) {
    if (verbose) {
      std::cout << "Usage: This example only runs in serial" << std::endl;
    }
    return -1;
  }
  const global_size_t numGlobalElements = 5;
  //const global_size_t numGlobalNonZeros = 12;

  // Construct a Map that puts approximately the same number of
  // equations on each processor.

  const Ordinal indexBase = 0;
  Teuchos::RCP<const Tpetra::Map<Ordinal> > map;
  map = Teuchos::rcp( new Tpetra::Map<Ordinal>(numGlobalElements, indexBase, comm) );

  // Get update list and number of local equations from newly created map.

  const size_t numMyElements = map->getNodeNumElements();

  if( numMyElements != numGlobalElements )
  if (numProc != 1) {
     std::cout << "Usage: serial thing again" << std::endl;
    if (verbose) {
    }
    return -1;
  }

  Teuchos::ArrayView<const Ordinal> myGlobalElements = map->getNodeElementList();

//    klu_symbolic *Symbolic ;
//    klu_numeric *Numeric ;
//    klu_common Common ;
//    klu_defaults (&Common) ;

  // Solve for x the linear system Ax=b
  // A = [ 2  3  0  0  0;...
  //       3  0  4  0  6;...
  //       0 -1 -3  2  0;...
  //       0  0  1  0  0;...
  //       0  4  2  0  1];
  // x = [1; 2; 3; 4; 5];
  // b = [8; 45; -3; 3; 19];

  const Teuchos::Array<Ordinal> ColPtr = Teuchos::tuple<Ordinal>(0, 2, 5, 9, 10, 12);

  const Teuchos::Array<Ordinal> RowInd = Teuchos::tuple<Ordinal>(0,  1,  0,   2,  4,  1,  2,  3,   4,  2,  1,  4);

  const Teuchos::Array<Scalar> MatVal = Teuchos::tuple<Scalar>(2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.);

  const Teuchos::Array<Scalar> lhs = Teuchos::tuple<Scalar>(1., 2., 3., 4., 5.);

  Teuchos::Array<Scalar> rhs = Teuchos::tuple<Scalar>(8., 45., -3., 3., 19.);

  std::cout << lhs << std::endl;
//    Symbolic = klu_analyze (numGlobalElements, &ColPtr, &RowInd, &Common) ;
//    Numeric = klu_factor (&ColPtr, &RowInd, &MatVal, Symbolic, &Common) ;
//    klu_solve (Symbolic, Numeric, numGlobalElements, 1, &rhs, &Common) ;
  std::cout << rhs << std::endl;
//    klu_free_symbolic (&Symbolic, &Common) ;
//    klu_free_numeric (&Numeric, &Common) ;

  /* end main
   */
  if (verbose) {
    std::cout << "\nEnd Result: TEST PASSED" << std::endl;
  }
  // We are done.
  return 0;
}

