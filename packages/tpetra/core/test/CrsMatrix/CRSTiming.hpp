/*
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
*/

#ifndef CRSTIMING_HPP_
#define CRSTIMING_HPP_

#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Vector.hpp"
#include <Teuchos_TimeMonitor.hpp>

template <class Node>
void CRSTiming(const Teuchos::RCP<const Tpetra::CrsMatrix<double,int,int,Node> > &A) {
  typedef Tpetra::Map<int,int,Node>                Map;
  typedef Tpetra::Vector<double,int,int,Node>      Vector;

  Teuchos::RCP<const Map> map = A->getRowMap();
  Teuchos::RCP<const Teuchos::Comm<int> > comm = map->getComm();

  const int NUM_ITERS = 10;
  const bool speak = (comm->getRank() == 0);

  if (speak) {
    std::cout << "CRSTiming<" << Teuchos::TypeNameTraits<Node>::name() << ">" << std::endl;
  }

  Teuchos::RCP<Vector> x = Tpetra::createVector<double>(map),
                       y = Tpetra::createVector<double>(map);
  Teuchos::RCP<Teuchos::Time> timeMatVec = Teuchos::TimeMonitor::getNewTimer("CRS Mat-Vec");
  {
    Teuchos::TimeMonitor lcltimer(*timeMatVec);
    for (int i=0; i<NUM_ITERS; ++i) {
      A->apply(*x, *y);
    }
  }
  if (speak) {
    std::cout << timeMatVec->name() << ": " << timeMatVec->totalElapsedTime()/(double)(NUM_ITERS) << std::endl;
  }
}

//   TEUCHOS_UNIT_TEST_TEMPLATE_2_DECL( CrsTiming, PowerTriDiag, GRPH, MAT )
//   {
//     typedef typename MAT::ScalarType        Scalar;
//     typedef typename MAT::OrdinalType       Ordinal;
//     typedef typename MAT::NodeType          Node;
//     typedef MultiVector<Scalar,Node>                                 MV;
//     typedef Teuchos::ScalarTraits<Scalar>                            ST;
//
//
//     out << "Testing " << Teuchos::TypeNameTraits<GRPH>::name() << std::endl
//         << "    and " << Teuchos::TypeNameTraits<MAT>::name()  << std::endl;
//
//     RCP<Node> node = getNode<Node>();
//
//     GRPH G(Test::numRows,node);
//     MAT  A(G);
//     //
//     // allocate buffers and fill the graph and matrix
//     //
//     {
//       const size_t totalNNZ = 3*Test::numRows - 2;
//       ArrayRCP<size_t> offsets(Test::numRows+1);
//       ArrayRCP<Ordinal>   inds(totalNNZ);
//       ArrayRCP<Scalar>    vals(totalNNZ);
//       size_t NNZsofar = 0;
//       offsets[0] = NNZsofar;
//       inds[NNZsofar] = 0; inds[NNZsofar+1] =  1;
//       vals[NNZsofar] = 2; vals[NNZsofar+1] = -1;
//       NNZsofar += 2;
//       for (int i=1; i != Test::numRows-1; ++i) {
//         offsets[i] = NNZsofar;
//         inds[NNZsofar] = i-1; inds[NNZsofar+1] = i; inds[NNZsofar+2] = i+1;
//         vals[NNZsofar] =  -1; vals[NNZsofar+1] = 3; vals[NNZsofar+2] =  -1;
//         NNZsofar += 3;
//       }
//       offsets[Test::numRows-1] = NNZsofar;
//       inds[NNZsofar] = Test::numRows-2; inds[NNZsofar+1] = Test::numRows-1;
//       vals[NNZsofar] =  -1;           vals[NNZsofar+1] = 2;
//       NNZsofar += 2;
//       offsets[Test::numRows]   = NNZsofar;
//       TEUCHOS_TEST_FOR_EXCEPT(NNZsofar != totalNNZ);
//       G.set1DStructure(inds, offsets, offsets.persistingView(1,Test::numRows));
//       offsets = Teuchos::null;
//       inds    = Teuchos::null;
//       A.set1DValues(vals);
//       vals    = Teuchos::null;
//       A.finalize(true);
//     }
//     //
//     // fill the matvec
//     //
//     typename SparseOps::template rebind<Scalar>::other matvec(node);
//     matvec.initializeStructure(G);
//     matvec.initializeValues(A);
//     //
//     // time the matvec
//     //
//     MV X(node), Y(node);
//     X.initializeValues( Test::numRows,1, node->template allocBuffer<Scalar>(Test::numRows), Test::numRows);
//     Y.initializeValues( Test::numRows,1, node->template allocBuffer<Scalar>(Test::numRows), Test::numRows);
//     DefaultArithmetic<MV>::Init( X, ST::one() );
//     DefaultArithmetic<MV>::Init( Y, ST::one() );
//     Teuchos::RCP<Teuchos::Time> matvectime = Teuchos::TimeMonitor::getNewTimer("LocalTimer");
//     {
//       Teuchos::TimeMonitor lcltimer(*matvectime);
//       for (int i=0; i<Test::numIters; ++i) {
//         if ( (i%2) == 0 ) {
//           matvec.multiply(Teuchos::NO_TRANS,ST::one(),X,Y);
//         }
//         else {
//           matvec.multiply(Teuchos::NO_TRANS,ST::one(),Y,X);
//         }
//       }
//     }
//     out << "Time is " << matvectime->totalElapsedTime() << std::endl;
//   }

#endif // CRSTIMING_HPP_
