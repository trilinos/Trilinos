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

#ifndef GEMMTIMING_HPP_
#define GEMMTIMING_HPP_

#include "Tpetra_Version.hpp"
#include "Tpetra_MultiVector.hpp"
#include <Teuchos_TimeMonitor.hpp>

template <class Node> 
void GEMMTiming(int M, int N, const Teuchos::RCP<const Teuchos::Comm<int> > &comm, const Teuchos::RCP<Node> &node) {
  typedef Tpetra::Map<int,int,Node>                Map;
  typedef Tpetra::MultiVector<double,int,int,Node> MultiVector;
  const bool speak = (comm->getRank() == 0);
  if (speak) {
    std::cout << "GEMMTiming<" << Teuchos::typeName(*node) << ">" << std::endl;
  }

  Teuchos::RCP<const Map> dstmap = Tpetra::createUniformContigMapWithNode<int,int,Node>(M, comm, node);
  Teuchos::RCP<MultiVector> X = Tpetra::createMultiVector<double>(dstmap,N),
                            Y = Tpetra::createMultiVector<double>(dstmap,N);
  Teuchos::RCP<const Map> lclmap = Tpetra::createLocalMapWithNode<int,int,Node>(N, comm, node);
  Teuchos::RCP<MultiVector> Z = Tpetra::createMultiVector<double>(lclmap,N);
  X->putScalar(0.0);
  Y->putScalar(0.0);
  Z->putScalar(0.0);
  Teuchos::RCP<Teuchos::Time> timeReduce = Teuchos::TimeMonitor::getNewTimer("Z_lcl = X^T * Y"),
                              timeLocal  = Teuchos::TimeMonitor::getNewTimer("Y = X * Z_lcl");
  {
    Teuchos::TimeMonitor lcltimer(*timeReduce);
    for (int i=0; i<10; ++i) {
      Z->multiply(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, 1.0, *X, *Y, 0.0);
    }
  }
  { 
    Teuchos::TimeMonitor lcltimer(*timeLocal);
    for (int i=0; i<10; ++i) {
      Y->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *X, *Z, 0.0);
    }
  }
  if (speak) {
    std::cout << timeReduce->name() << ": " << timeReduce->totalElapsedTime()/10.0 << std::endl;
    std::cout << timeLocal->name() << ": " << timeLocal->totalElapsedTime()/10.0 << std::endl;
  }
}

#endif
