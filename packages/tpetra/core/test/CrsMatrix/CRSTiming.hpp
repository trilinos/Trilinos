// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

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

#endif // CRSTIMING_HPP_
