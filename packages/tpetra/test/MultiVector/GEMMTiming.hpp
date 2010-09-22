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
    for (int i=0; i<3; ++i) {
      Z->multiply(Teuchos::CONJ_TRANS, Teuchos::NO_TRANS, 1.0, *X, *Y, 0.0);
    }
  }
  { 
    Teuchos::TimeMonitor lcltimer(*timeLocal);
    for (int i=0; i<3; ++i) {
      Y->multiply(Teuchos::NO_TRANS, Teuchos::NO_TRANS, 1.0, *X, *Z, 0.0);
    }
  }
  if (speak) {
    std::cout << timeReduce->name() << ": " << timeReduce->totalElapsedTime()/3.0 << std::endl;
    std::cout << timeLocal->name() << ": " << timeLocal->totalElapsedTime()/3.0 << std::endl;
  }
}

#endif
