#ifndef _build_precond_hpp_
#define _build_precond_hpp_

#include <iostream>

#include "Ifpack2_Factory.hpp"

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
Teuchos::RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
build_precond(Teuchos::ParameterList& test_params,
              const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> >& A)
{
  Teuchos::Time timer("precond");

  typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tprec;
  Teuchos::RCP<Tprec> prec;
  Ifpack2::Factory factory;

  std::string prec_name("not specified");
  Ifpack2::getParameter(test_params, "Ifpack2::Preconditioner", prec_name);
  prec = factory.create(prec_name, A);

  Teuchos::ParameterList tif_params;
  if (test_params.isSublist("Ifpack2")) {
    tif_params = test_params.sublist("Ifpack2");
  }

  if (A->getRowMap()->getComm()->getRank() == 0) {
    std::cout << "Configuring/Initializing/Computing Ifpack2 preconditioner..."
       << std::endl;
  }

  prec->setParameters(tif_params);
  prec->initialize();
  timer.start();
  prec->compute();
  timer.stop();

  if (A->getRowMap()->getComm()->getRank() == 0) {
    std::cout << "... Finished Computing Ifpack2 preconditioner (time: "<<timer.totalElapsedTime() << "s)"
       << std::endl;
  }

  // typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  // magnitudeType condest = prec->computeCondEst(Ifpack2::Cheap);
  // if (A->getRowMap()->getComm()->getRank() == 0) {
  //   std::cout << "Condition estimate(cheap) for preconditioner on proc 0: "
  //             << condest << std::endl;
  // }
  return prec;
}

#endif

