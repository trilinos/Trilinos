#ifndef _build_precond_hpp_
#define _build_precond_hpp_

#include <iostream>

#include "Tifpack_Factory.hpp"

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatVec,class LocalMatSolve>
Teuchos::RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
build_precond(Teuchos::ParameterList& test_params,
              const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatVec,LocalMatSolve> >& A)
{
  Teuchos::Time timer("precond");

  typedef Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tprec;
  Teuchos::RCP<Tprec> prec;
  Tifpack::Factory factory;

  std::string prec_name("not specified");
  Tifpack::GetParameter(test_params, "Tifpack::Preconditioner", prec_name);
  prec = factory.create(prec_name, A);

  Teuchos::ParameterList tif_params;
  if (test_params.isSublist("Tifpack")) {
    tif_params = test_params.sublist("Tifpack");
  }

  if (A->getRowMap()->getComm()->getRank() == 0) {
    std::cout << "Configuring/Initializing/Computing Tifpack preconditioner..."
       << std::endl;
  }

  prec->setParameters(tif_params);
  prec->initialize();
  timer.start();
  prec->compute();
  timer.stop();

  if (A->getRowMap()->getComm()->getRank() == 0) {
    std::cout << "... Finished Computing Tifpack preconditioner (time: "<<timer.totalElapsedTime() << "s)"
       << std::endl;
  }

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  magnitudeType condest = prec->computeCondEst(Tifpack::Cheap);
  if (A->getRowMap()->getComm()->getRank() == 0) {
    std::cout << "Condition estimate(cheap) for preconditioner on proc 0: "
              << condest << std::endl;
  }
  return prec;
}

#endif

