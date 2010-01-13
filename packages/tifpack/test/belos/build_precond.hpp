#ifndef _build_precond_hpp_
#define _build_precond_hpp_

#include <iostream>

#include "Tpetra_FlipOp.hpp"
#include "Tifpack_Factory.hpp"

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
build_precond(Teuchos::ParameterList& test_params,
              const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A)
{
  Teuchos::Time timer("precond");

  typedef Tifpack::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tprec;
  Teuchos::RCP<Tprec> prec;
  Tifpack::Factory<Scalar,LocalOrdinal,GlobalOrdinal,Node> factory;

  std::string prec_name("not specified");
  Tifpack::GetParameter(test_params, "Tifpack::Preconditioner", prec_name);
  prec = factory.Create(prec_name, A);

  Teuchos::ParameterList tif_params;
  if (test_params.isSublist("Tifpack")) {
    tif_params = test_params.sublist("Tifpack");
  }

  if (A->getRowMap()->getComm()->getRank() == 0) {
    std::cout << "Configuring/Initializing/Computing Tifpack preconditioner..."
       << std::endl;
  }

  prec->SetParameters(tif_params);
  prec->Initialize();
  timer.start();
  prec->Compute();
  timer.stop();

  if (A->getRowMap()->getComm()->getRank() == 0) {
    std::cout << "... Finished Computing Tifpack preconditioner (time: "<<timer.totalElapsedTime() << "s)"
       << std::endl;
  }

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  magnitudeType condest = prec->Condest(Tifpack::Cheap);
  if (A->getRowMap()->getComm()->getRank() == 0) {
    std::cout << "Condition estimate(cheap) for preconditioner on proc 0: "
        << condest << std::endl;
  }

  Teuchos::RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec_op =
    Teuchos::rcp(new Tpetra::FlipOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>(prec));
  return prec_op;
}

#endif

