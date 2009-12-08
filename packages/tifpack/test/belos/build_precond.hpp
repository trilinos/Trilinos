#ifndef _build_precond_hpp_
#define _build_precond_hpp_

#include "Tpetra_FlipOp.hpp"
#include "Tifpack_Factory.hpp"

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node>
Teuchos::RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
build_precond(Teuchos::ParameterList& test_params,
              const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> >& A)
{
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

  prec->SetParameters(tif_params);
  prec->Initialize();
  prec->Compute();

  Teuchos::RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> > prec_op =
    Teuchos::rcp(new Tpetra::FlipOp<Scalar,LocalOrdinal,GlobalOrdinal,Node>(prec));
  return prec_op;
}

#endif

