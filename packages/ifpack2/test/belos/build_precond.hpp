#ifndef _build_precond_hpp_
#define _build_precond_hpp_

#include <iostream>

#include "Ifpack2_Factory.hpp"

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class LocalMatOps>
Teuchos::RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
build_precond(Teuchos::ParameterList& test_params,
              const Teuchos::RCP<const Tpetra::CrsMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node,LocalMatOps> >& A)
{
  using Teuchos::FancyOStream;
  using Teuchos::getFancyOStream;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  using std::cout;
  using std::endl;
  Teuchos::Time timer("precond");
  const int myRank = A->getRowMap ()->getComm ()->getRank ();

  RCP<FancyOStream> out = getFancyOStream (rcpFromRef (cout));

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

  if (myRank == 0) {
    *out << "Configuring, initializing, and computing Ifpack2 preconditioner" << endl;
  }
  {
    OSTab tab (*out);
    prec->setParameters (tif_params);
    prec->initialize ();
    {
      Teuchos::TimeMonitor timeMon (timer);
      prec->compute ();
    }
    if (myRank == 0) {
      *out << "Finished computing Ifpack2 preconditioner" << endl;
      OSTab tab2 (*out);
      *out << "Time (s): " << timer.totalElapsedTime () << endl;
    }
  }
  if (myRank == 0) {
    *out << "Preconditioner attributes:" << endl;
    OSTab tab (*out);
    prec->describe (*out, Teuchos::VERB_LOW);
  }

  typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitudeType;
  const magnitudeType condest = prec->computeCondEst (Ifpack2::Cheap);
  if (myRank == 0) {
    *out << "Condition number estimate (cheap): " << condest << endl;
  }
  return prec;
}

#endif

