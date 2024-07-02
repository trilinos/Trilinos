// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef _build_precond_hpp_
#define _build_precond_hpp_

#include <iostream>

#include "Ifpack2_Factory.hpp"

template<class Scalar,class LocalOrdinal,class GlobalOrdinal,class Node,class MatrixType>
Teuchos::RCP<Tpetra::Operator<Scalar,LocalOrdinal,GlobalOrdinal,Node> >
build_precond (Teuchos::ParameterList& test_params,
               const Teuchos::RCP<const MatrixType>& A)
{
  using Teuchos::FancyOStream;
  using Teuchos::getFancyOStream;
  using Teuchos::OSTab;
  using Teuchos::RCP;
  using Teuchos::rcpFromRef;
  using std::cout;
  using std::endl;
  typedef Tpetra::RowMatrix<Scalar,LocalOrdinal,GlobalOrdinal,Node> row_matrix_type;
  Teuchos::Time timer_init("Init preconditioner");
  Teuchos::Time timer("Compute preconditioner");
  Teuchos::Time timer2("Compute preconditioner (reuse)");
  auto comm = A->getRowMap ()->getComm ();
  const int myRank = comm->getRank ();

  RCP<FancyOStream> out = getFancyOStream (rcpFromRef (cout));

  typedef Ifpack2::Preconditioner<Scalar,LocalOrdinal,GlobalOrdinal,Node> Tprec;
  Teuchos::RCP<Tprec> prec;
  Ifpack2::Factory factory;

  std::string prec_name("not specified");
  Ifpack2::getParameter(test_params, "Ifpack2::Preconditioner", prec_name);
  prec = factory.create<row_matrix_type> (prec_name, A);

  Teuchos::ParameterList tif_params;
  if (test_params.isSublist("Ifpack2")) {
    tif_params = test_params.sublist("Ifpack2");
  }
  bool reuse_pattern = false;
  if (test_params.isParameter("Reuse Pattern"))
    {
      reuse_pattern = test_params.get<bool>("Reuse Pattern");
    }

  if (myRank == 0) {
    *out << "Configuring, initializing, and computing Ifpack2 preconditioner" << endl;
  }
  comm->barrier();
  {
    OSTab tab (*out);
    prec->setParameters (tif_params);
    {
      Teuchos::TimeMonitor timeMon (timer_init);
      prec->initialize ();
      comm->barrier();
    }
     if (myRank == 0) {
      *out << "Time Init: " << timer_init.totalElapsedTime() << endl;
     }
    {
      Teuchos::TimeMonitor timeMon (timer);
      prec->compute ();
      comm->barrier();
    }
    if (myRank == 0) {
      *out << "Finished computing Ifpack2 preconditioner" << endl;
      OSTab tab2 (*out);
      *out << "Time (s): " << timer.totalElapsedTime () << endl;
    }
    if (reuse_pattern == true)
      {
	{
	  Teuchos::TimeMonitor timeMon (timer2);
	  prec->compute ();
          comm->barrier();
	}
	if (myRank == 0) {
	  *out << "Finished recomputing Ifpack2 preconditioner" 
	       << endl;
	  OSTab tab2 (*out);
	  *out << "Time (s): " << timer2.totalElapsedTime () 
	       << endl;
	}
      }
  }

  return prec;
}

#endif

