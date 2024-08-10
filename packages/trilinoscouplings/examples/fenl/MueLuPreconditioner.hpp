// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_MUELU_PRECONDITIONER_HPP
#define STOKHOS_MUELU_PRECONDITIONER_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Tpetra_MultiVector.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"

namespace Kokkos {
namespace Example {

  template<class S, class LO, class GO, class N>
  Teuchos::RCP<Tpetra::Operator<S,LO,GO,N> >
  build_muelu_preconditioner(const Teuchos::RCP<Tpetra::CrsMatrix<S,LO,GO,N> >& A,
                             const Teuchos::RCP<Teuchos::ParameterList>& mueluParams,
                             const Teuchos::RCP<Tpetra::MultiVector<double,LO,GO,N> >& coords)
  {
    using Teuchos::RCP;
    using OperatorType       = Tpetra::Operator<S,LO,GO,N>;
    using PreconditionerType = MueLu::TpetraOperator<S,LO,GO,N>;

    RCP<OperatorType> A_op = A;

    Teuchos::ParameterList& userList = mueluParams->sublist("user data");
    if (coords != Teuchos::null) {
      userList.set<RCP<Tpetra::MultiVector<double,LO,GO,N> > >("Coordinates", coords);
    }
    RCP<PreconditionerType> mueluPreconditioner
      = MueLu::CreateTpetraPreconditioner(A_op, *mueluParams);

    return mueluPreconditioner;
  }

}
}

#endif // STOKHOS_MUELU_PRECONDITIONER_HPP
