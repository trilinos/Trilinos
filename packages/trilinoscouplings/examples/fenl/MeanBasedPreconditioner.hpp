// @HEADER
// *****************************************************************************
//           Trilinos: An Object-Oriented Solver Framework
//
// Copyright 2001-2024 NTESS and the Trilinos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STOKHOS_MEAN_BASED_PRECONDITIONER_HPP
#define STOKHOS_MEAN_BASED_PRECONDITIONER_HPP

#include "MueLuPreconditioner.hpp"

namespace Kokkos {
namespace Example {

  // Default implementation just uses the original matrix
  template<class Scalar, class LO, class GO, class N>
  Teuchos::RCP<Tpetra::Operator<Scalar,LO,GO,N> >
  build_mean_based_muelu_preconditioner(
    const Teuchos::RCP<Tpetra::CrsMatrix<Scalar,LO,GO,N> >& A,
    const Teuchos::RCP<Teuchos::ParameterList>& precParams,
    const Teuchos::RCP<Tpetra::MultiVector<double,LO,GO,N> >& coords)
  {
    return build_muelu_preconditioner(A, precParams, coords);
  }

}
}

#endif // STOKHOS_MEAN_BASED_PRECONDITIONER_HPP
