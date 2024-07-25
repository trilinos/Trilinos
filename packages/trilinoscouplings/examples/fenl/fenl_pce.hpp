//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER

// Use "scalar" version of mean-based preconditioner (i.e., a preconditioner
// with double as the scalar type).  This is currently necessary to get the
// MueLu tests to pass on OpenMP and Cuda due to various kernels that don't
// work with the PCE scalar type.
#define USE_SCALAR_MEAN_BASED_PREC 1

#include "Stokhos_Tpetra_UQ_PCE.hpp"
#include "Stokhos_Sacado_Kokkos_MP_Vector.hpp"

#if defined( HAVE_STOKHOS_BELOS )
#include "Belos_TpetraAdapter_UQ_PCE.hpp"
#endif

#if defined( HAVE_STOKHOS_MUELU )
#include "Stokhos_MueLu_UQ_PCE.hpp"
#endif

#include <Kokkos_Core.hpp>
#include <HexElement.hpp>
#include <fenl.hpp>
#include <fenl_functors_pce.hpp>

#if defined(HAVE_TRILINOSCOUPLINGS_BELOS) && defined(HAVE_TRILINOSCOUPLINGS_MUELU)
#include <MeanBasedPreconditioner.hpp>
#include "Stokhos_Tpetra_Utilities.hpp"

namespace Kokkos {
namespace Example {

  // Overload of build_mean_based_muelu_preconditioner() for UQ::PCE scalar type
  template<class S, class LO, class GO, class N>
  Teuchos::RCP<Tpetra::Operator<Sacado::UQ::PCE<S>,LO,GO,N> >
  build_mean_based_muelu_preconditioner(
    const Teuchos::RCP<Tpetra::CrsMatrix<Sacado::UQ::PCE<S>,LO,GO,N> >& A,
    const Teuchos::RCP<Teuchos::ParameterList>& precParams,
    const Teuchos::RCP<Tpetra::MultiVector<double,LO,GO,N> >& coords)
  {
    typedef Sacado::UQ::PCE<S> Scalar;
    typedef typename Scalar::value_type BaseScalar;

    using Teuchos::RCP;
    using Teuchos::rcp;

#if USE_SCALAR_MEAN_BASED_PREC
    RCP<Tpetra::CrsMatrix<BaseScalar,LO,GO,N> > mean_scalar =
      build_mean_scalar_matrix(*A);
    RCP<Tpetra::Operator<BaseScalar,LO,GO,N> > prec_scalar =
      build_muelu_preconditioner(mean_scalar, precParams, coords);
    RCP<Tpetra::Operator<Scalar,LO,GO,N> > prec =
      rcp(new Stokhos::MeanBasedTpetraOperator<Scalar,LO,GO,N>(prec_scalar));
#else
    RCP<Tpetra::CrsMatrix<Scalar,LO,GO,N> > mean =
      Stokhos::build_mean_matrix(*A);
    RCP<Tpetra::Operator<Scalar,LO,GO,N> > prec =
      build_muelu_preconditioner(mean, precParams, coords);
#endif

    return prec;
  }

}
}

#endif

#include <fenl_impl.hpp>
