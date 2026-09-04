// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOSTPETRAETIHELPERS_HPP
#define BELOSTPETRAETIHELPERS_HPP

#ifdef HAVE_BELOS_TPETRA

#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Operator.hpp"
#include "TpetraCore_ETIHelperMacros.h"


#define BELOS_TPETRA_CALL(CLASS, SC, LO, GO, NT)                               \
  template class CLASS<SC, Tpetra::MultiVector<SC, LO, GO, NT>,                \
                       Tpetra::Operator<SC, LO, GO, NT>,                       \
                       Teuchos::SerialDenseMatrix<int, SC>>;                   \
  template class CLASS<                                                        \
      SC, Tpetra::MultiVector<SC, LO, GO, NT>,                                 \
      Tpetra::Operator<SC, LO, GO, NT>,                                        \
      typename Tpetra::MultiVector<SC, LO, GO,                                 \
                                   NT>::wrapped_dual_view_type::DVT>;

#define BELOS_TPETRA_EXTERN_CALL(CLASS, SC, LO, GO, NT)                        \
  extern template class CLASS<SC, Tpetra::MultiVector<SC, LO, GO, NT>,         \
                              Tpetra::Operator<SC, LO, GO, NT>,                \
                              Teuchos::SerialDenseMatrix<int, SC>>;            \
  extern template class CLASS<                                                 \
      SC, Tpetra::MultiVector<SC, LO, GO, NT>,                                 \
      Tpetra::Operator<SC, LO, GO, NT>,                                        \
      typename Tpetra::MultiVector<SC, LO, GO,                                 \
                                   NT>::wrapped_dual_view_type::DVT>;

TPETRA_ETI_MANGLING_TYPEDEFS()

#endif

#endif
