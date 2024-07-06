// @HEADER
// *****************************************************************************
//          PyTrilinos2: Automatic Python Interfaces to Trilinos Packages
//
// Copyright 2022 NTESS and the PyTrilinos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PYTRILINOS2_MUELU_ETI
#define PYTRILINOS2_MUELU_ETI

#include <MueLu_CreateTpetraPreconditioner.hpp>

#define BINDER_MUELU_CREATETPETRAPRECONDITIONER_INSTANT(SCALAR,LO,GO,NO) \
  template Teuchos::RCP<MueLu::TpetraOperator<SCALAR, LO, GO, NO> > CreateTpetraPreconditioner<SCALAR, LO, GO, NO>(const Teuchos::RCP<Tpetra::Operator<SCALAR, LO, GO, NO> > &inA, Teuchos::ParameterList& inParamList);

namespace MueLu {

    template <typename T>
    void initiate(T) {};

  BINDER_MUELU_CREATETPETRAPRECONDITIONER_INSTANT(Tpetra::Details::DefaultTypes::scalar_type, Tpetra::Details::DefaultTypes::local_ordinal_type, Tpetra::Details::DefaultTypes::global_ordinal_type, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType)
}

#endif // PYTRILINOS2_MUELU_ETI
