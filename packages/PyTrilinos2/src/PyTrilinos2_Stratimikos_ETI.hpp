// @HEADER
// *****************************************************************************
//          PyTrilinos2: Automatic Python Interfaces to Trilinos Packages
//
// Copyright 2022 NTESS and the PyTrilinos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PYTRILINOS2_STRATIMIKOS_ETI
#define PYTRILINOS2_STRATIMIKOS_ETI

#include <PyTrilinos2_config.hpp>
#include <Stratimikos_LinearSolverBuilder.hpp>
#ifdef HAVE_PYTRILINOS2_MUELU
#include <Stratimikos_MueLuHelpers.hpp>
#endif


#define BINDER_STRATIMIKOS_LINEARSOLVERBUILDER_INSTANT(SC) \
  template class Stratimikos::LinearSolverBuilder<SC>;

#define BINDER_STRATIMIKOS_MUELU_INSTANT(SCALAR,LO,GO,NO) \
  template void enableMueLu<SCALAR, LO, GO, NO>(Stratimikos::LinearSolverBuilder<SCALAR>& builder, const std::string& stratName); \
  template void enableMueLuRefMaxwell<SCALAR, LO, GO, NO>(Stratimikos::LinearSolverBuilder<SCALAR>& builder, const std::string& stratName); \
  template void enableMueLuMaxwell1<SCALAR, LO, GO, NO>(Stratimikos::LinearSolverBuilder<SCALAR>& builder, const std::string& stratName);


namespace Stratimikos {

  template <typename T>
    void initiate(T) {};

  BINDER_STRATIMIKOS_LINEARSOLVERBUILDER_INSTANT(Tpetra::Details::DefaultTypes::scalar_type)

#ifdef HAVE_PYTRILINOS2_MUELU
  BINDER_STRATIMIKOS_MUELU_INSTANT(Tpetra::Details::DefaultTypes::scalar_type, Tpetra::Details::DefaultTypes::local_ordinal_type, Tpetra::Details::DefaultTypes::global_ordinal_type, Tpetra::KokkosClassic::DefaultNode::DefaultNodeType)
#endif

}

#endif // PYTRILINOS2_STRATIMIKOS_ETI
