// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_TEKOMULTIPHYSBLOCKPRECONDITIONERFACTORY_HPP
#define MUELU_TEKOMULTIPHYSBLOCKPRECONDITIONERFACTORY_HPP

#include "MueLu_ConfigDefs.hpp"

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA) && defined(HAVE_MUELU_THYRATPETRAADAPTERS) && defined(HAVE_MUELU_TEKO)

#include "Teko_BlockPreconditionerFactory.hpp"
#include "Teko_BlockImplicitLinearOp.hpp"
#include "Teko_TpetraHelpers.hpp"
#include "Teko_Utilities.hpp"

#include "MueLu_MultiPhys_decl.hpp"

#include "Xpetra_BlockedCrsMatrix_decl.hpp"
#include "Xpetra_MapExtractor_fwd.hpp"
#include "Xpetra_MultiVector_decl.hpp"

#include "Teuchos_ParameterList.hpp"
#include "Teuchos_RCP.hpp"

namespace MueLu {

class TekoMultiPhysBlockPreconditionerState : public Teko::BlockPreconditionerState {
 public:
  using ST                  = Teko::ST;
  using LO                  = Teko::LO;
  using GO                  = Teko::GO;
  using NT                  = Teko::NT;
  using TpetraMatrix        = Tpetra::CrsMatrix<ST, LO, GO, NT>;
  using XpetraBlockedMatrix = Xpetra::BlockedCrsMatrix<ST, LO, GO, NT>;
  using MultiPhysOperator   = MueLu::MultiPhys<ST, LO, GO, NT>;

  TekoMultiPhysBlockPreconditionerState();

  Teuchos::RCP<MultiPhysOperator> multiphysPreconditioner;
  Teuchos::RCP<TpetraMatrix> monolithicMatrix;
  Teuchos::RCP<XpetraBlockedMatrix> blockedMatrix;
  int numBlocks;
};

class TekoMultiPhysBlockPreconditionerFactory : public Teko::BlockPreconditionerFactory {
 public:
  using ST                  = Teko::ST;
  using LO                  = Teko::LO;
  using GO                  = Teko::GO;
  using NT                  = Teko::NT;
  using XpetraMultiVector   = Xpetra::MultiVector<ST, LO, GO, NT>;
  using XpetraMatrix        = Xpetra::Matrix<ST, LO, GO, NT>;
  using XpetraBlockedMatrix = Xpetra::BlockedCrsMatrix<ST, LO, GO, NT>;
  using XpetraMapExtractor  = Xpetra::MapExtractor<ST, LO, GO, NT>;
  using MultiPhysOperator   = MueLu::MultiPhys<ST, LO, GO, NT>;

  TekoMultiPhysBlockPreconditionerFactory();

  Teko::LinearOp buildPreconditionerOperator(Teko::BlockedLinearOp& blo,
                                             Teko::BlockPreconditionerState& state) const override;

  Teuchos::RCP<Teko::PreconditionerState> buildPreconditionerState() const override;

 public:
  void initializeFromParameterList(const Teuchos::ParameterList& pl) override;

 private:
  Teuchos::RCP<XpetraMultiVector> convertToXpetraMultiVector(const Teuchos::ParameterList& pl,
                                                             const std::string& name) const;

  Teuchos::RCP<TekoMultiPhysBlockPreconditionerState> ensureConcreteState(
      Teko::BlockPreconditionerState& state) const;

  Teuchos::RCP<Teuchos::ParameterList> amgSettings_;
  std::vector<Teuchos::RCP<XpetraMultiVector> > coordinates_;
  std::vector<Teuchos::RCP<XpetraMultiVector> > material_;
  std::vector<Teuchos::RCP<XpetraMultiVector> > nullspace_;
  bool useBlockSmoothers_;
  bool printComplexity_;
};

void registerTekoMultiPhysBlockPreconditionerFactory(const std::string& name = "Block AMG");

}  // namespace MueLu

#endif

#endif