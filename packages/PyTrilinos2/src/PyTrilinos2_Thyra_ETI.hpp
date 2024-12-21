// @HEADER
// *****************************************************************************
//          PyTrilinos2: Automatic Python Interfaces to Trilinos Packages
//
// Copyright 2022 NTESS and the PyTrilinos2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PYTRILINOS2_THYRA_ETI
#define PYTRILINOS2_THYRA_ETI

#include <Thyra_TpetraThyraWrappers.hpp>
#include <Thyra_TpetraEuclideanScalarProd.hpp>
#include <Thyra_TpetraVectorSpace.hpp>
#include <Thyra_TpetraMultiVector.hpp>
#include <Thyra_TpetraVector.hpp>
#include <Thyra_TpetraLinearOp.hpp>

#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearSolverBuilderBase.hpp>
#include <Thyra_VectorBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_OperatorVectorTypes.hpp>


#define BINDER_THYRA_WRAPPERS_INSTANT(SC,LO,GO,NO) \
  template Teuchos::RCP<Thyra::VectorBase<SC> > createVector<SC, LO, GO, NO>(const Teuchos::RCP<Tpetra::Vector<SC, LO, GO, NO> > &tpetraVector, \
                                                                             const Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > space); \
  template Teuchos::RCP<Thyra::MultiVectorBase<SC> > createMultiVector<SC, LO, GO, NO>(const Teuchos::RCP<Tpetra::MultiVector<SC, LO, GO, NO> > &tpetraMultiVector, \
                                                                                       const Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > rangeSpace, \
                                                                                       const Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > domainSpace);

#define BINDER_THYRA_EUCLIDEANSCALARPROD_INSTANT(SC,LO,GO,NO) \
  template class Thyra::TpetraEuclideanScalarProd<SC, LO, GO, NO>;

#define BINDER_THYRA_VECTORSPACE_INSTANT(SC,LO,GO,NO) \
  template class Thyra::TpetraVectorSpace<SC, LO, GO, NO>;

#define BINDER_THYRA_MULTIVECTOR_INSTANT(SC,LO,GO,NO) \
  template class Thyra::TpetraMultiVector<SC, LO, GO, NO>;

#define BINDER_THYRA_VECTOR_INSTANT(SC, LO, GO, NO) \
  template class Thyra::TpetraVector<SC, LO, GO, NO>;

#define BINDER_THYRA_LINEAROP_INSTANT(SC,LO,GO,NO) \
  template class Thyra::TpetraLinearOp<SC, LO, GO, NO>; \
  template Teuchos::RCP<Thyra::TpetraLinearOp<SC, LO, GO, NO> > tpetraLinearOp(const Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > &rangeSpace, \
                                                                               const Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > &domainSpace, \
                                                                               const Teuchos::RCP<Tpetra::Operator<SC, LO, GO, NO> > &tpetraOperator); \
  template Teuchos::RCP<const Thyra::TpetraLinearOp<SC, LO, GO, NO> > constTpetraLinearOp(const Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > &rangeSpace, \
                                                                                          const Teuchos::RCP<const Thyra::VectorSpaceBase<SC> > &domainSpace, \
                                                                                          const Teuchos::RCP<const Tpetra::Operator<SC, LO, GO, NO> > &tpetraOperator);

#define BINDER_THYRA_LINEARSOLVER_INSTANT(SC) \
  template class Thyra::LinearOpWithSolveFactoryBase<SC>; \
  template Teuchos::RCP<Thyra::LinearOpWithSolveFactoryBase<SC> > Thyra::createLinearSolveStrategy(const Thyra::LinearSolverBuilderBase<SC> &linearSolverBuilder, const std::string &linearSolveStrategyName); \
  template Teuchos::RCP<Thyra::LinearOpWithSolveBase<SC> > linearOpWithSolve(const Thyra::LinearOpWithSolveFactoryBase<SC> &lowsFactory, \
                                                                             const Teuchos::RCP<const Thyra::LinearOpBase<SC> > &fwdOp, \
                                                                             const Thyra::ESupportSolveUse supportSolveUse); \
  template struct Thyra::SolveCriteria<SC>; \
  template struct Thyra::SolveStatus<SC>;

#define BINDER_TPETRA_OPERATOR_INSTANT(SC,LO,GO,NO) \
   template class Tpetra::Operator<SC, LO, GO, NO>;

namespace Tpetra {
  BINDER_TPETRA_OPERATOR_INSTANT(Tpetra::Details::DefaultTypes::scalar_type, Tpetra::Details::DefaultTypes::local_ordinal_type, Tpetra::Details::DefaultTypes::global_ordinal_type, Tpetra::Details::DefaultTypes::node_type)
}


namespace Thyra {

    template <typename T>
    void initiate(T) {};

  BINDER_THYRA_WRAPPERS_INSTANT(Tpetra::Details::DefaultTypes::scalar_type, Tpetra::Details::DefaultTypes::local_ordinal_type, Tpetra::Details::DefaultTypes::global_ordinal_type, Tpetra::Details::DefaultTypes::node_type)
  BINDER_THYRA_EUCLIDEANSCALARPROD_INSTANT(Tpetra::Details::DefaultTypes::scalar_type, Tpetra::Details::DefaultTypes::local_ordinal_type, Tpetra::Details::DefaultTypes::global_ordinal_type, Tpetra::Details::DefaultTypes::node_type)
  BINDER_THYRA_VECTORSPACE_INSTANT(Tpetra::Details::DefaultTypes::scalar_type, Tpetra::Details::DefaultTypes::local_ordinal_type, Tpetra::Details::DefaultTypes::global_ordinal_type, Tpetra::Details::DefaultTypes::node_type)
  BINDER_THYRA_MULTIVECTOR_INSTANT(Tpetra::Details::DefaultTypes::scalar_type, Tpetra::Details::DefaultTypes::local_ordinal_type, Tpetra::Details::DefaultTypes::global_ordinal_type, Tpetra::Details::DefaultTypes::node_type)
  BINDER_THYRA_VECTOR_INSTANT(Tpetra::Details::DefaultTypes::scalar_type, Tpetra::Details::DefaultTypes::local_ordinal_type, Tpetra::Details::DefaultTypes::global_ordinal_type, Tpetra::Details::DefaultTypes::node_type)
  BINDER_THYRA_LINEAROP_INSTANT(Tpetra::Details::DefaultTypes::scalar_type, Tpetra::Details::DefaultTypes::local_ordinal_type, Tpetra::Details::DefaultTypes::global_ordinal_type, Tpetra::Details::DefaultTypes::node_type)
  BINDER_THYRA_LINEARSOLVER_INSTANT(Tpetra::Details::DefaultTypes::scalar_type)
}

#endif // PYTRILINOS2_THYRA_ETI
