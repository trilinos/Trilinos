// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef MUELU_XPETRA_THYRALINEAROP_HPP
#define MUELU_XPETRA_THYRALINEAROP_HPP

#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Operator.hpp>
#include <Xpetra_MultiVector.hpp>

// Stratimikos
#include <Thyra_VectorBase.hpp>
#include <Thyra_SolveSupportTypes.hpp>
#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Teuchos_AbstractFactoryStd.hpp>
#include <Stratimikos_LinearSolverBuilder.hpp>
#include <Thyra_MueLuPreconditionerFactory.hpp>

namespace MueLu {

/*!  @brief Wraps an existing Thyra::LinearOp as a Xpetra::Operator.
 */
template <class Scalar        = DefaultScalar,
          class LocalOrdinal  = DefaultLocalOrdinal,
          class GlobalOrdinal = DefaultGlobalOrdinal,
          class Node          = DefaultNode>
class XpetraThyraLinearOp : public Xpetra::Operator<Scalar, LocalOrdinal, GlobalOrdinal, Node> {
 protected:
  XpetraThyraLinearOp() = default;

 public:
  //! @name Constructor/Destructor
  //@{

  //! Constructor
  XpetraThyraLinearOp(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A,
                      RCP<ParameterList> params)
    : A_(A) {
    throw Exceptions::RuntimeError("Interface not supported");
  };

  //! Destructor.
  ~XpetraThyraLinearOp() = default;

  //@}

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> getDomainMap() const {
    throw Exceptions::RuntimeError("Interface not supported");
  }

  // //! Returns the Tpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> getRangeMap() const {
    throw Exceptions::RuntimeError("Interface not supported");
  }

  //! Returns in Y the result of a Xpetra::Operator applied to a Xpetra::MultiVector X.
  /*!
    \param[in]  X - Xpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param[out] Y - Xpetra::MultiVector of dimension NumVectors containing result.
  */
  void apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
             Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta           = Teuchos::ScalarTraits<Scalar>::one()) const {
    throw Exceptions::RuntimeError("Interface not supported");
  }

  //! Indicates whether this operator supports applying the adjoint operator.
  bool hasTransposeApply() const { return false; }

  //! Compute a residual R = B - (*this) * X
  void residual(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& R) const {
    throw Exceptions::RuntimeError("Interface not supported");
  }

 private:
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A_;
};

// Partial specialization for Scalar == double.
// Allows to avoid issues with Stokhos instantiating Thyra objects.
template <class LocalOrdinal,
          class GlobalOrdinal,
          class Node>
class XpetraThyraLinearOp<double, LocalOrdinal, GlobalOrdinal, Node> : public Xpetra::Operator<double, LocalOrdinal, GlobalOrdinal, Node> {
  using Scalar = double;

 protected:
  XpetraThyraLinearOp() = default;

 public:
  //! @name Constructor/Destructor
  //@{

  //! Constructor
  XpetraThyraLinearOp(RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A,
                      RCP<ParameterList> params)
    : A_(A) {
    // Build Thyra linear algebra objects
    RCP<const Thyra::LinearOpBase<Scalar>> thyraA = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyra(Teuchos::rcp_dynamic_cast<Xpetra::CrsMatrixWrap<Scalar, LocalOrdinal, GlobalOrdinal, Node>>(A)->getCrsMatrix());

    Stratimikos::LinearSolverBuilder<Scalar> linearSolverBuilder;
    typedef Thyra::PreconditionerFactoryBase<Scalar> Base;
    typedef Thyra::MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> ImplMueLu;
    linearSolverBuilder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, ImplMueLu>(), "MueLu");

    linearSolverBuilder.setParameterList(params);

    // Build a new "solver factory" according to the previously specified parameter list.
    // RCP<Thyra::LinearOpWithSolveFactoryBase<Scalar> > solverFactory = Thyra::createLinearSolveStrategy(linearSolverBuilder);

    auto precFactory = Thyra::createPreconditioningStrategy(linearSolverBuilder);
    auto prec        = precFactory->createPrec();

    precFactory->initializePrec(Thyra::defaultLinearOpSource(thyraA), prec.get(), Thyra::SUPPORT_SOLVE_UNSPECIFIED);
    prec_ = prec;
  };

  //! Destructor.
  ~XpetraThyraLinearOp() = default;

  //@}

  //! Returns the Tpetra::Map object associated with the domain of this operator.
  const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> getDomainMap() const {
    return A_->getDomainMap();
  }

  // //! Returns the Tpetra::Map object associated with the range of this operator.
  const Teuchos::RCP<const Xpetra::Map<LocalOrdinal, GlobalOrdinal, Node>> getRangeMap() const {
    return A_->getRangeMap();
  }

  //! Returns in Y the result of a Xpetra::Operator applied to a Xpetra::MultiVector X.
  /*!
    \param[in]  X - Xpetra::MultiVector of dimension NumVectors to multiply with matrix.
    \param[out] Y - Xpetra::MultiVector of dimension NumVectors containing result.
  */
  void apply(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
             Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& Y,
             Teuchos::ETransp mode = Teuchos::NO_TRANS,
             Scalar alpha          = Teuchos::ScalarTraits<Scalar>::one(),
             Scalar beta           = Teuchos::ScalarTraits<Scalar>::one()) const {
    RCP<const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rcpX = Teuchos::rcpFromRef(X);
    RCP<const Thyra::MultiVectorBase<Scalar>> thyraX                               = Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyraMultiVector(rcpX);

    RCP<Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>> rcpY = Teuchos::rcpFromRef(Y);
    RCP<Thyra::MultiVectorBase<Scalar>> thyraY                               = Teuchos::rcp_const_cast<Thyra::MultiVectorBase<Scalar>>(Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toThyraMultiVector(rcpY));

    prec_->getUnspecifiedPrecOp()->apply(Thyra::NOTRANS, *thyraX, thyraY.ptr(), alpha, beta);
    Y = *Xpetra::ThyraUtils<Scalar, LocalOrdinal, GlobalOrdinal, Node>::toXpetra(thyraY, Y.getMap()->getComm());
  }

  //! Indicates whether this operator supports applying the adjoint operator.
  bool hasTransposeApply() const { return false; }

  //! Compute a residual R = B - (*this) * X
  void residual(const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& X,
                const Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& B,
                Xpetra::MultiVector<Scalar, LocalOrdinal, GlobalOrdinal, Node>& R) const {
    using STS = Teuchos::ScalarTraits<Scalar>;
    R.update(STS::one(), B, STS::zero());
    this->apply(X, R, Teuchos::NO_TRANS, -STS::one(), STS::one());
  }

 private:
  RCP<Xpetra::Matrix<Scalar, LocalOrdinal, GlobalOrdinal, Node>> A_;
  RCP<Thyra::PreconditionerBase<Scalar>> prec_;
};

}  // namespace MueLu

#endif  // defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)

#endif  // MUELU_XPETRA_THYRALINEAROP_HPP
