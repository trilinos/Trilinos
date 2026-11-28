// @HEADER
// *****************************************************************************
//        MueLu: A package for multigrid based preconditioning
//
// Copyright 2012 NTESS and the MueLu contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef STRATIMIKOS_MUELU_TPETRA_HELPERS_HPP
#define STRATIMIKOS_MUELU_TPETRA_HELPERS_HPP

#include "MueLu_Details_DefaultTypes.hpp"
#include "Stratimikos_LinearSolverBuilder.hpp"

#include "Thyra_MueLuPreconditionerFactory.hpp"
#include "Thyra_MueLuRefMaxwellPreconditionerFactory.hpp"
#include "Thyra_MueLuMaxwell1PreconditionerFactory.hpp"

#if defined(HAVE_MUELU_EXPERIMENTAL) && defined(HAVE_MUELU_TEKO)
#include "Thyra_MueLuTpetraQ2Q1PreconditionerFactory.hpp"
#endif

// Support for injecting a MueLu preconditioner into the innermost Belos solve
// of Stratimikos' "BelosPrecTpetra" preconditioner.  Guarded so the MueLu
// Stratimikos adapter still builds when Stratimikos is configured without the
// Belos/Tpetra pieces that supply the factory below.
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_BELOS) && defined(HAVE_MUELU_THYRA)
#include "Thyra_BelosTpetraPreconditionerFactory_decl.hpp"
#include "Thyra_BelosTpetraPreconditionerFactory_def.hpp"
#include "MueLu_CreateTpetraPreconditioner.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_Operator.hpp"
#endif

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include <string>

namespace Stratimikos {

template <typename Scalar = MueLu::DefaultScalar, typename LocalOrdinal = MueLu::DefaultLocalOrdinal, typename GlobalOrdinal = MueLu::DefaultGlobalOrdinal, typename Node = MueLu::DefaultNode>
void enableMueLu(LinearSolverBuilder<Scalar>& builder, const std::string& stratName = "MueLu") {
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
  const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

  TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                             "Stratimikos::enableMueLu cannot add \"" + stratName + "\" because it is already included in builder!");

  typedef Thyra::PreconditionerFactoryBase<Scalar> Base;
  typedef Thyra::MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> Impl;

  builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
#endif
}

template <typename Scalar = MueLu::DefaultScalar, typename LocalOrdinal = MueLu::DefaultLocalOrdinal, typename GlobalOrdinal = MueLu::DefaultGlobalOrdinal, typename Node = MueLu::DefaultNode>
void enableMueLuRefMaxwell(LinearSolverBuilder<Scalar>& builder, const std::string& stratName = "MueLuRefMaxwell") {
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
  const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

  TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                             "Stratimikos::enableMueLuRefMaxwell cannot add \"" + stratName + "\" because it is already included in builder!");

  typedef Thyra::PreconditionerFactoryBase<Scalar> Base;
  typedef Thyra::MueLuRefMaxwellPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> Impl;

  builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
#endif
}

template <typename Scalar = MueLu::DefaultScalar, typename LocalOrdinal = MueLu::DefaultLocalOrdinal, typename GlobalOrdinal = MueLu::DefaultGlobalOrdinal, typename Node = MueLu::DefaultNode>
void enableMueLuMaxwell1(LinearSolverBuilder<Scalar>& builder, const std::string& stratName = "MueLuMaxwell1") {
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
  const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

  TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                             "Stratimikos::enableMueLuRefMaxwell cannot add \"" + stratName + "\" because it is already included in builder!");

  typedef Thyra::PreconditionerFactoryBase<Scalar> Base;
  typedef Thyra::MueLuMaxwell1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> Impl;

  builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
#endif
}

#if defined(HAVE_MUELU_EXPERIMENTAL) && defined(HAVE_MUELU_TEKO)
#if 0
  // Dynamically register MueLu Tpetra adapters in Stratimikos
  void enableMueLuTpetraQ2Q1(DefaultLinearSolverBuilder &builder, const std::string &stratName = "MueLu");
#endif

template <typename Scalar = MueLu::DefaultScalar, typename LocalOrdinal = MueLu::DefaultLocalOrdinal, typename GlobalOrdinal = MueLu::DefaultGlobalOrdinal, typename Node = MueLu::DefaultNode>
void enableMueLuTpetraQ2Q1(LinearSolverBuilder<Scalar>& builder, const std::string& stratName = "MueLu") {
  const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

  TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                             "Stratimikos::enableMueLuTpetraQ2Q1 cannot add \"" + stratName + "\" because it is already included in builder!");

  typedef Thyra::PreconditionerFactoryBase<Scalar> Base;
  typedef Thyra::MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> Impl;

  builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
}

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
MUELU_DEPRECATED void enableMueLuTpetraQ2Q1(LinearSolverBuilder<double>& builder, const std::string& stratName = "MueLu") {
  enableMueLuTpetraQ2Q1<double, LocalOrdinal, GlobalOrdinal, Node>(builder, stratName);
}
#endif

/*! \brief Register a MueLu preconditioner as an inner algebraic preconditioner
    for Stratimikos' "BelosPrecTpetra" preconditioner (the innermost Belos solve).

    Stratimikos lives \e upstream of MueLu (MueLu optionally depends on
    Stratimikos), so the Belos-as-preconditioner factory cannot call MueLu
    directly without creating a circular package dependency.  Instead the factory
    exposes a small builder registry, and MueLu injects itself here (dependency
    inversion) via \c setInnerPrecBuilder.

    Call this once, before solving, on the same \c MatrixType used with
    \c enableBelosPrecTpetra (default \c Tpetra::CrsMatrix<Scalar,LO,GO,Node>).
    Afterwards, an XML block such as

    \code{.xml}
      <ParameterList name="BelosPrecTpetra">
        <Parameter name="BelosPrec Solver Type" type="string" value="GmresPoly"/>
        <Parameter name="Inner Prec Type"       type="string" value="MueLu"/>
        <Parameter name="Inner Prec Side"       type="string" value="left"/>
        <ParameterList name="Inner Prec Params"> ... MueLu parameters ... </ParameterList>
      </ParameterList>
    \endcode

    builds a MueLu hierarchy on the inner Belos operator and applies it as a
    left/right preconditioner of that inner solve.

    \param precName  Value of "Inner Prec Type" this builder answers to (default "MueLu").
*/
template <typename Scalar        = MueLu::DefaultScalar,
          typename LocalOrdinal  = MueLu::DefaultLocalOrdinal,
          typename GlobalOrdinal = MueLu::DefaultGlobalOrdinal,
          typename Node          = MueLu::DefaultNode>
void enableBelosPrecMueLu(const std::string& precName = "MueLu") {
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_BELOS) && defined(HAVE_MUELU_THYRA)
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> matrix_type;
  typedef Thyra::BelosTpetraPreconditionerFactory<matrix_type> factory_type;
  // Tpetra::Operator<Scalar,LO,GO,Node>: the inner solve's forward operator and
  // the returned preconditioner operator type.
  typedef typename factory_type::inner_prec_op_type op_type;

  factory_type::setInnerPrecBuilder(
      precName,
      [](const Teuchos::RCP<const op_type>& fwdOp,
         Teuchos::ParameterList& params) -> Teuchos::RCP<const op_type> {
        // const-operator overload of CreateTpetraPreconditioner; the returned
        // MueLu::TpetraOperator is-a Tpetra::Operator and upcasts to op_type.
        return MueLu::CreateTpetraPreconditioner<Scalar, LocalOrdinal, GlobalOrdinal, Node>(fwdOp, params);
      });
#endif
}

/*! \brief Register a MueLu preconditioner as the inner algebraic preconditioner
    of a \e half-precision "BelosPrecTpetra" inner Belos solve
    ("half precision"=true).

    This is the half-precision analog of enableBelosPrecMueLu().  It is a
    \e separate, opt-in call because it forces MueLu to be instantiated for the
    half scalar (e.g. \c float when \c Scalar is \c double): calling it in a
    build whose MueLu ETI does not include the half scalar will fail to link.
    enableBelosPrecMueLu() (full precision) never pulls in that instantiation, so
    it stays safe in double-only MueLu builds.

    The whole body compiles to nothing unless a half-precision scalar is
    available in the build (THYRA_BELOS_PREC_ENABLE_HALF_PRECISION), matching the
    factory's own half-precision path.  Call it (in addition to, or instead of,
    enableBelosPrecMueLu()) before solving with "half precision"=true and
    "Inner Prec Type"="MueLu".

    \param precName  Value of "Inner Prec Type" this builder answers to (default "MueLu").
*/
template <typename Scalar        = MueLu::DefaultScalar,
          typename LocalOrdinal  = MueLu::DefaultLocalOrdinal,
          typename GlobalOrdinal = MueLu::DefaultGlobalOrdinal,
          typename Node          = MueLu::DefaultNode>
void enableBelosPrecMueLuHalf(const std::string& precName = "MueLu") {
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_BELOS) && defined(HAVE_MUELU_THYRA) && \
    defined(THYRA_BELOS_PREC_ENABLE_HALF_PRECISION)
  typedef Tpetra::CrsMatrix<Scalar, LocalOrdinal, GlobalOrdinal, Node> matrix_type;
  typedef Thyra::BelosTpetraPreconditionerFactory<matrix_type> factory_type;
  // Half-precision scalar (e.g. float) and the matching inner operator type.
  typedef typename factory_type::half_scalar_type half_scalar_type;
  typedef typename factory_type::inner_prec_op_half_type op_half_type;

  factory_type::setInnerPrecBuilderHalf(
      precName,
      [](const Teuchos::RCP<const op_half_type>& fwdOp,
         Teuchos::ParameterList& params) -> Teuchos::RCP<const op_half_type> {
        // MueLu built at the half scalar; requires MueLu ETI for half_scalar_type.
        return MueLu::CreateTpetraPreconditioner<half_scalar_type, LocalOrdinal, GlobalOrdinal, Node>(fwdOp, params);
      });
#endif
}

}  // namespace Stratimikos

#endif
