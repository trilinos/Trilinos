// @HEADER
// *****************************************************************************
//         Stratimikos: Thyra-based strategies for linear solvers
//
// Copyright 2006 NTESS and the Stratimikos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_BELOS_TPETRA_PRECONDITIONERFACTORY_DECL_HPP
#define THYRA_BELOS_TPETRA_PRECONDITIONERFACTORY_DECL_HPP

#include "Thyra_PreconditionerFactoryBase.hpp"

#include "Tpetra_Operator.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include <functional>
#include <map>
#include <string>

namespace Thyra {

/** \brief Concrete preconditioner factory subclass based on Belos.
 * (Yes, Belos solvers can also be used as preconditioners!)
 *
 * The inner Belos solve can optionally be given its own algebraic
 * preconditioner, selected from XML via the "Inner Prec Type" parameter.
 * Builders for concrete preconditioner packages are looked up in a static
 * registry.  The Ifpack2 builder is registered automatically when Stratimikos
 * is configured with Ifpack2.  Packages that live *downstream* of Stratimikos
 * (e.g. MueLu, which optionally depends on Stratimikos) cannot be called
 * directly from here without creating a circular dependency; instead they
 * inject their own builder via setInnerPrecBuilder() -- see
 * Stratimikos::enableBelosPrecMueLu().
 */
template <typename MatrixType>
class BelosTpetraPreconditionerFactory :
  public PreconditionerFactoryBase<typename MatrixType::scalar_type> {
public:
  /** \brief . */
  typedef typename MatrixType::scalar_type scalar_type;
  /** \brief . */
  typedef typename MatrixType::local_ordinal_type local_ordinal_type;
  /** \brief . */
  typedef typename MatrixType::global_ordinal_type global_ordinal_type;
  /** \brief . */
  typedef typename MatrixType::node_type node_type;

  //! Tpetra operator type of the inner algebraic preconditioner.
  typedef Tpetra::Operator<scalar_type, local_ordinal_type, global_ordinal_type, node_type>
    inner_prec_op_type;

  //! Half-of-standard-precision scalar (e.g. float when scalar_type is double).
  typedef typename Teuchos::ScalarTraits<scalar_type>::halfPrecision half_scalar_type;

  //! Tpetra operator type of a *half-precision* inner algebraic preconditioner.
  typedef Tpetra::Operator<half_scalar_type, local_ordinal_type, global_ordinal_type, node_type>
    inner_prec_op_half_type;

  /** \brief Callback that builds an inner algebraic preconditioner operator.
   *
   * Arguments: the (const) forward Tpetra operator of the inner solve, and the
   * "Inner Prec Params" parameter list.  Returns the preconditioner operator to
   * install on the inner Belos LinearProblem.
   */
  typedef std::function<
    Teuchos::RCP<const inner_prec_op_type>(
      const Teuchos::RCP<const inner_prec_op_type>&,
      Teuchos::ParameterList&)
    > InnerPrecBuilder;

  /** \brief Half-precision analog of InnerPrecBuilder, used when the inner Belos
   * solve is built at half_scalar_type ("half precision"=true).  The forward
   * operator and the returned preconditioner are both at half_scalar_type. */
  typedef std::function<
    Teuchos::RCP<const inner_prec_op_half_type>(
      const Teuchos::RCP<const inner_prec_op_half_type>&,
      Teuchos::ParameterList&)
    > InnerPrecBuilderHalf;

  /** \brief Register (or replace) the inner-preconditioner builder used for a
   * given "Inner Prec Type" name.  Intended for downstream packages to inject
   * their own preconditioners (e.g. MueLu). */
  static void setInnerPrecBuilder(const std::string& precType, InnerPrecBuilder builder);

  /** \brief Whether a builder is registered for the given "Inner Prec Type". */
  static bool hasInnerPrecBuilder(const std::string& precType);

  /** \brief Register (or replace) the *half-precision* inner-preconditioner
   * builder used for a given "Inner Prec Type" name, applied when the inner
   * Belos solve runs at half precision.  Downstream packages inject their own
   * (e.g. Stratimikos::enableBelosPrecMueLuHalf()); note this requires the
   * package to be ETI-instantiated for half_scalar_type. */
  static void setInnerPrecBuilderHalf(const std::string& precType, InnerPrecBuilderHalf builder);

  /** \brief Whether a half-precision builder is registered for the given name. */
  static bool hasInnerPrecBuilderHalf(const std::string& precType);

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  BelosTpetraPreconditionerFactory();
  //@}

  /** @name Overridden from PreconditionerFactoryBase */
  //@{

  /** \brief . */
  bool isCompatible(const LinearOpSourceBase<scalar_type> &fwdOp) const;

  /** \brief . */
  Teuchos::RCP<PreconditionerBase<scalar_type> > createPrec() const;

  /** \brief . */
  void initializePrec(
    const Teuchos::RCP<const LinearOpSourceBase<scalar_type> > &fwdOp,
    PreconditionerBase<scalar_type> *prec,
    const ESupportSolveUse supportSolveUse
    ) const;

  /** \brief . */
  void uninitializePrec(
    PreconditionerBase<scalar_type> *prec,
    Teuchos::RCP<const LinearOpSourceBase<scalar_type> > *fwdOp,
    ESupportSolveUse *supportSolveUse
    ) const;

  //@}

  /** @name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  /** \brief . */
  void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> &paramList);
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();
  /** \brief . */
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;
  /** \brief . */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

  /** \name Public functions overridden from Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

private:

  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  //! Registry of inner-preconditioner builders, keyed by "Inner Prec Type".
  static std::map<std::string, InnerPrecBuilder>& innerPrecRegistry();

  //! Half-precision analog of innerPrecRegistry().
  static std::map<std::string, InnerPrecBuilderHalf>& innerPrecRegistryHalf();

};

} // namespace Thyra

#endif // THYRA_BELOS_TPETRA_PRECONDITIONERFACTORY_DECL_HPP
