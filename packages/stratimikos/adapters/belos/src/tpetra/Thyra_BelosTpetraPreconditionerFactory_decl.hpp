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

namespace Thyra {

/** \brief Concrete preconditioner factory subclass based on Belos.
 * (Yes, Belos solvers can also be used as preconditioners!)
 */
template <typename MatrixType>
class BelosTpetraPreconditionerFactory :
  public PreconditionerFactoryBase<typename MatrixType::scalar_type> {
public:
  /** \brief . */
  typedef typename MatrixType::scalar_type scalar_type;

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

};

} // namespace Thyra

#endif // THYRA_BELOS_TPETRA_PRECONDITIONERFACTORY_DECL_HPP
