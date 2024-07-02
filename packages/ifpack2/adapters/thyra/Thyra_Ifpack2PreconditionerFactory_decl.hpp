// @HEADER
// *****************************************************************************
//       Ifpack2: Templated Object-Oriented Algebraic Preconditioner Package
//
// Copyright 2009 NTESS and the Ifpack2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER
#ifndef THYRA_IFPACK2_PRECONDITIONERFACTORY_DECL_HPP
#define THYRA_IFPACK2_PRECONDITIONERFACTORY_DECL_HPP

#include "Thyra_PreconditionerFactoryBase.hpp"

namespace Thyra {

/** \brief Concrete preconditioner factory subclass based on Ifpack2.
 */
template <typename MatrixType>
class Ifpack2PreconditionerFactory :
  public PreconditionerFactoryBase<typename MatrixType::scalar_type> {
public:
  /** \brief . */
  typedef MatrixType fwd_op_type;

  /** \brief . */
  typedef typename MatrixType::scalar_type scalar_type;

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief . */
  Ifpack2PreconditionerFactory();
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

  // ToDo: Add an override of describe(...) to give more detail!

  //@}

private:

  Teuchos::RCP<Teuchos::ParameterList> paramList_;

};

} // namespace Thyra

#endif // THYRA_IFPACK2_PRECONDITIONERFACTORY_DECL_HPP
