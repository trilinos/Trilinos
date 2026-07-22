// @HEADER
// *****************************************************************************
//      Teko: A package for block and physics based preconditioning
//
// Copyright 2010 NTESS and the Teko contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Teko_NeumannSeriesPreconditionerFactoryDecl_hpp__
#define __Teko_NeumannSeriesPreconditionerFactoryDecl_hpp__

#include "Teuchos_ParameterListAcceptor.hpp"

// Thyra includes
#include "Thyra_PreconditionerFactoryBase.hpp"

// Teko includes
#include "Teko_Utilities.hpp"

namespace Teko {

using Teuchos::RCP;

template <typename ScalarT>
class NeumannSeriesPreconditionerFactory
    : public virtual Thyra::PreconditionerFactoryBase<ScalarT> {
 public:
  NeumannSeriesPreconditionerFactory();

  //! is this operator compatiable with the preconditioner factory?
  bool isCompatible(const Thyra::LinearOpSourceBase<ScalarT> &fwdOpSrc) const;

  //! create an instance of the preconditioner
  RCP<Thyra::PreconditionerBase<ScalarT> > createPrec() const;

  /** \brief initialize a newly created preconditioner object
   *
   * Initialize a newly created preconditioner object. For use with
   * nonlinear solvers.
   *
   * \param[in] fwdOpSrc Forward operator to be preconditioned
   * \param[in] solnVec Vector associated with this linear operator.
   * \param[in,out] precOp Return location for the preconditioner
   * \param[in] supportSolveUse Thyra information (?)
   */
  void initializePrec(const RCP<const Thyra::LinearOpSourceBase<ScalarT> > &fwdOpSrc,
                      const RCP<const Thyra::MultiVectorBase<ScalarT> > &solnVec,
                      Thyra::PreconditionerBase<ScalarT> *precOp,
                      const Thyra::ESupportSolveUse supportSolveUse) const;

  /** \brief initialize a newly created preconditioner object
   *
   * Initialize a newly created preconditioner object.
   *
   * \param[in] fwdOpSrc Forward operator to be preconditioned
   * \param[in,out] precOp Return location for the preconditioner
   * \param[in] supportSolveUse Thyra information (?)
   */
  void initializePrec(const RCP<const Thyra::LinearOpSourceBase<ScalarT> > &fwdOpSrc,
                      Thyra::PreconditionerBase<ScalarT> *precOp,
                      const Thyra::ESupportSolveUse supportSolveUse) const;

  //! wipe clean a already initialized preconditioner object
  void uninitializePrec(Thyra::PreconditionerBase<ScalarT> *prec,
                        RCP<const Thyra::LinearOpSourceBase<ScalarT> > *fwdOpSrc,
                        Thyra::ESupportSolveUse *supportSolveUse) const;

  /** @name Overridden from Teuchos::ParameterListAcceptor */
  //@{

  //! \brief Set parameters from a parameter list
  void setParameterList(Teuchos::RCP<Teuchos::ParameterList> const &paramList);

  //! \brief Get the parameter list that was set using setParameterList().
  Teuchos::RCP<Teuchos::ParameterList> getNonconstParameterList();

  //! \brief Unset the parameter list that was set using setParameterList().
  Teuchos::RCP<Teuchos::ParameterList> unsetParameterList();

  //! \brief Get the parameter list that was set using setParameterList().
  Teuchos::RCP<const Teuchos::ParameterList> getParameterList() const;

  /** \brief Get the valid parameters */
  Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;
  //@}

  /** \name Public functions overridden from Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

 protected:
  //! for ParameterListAcceptor
  Teuchos::RCP<Teuchos::ParameterList> paramList_;

  int numberOfTerms_;
  Teko::DiagonalType scalingType_;
};

}  // end namespace Teko

#endif
