// @HEADER
// *****************************************************************************
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//
// Copyright 2004 NTESS and the Thyra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
#define THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP


#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_LinearOpSourceBase.hpp"


namespace Thyra {


/** \brief General delayed construction LinearOpWithSolveFactoryBase subclass.
 *
 * This simple decorator class allows for the delayed construction of the
 * linear solver until the last possible moment when the linear solver is
 * needed.  This class creates LinearOpWithSolveBase objects of type
 * DelayedLinearOpWithSolve.  The class object DelayedLinearOpWithSolve
 * actually implements the delayed linear solver construction.
 */
template<class Scalar>
class DelayedLinearOpWithSolveFactory
  : virtual public LinearOpWithSolveFactoryBase<Scalar>
{
public:

  /** @name Overridden from Constructors/Initializers/Accessors */
  //@{
  
  /** \brief . */
  DelayedLinearOpWithSolveFactory(
    const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf
    );

  /** \brief . */
  RCP<LinearOpWithSolveFactoryBase<Scalar> > getUnderlyingLOWSF();

  /** \brief . */
  RCP<const LinearOpWithSolveFactoryBase<Scalar> > getUnderlyingLOWSF() const;

  //@}

  /** \name Overridden from Teuchos::Describable. */
  //@{

  /** \brief . */
  std::string description() const;

  //@}

  /** @name Overridden from ParameterListAcceptor (simple forwarding functions) */
  //@{

  /** \brief . */
  void setParameterList(RCP<ParameterList> const& paramList);
  /** \brief . */
  RCP<ParameterList> getNonconstParameterList();
  /** \brief . */
  RCP<ParameterList> unsetParameterList();
  /** \brief . */
  RCP<const ParameterList> getParameterList() const;
  /** \brief . */
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from LinearOpWithSolveFactoyBase */
  //@{
  
  /** \brief . */
  virtual bool acceptsPreconditionerFactory() const;

  /** \brief . */
  virtual void setPreconditionerFactory(
    const RCP<PreconditionerFactoryBase<Scalar> > &precFactory,
    const std::string &precFactoryName
    );

  /** \brief . */
  virtual RCP<PreconditionerFactoryBase<Scalar> >
  getPreconditionerFactory() const;

  /** \brief . */
  virtual void unsetPreconditionerFactory(
    RCP<PreconditionerFactoryBase<Scalar> > *precFactory,
    std::string *precFactoryName
    );

  /** \brief . */
  virtual bool isCompatible(
    const LinearOpSourceBase<Scalar> &fwdOpSrc
    ) const;

  /** \brief . */
  virtual RCP<LinearOpWithSolveBase<Scalar> > createOp() const;

  /** \brief . */
  virtual void initializeOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;

  /** \brief . */
  virtual void initializeAndReuseOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op
    ) const;

  /** \brief . */
  virtual void uninitializeOp(
    LinearOpWithSolveBase<Scalar> *Op,
    RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
    RCP<const PreconditionerBase<Scalar> > *prec,
    RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
    ESupportSolveUse *supportSolveUse
    ) const;
 
  /** \brief . */
  virtual bool supportsPreconditionerInputType(
    const EPreconditionerInputType precOpType
    ) const;

  /** \brief . */
  virtual void initializePreconditionedOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const RCP<const PreconditionerBase<Scalar> > &prec,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;

  /** \brief . */
  virtual void initializeApproxPreconditionedOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse
    ) const;

  //@}

protected:

  /** \brief Overridden from Teuchos::VerboseObjectBase */
  //@{

  /** \brief . */
  void informUpdatedVerbosityState() const;

  //@}

private:

  RCP<LinearOpWithSolveFactoryBase<Scalar> > lowsf_;

  // Not defined and not to be called
  DelayedLinearOpWithSolveFactory();

};


} // namespace Thyra


#endif // THYRA_DELAYED_LINEAR_OP_WITH_SOLVE_FACTORY_DECL_HPP
