//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Thyra_ScaledIdentityLinearOpWithSolveFactory_hpp
#define Thyra_ScaledIdentityLinearOpWithSolveFactory_hpp

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_MultiVectorLinearOp.hpp"
#include "Thyra_MultiVectorPreconditioner.hpp"
#include "Thyra_MultiVectorPreconditionerFactory.hpp"
#include "Thyra_ScaledIdentityLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"

namespace Thyra {

/** \brief Create a LinearOpWithSolveFactory for a flattened-out multi-vector.
 */
template <class Scalar>
class ScaledIdentityLinearOpWithSolveFactory
  : virtual public LinearOpWithSolveFactoryBase<Scalar> {
 public:
  /** @name Overridden from Constructors/Initializers/Accessors */
  //@{

  /** \brief Construct to uninitialized. */
  ScaledIdentityLinearOpWithSolveFactory() {}

  /** \brief Initialize
   */
  void initialize(const RCP<const VectorSpaceBase<Scalar> > &space,
                  const Scalar &s)
  {
    space_ = space;
    s_     = s;
  }

  //@}

  /** @name Overridden from ParameterListAcceptor */
  //@{

  void setParameterList(RCP<ParameterList> const & /* paramList */) {}
  RCP<ParameterList> getNonconstParameterList()
  {
    return Teuchos::parameterList();
  }
  RCP<ParameterList> unsetParameterList() { return Teuchos::parameterList(); }
  RCP<const ParameterList> getParameterList() const
  {
    return Teuchos::parameterList();
  }
  RCP<const ParameterList> getValidParameters() const
  {
    return Teuchos::parameterList();
  }

  //@}

  /** \name Overridden from LinearOpWithSolveFactoyBase */
  //@{

  /** \brief returns false. */
  virtual bool acceptsPreconditionerFactory() const { return false; }

  /** \brief Throws exception. */
  virtual void setPreconditionerFactory(
      const RCP<PreconditionerFactoryBase<Scalar> > & /* precFactory */,
      const std::string & /* precFactoryName */
  )
  {
  }

  /** \brief Returns null . */
  virtual RCP<PreconditionerFactoryBase<Scalar> > getPreconditionerFactory()
      const
  {
    return Teuchos::null;
  }

  /** \brief Throws exception. */
  virtual void unsetPreconditionerFactory(
      RCP<PreconditionerFactoryBase<Scalar> > * /* precFactory */,
      std::string * /* precFactoryName */
  )
  {
  }

  virtual bool isCompatible(const LinearOpSourceBase<Scalar> &fwdOpSrc) const
  {
    return !is_null(
        Teuchos::rcp_dynamic_cast<
            const ScaledIdentityLinearOpWithSolve<Scalar> >(fwdOpSrc.getOp()));
  }

  virtual RCP<LinearOpWithSolveBase<Scalar> > createOp() const
  {
    return scaledIdentity(space_, s_);
  }

  virtual void initializeOp(
      const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
      LinearOpWithSolveBase<Scalar> *Op,
      const ESupportSolveUse supportSolveUse) const;

  virtual void initializeAndReuseOp(
      const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
      LinearOpWithSolveBase<Scalar> *Op) const
  {
    initializeOp(fwdOpSrc, Op, SUPPORT_SOLVE_UNSPECIFIED);
  }

  virtual void uninitializeOp(
      LinearOpWithSolveBase<Scalar> *Op,
      RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
      RCP<const PreconditionerBase<Scalar> > *prec,
      RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
      ESupportSolveUse *supportSolveUse) const;

  virtual bool supportsPreconditionerInputType(
      const EPreconditionerInputType /* precOpType */
  ) const
  {
    return false;
  }

  virtual void initializePreconditionedOp(
      const RCP<const LinearOpSourceBase<Scalar> > & /* fwdOpSrc */,
      const RCP<const PreconditionerBase<Scalar> > & /* prec */,
      LinearOpWithSolveBase<Scalar> * /* Op */,
      const ESupportSolveUse /* supportSolveUse */
  ) const
  {
  }

  virtual void initializeApproxPreconditionedOp(
      const RCP<const LinearOpSourceBase<Scalar> > & /* fwdOpSrc */,
      const RCP<const LinearOpSourceBase<Scalar> > & /* approxFwdOpSrc */,
      LinearOpWithSolveBase<Scalar> * /* Op */,
      const ESupportSolveUse /* supportSolveUse */
  ) const
  {
  }

  //@}

 protected:
  /** \brief Overridden from Teuchos::VerboseObjectBase */
  //@{

  void informUpdatedVerbosityState() const {}

  //@}

 private:
  RCP<const VectorSpaceBase<Scalar> > space_;
  Scalar s_;
};

/** \brief Nonmember constructor.
 *
 * \relates ScaledIdentityLinearOpWithSolveFactory
 */
template <class Scalar>
RCP<ScaledIdentityLinearOpWithSolveFactory<Scalar> > scaledIdentitySolveFactory(
    const RCP<const VectorSpaceBase<Scalar> > &space, const Scalar &s)
{
  RCP<ScaledIdentityLinearOpWithSolveFactory<Scalar> > lowsf =
      Teuchos::rcp(new ScaledIdentityLinearOpWithSolveFactory<Scalar>);
  lowsf->initialize(space, s);
  return lowsf;
}

// Overridden from LinearOpWithSolveFactoyBase

template <class Scalar>
void ScaledIdentityLinearOpWithSolveFactory<Scalar>::initializeOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse /* supportSolveUse */
) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(0 == Op);
#endif

  const RCP<const LinearOpBase<Scalar> > tmpFwdOp = fwdOpSrc->getOp();
  RCP<const LinearOpBase<Scalar> > fwdOp;
  Scalar fwdOp_scalar = 0.0;
  EOpTransp fwdOp_transp;
  unwrap<Scalar>(tmpFwdOp, &fwdOp_scalar, &fwdOp_transp, &fwdOp);

  const RCP<const ScaledIdentityLinearOpWithSolve<Scalar> > fwdSi =
      rcp_dynamic_cast<const ScaledIdentityLinearOpWithSolve<Scalar> >(fwdOp,
                                                                       true);

  dyn_cast<ScaledIdentityLinearOpWithSolve<Scalar> >(*Op).initialize(
      fwdSi->space(), fwdSi->scale());
}

template <class Scalar>
void ScaledIdentityLinearOpWithSolveFactory<Scalar>::uninitializeOp(
    LinearOpWithSolveBase<Scalar> *Op,
    RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
    RCP<const PreconditionerBase<Scalar> > *prec,
    RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
    ESupportSolveUse * /* supportSolveUse */
) const
{
  using Teuchos::dyn_cast;
  using Teuchos::is_null;
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(0 == Op);
#else
  (void)Op;
#endif  // TEUCHOS_DEBUG
  if (fwdOpSrc) *fwdOpSrc = Teuchos::null;
  if (prec) *prec = Teuchos::null;
  if (approxFwdOpSrc) *approxFwdOpSrc = Teuchos::null;
}

}  // namespace Thyra

#endif
