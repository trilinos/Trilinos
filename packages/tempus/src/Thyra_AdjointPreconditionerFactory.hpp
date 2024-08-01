//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Thyra_AdjointPreconditionerFactory_hpp
#define Thyra_AdjointPreconditionerFactory_hpp

#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_AdjointPreconditioner.hpp"

namespace Thyra {

/** \brief Concrete <tt>PreconditionerFactoryBase</tt> subclass that
 * wraps a preconditioner in AdjointPreconditioner.
 */
template <class Scalar>
class AdjointPreconditionerFactory
  : virtual public PreconditionerFactoryBase<Scalar> {
 public:
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  AdjointPreconditionerFactory() {}

  void nonconstInitialize(
      const RCP<PreconditionerFactoryBase<Scalar> > &prec_fac)
  {
    validateInitialize(prec_fac);
    prec_fac_ = prec_fac;
  }

  void initialize(const RCP<const PreconditionerFactoryBase<Scalar> > &prec_fac)
  {
    validateInitialize(prec_fac);
    prec_fac_ = prec_fac;
  }

  RCP<PreconditionerFactoryBase<Scalar> > getNonconstPreconditionerFactory()
  {
    return prec_fac_.getNonconstObj();
  }

  RCP<const PreconditionerFactoryBase<Scalar> > getPreconditionerFactory() const
  {
    return prec_fac_.getConstObj();
  }

  void uninitialize() { prec_fac_.uninitialize(); }

  /** \name Overridden from Teuchos::Describable. */
  //@{

  std::string description() const
  {
    std::ostringstream oss;
    oss << this->Teuchos::Describable::description() << "{"
        << "prec_fac=";
    if (!is_null(prec_fac_.getConstObj()))
      oss << prec_fac_.getConstObj()->description();
    else
      oss << "NULL";
    oss << "}";
    return oss.str();
  }

  //@}

  /** @name Overridden from ParameterListAcceptor (simple forwarding functions)
   */
  //@{

  void setParameterList(RCP<ParameterList> const &paramList)
  {
    prec_fac_.getNonconstObj()->setParameterList(paramList);
  }

  RCP<ParameterList> getNonconstParameterList()
  {
    return prec_fac_.getNonconstObj()->getNonconstParameterList();
  }

  RCP<ParameterList> unsetParameterList()
  {
    return prec_fac_.getNonconstObj()->unsetParameterList();
  }

  RCP<const ParameterList> getParameterList() const
  {
    return prec_fac_.getConstObj()->getParameterList();
  }

  RCP<const ParameterList> getValidParameters() const
  {
    return prec_fac_.getConstObj()->getValidParameters();
  }

  //@}

  //@}

  /** @name Overridden from PreconditionerFactoryBase */
  //@{

  bool isCompatible(const LinearOpSourceBase<Scalar> &fwdOpSrc) const
  {
    return prec_fac_.getConstObj()->isCompatible(fwdOpSrc);
  }

  RCP<PreconditionerBase<Scalar> > createPrec() const
  {
    return nonconstAdjointPreconditioner(prec_fac_.getConstObj()->createPrec());
  }

  void initializePrec(
      const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
      PreconditionerBase<Scalar> *precOp,
      const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED) const
  {
    using Teuchos::dyn_cast;
    using Teuchos::rcp_dynamic_cast;

    typedef DefaultScaledAdjointLinearOp<Scalar> ALO;
    typedef AdjointPreconditioner<Scalar> AP;
    const RCP<const ALO> alo =
        rcp_dynamic_cast<const ALO>(fwdOpSrc->getOp().assert_not_null());
    AP &ap = dyn_cast<AP>(*precOp);
    prec_fac_.getConstObj()->initializePrec(
        defaultLinearOpSource<Scalar>(alo->getOp()),
        ap.getNonconstPreconditioner().get(), supportSolveUse);
  }

  void uninitializePrec(PreconditionerBase<Scalar> *precOp,
                        RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc = NULL,
                        ESupportSolveUse *supportSolveUse                = NULL) const
  {
    using Teuchos::dyn_cast;

#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(0 == precOp);
#endif
    typedef AdjointPreconditioner<Scalar> AP;
    AP &ap = dyn_cast<AP>(*precOp);
    RCP<const LinearOpSourceBase<Scalar> > inner_fwdOpSrc;
    prec_fac_.getConstObj()->uninitializePrec(
        ap.getNonconstPreconditioner().get(), fwdOpSrc ? &inner_fwdOpSrc : NULL,
        supportSolveUse);
    if (fwdOpSrc)
      *fwdOpSrc =
          defaultLinearOpSource<Scalar>(adjoint(inner_fwdOpSrc->getOp()));
  }

  //@}

 private:
  // //////////////////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<
      PreconditionerFactoryBase<Scalar> >
      CNPFB;

  // //////////////////////////////
  // Private data members

  CNPFB prec_fac_;

  // //////////////////////////////
  // Private member functions

  static void validateInitialize(
      const RCP<const PreconditionerFactoryBase<Scalar> > &prec_fac)
  {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(is_null(prec_fac));
#else
    (void)prec_fac;
#endif
  }
};

/** \brief Nonmember constructor function.
 *
 * \relates AdjointPreconditionerFactory
 */
template <class Scalar>
RCP<AdjointPreconditionerFactory<Scalar> > adjointPreconditionerFactory()
{
  return Teuchos::rcp(new AdjointPreconditionerFactory<Scalar>());
}

/** \brief Nonmember constructor function.
 *
 * \relates AdjointPreconditionerFactory
 */
template <class Scalar>
RCP<AdjointPreconditionerFactory<Scalar> > nonconstAdjointPreconditionerFactory(
    const RCP<PreconditionerFactoryBase<Scalar> > &prec_fac)
{
  RCP<AdjointPreconditionerFactory<Scalar> > afac =
      Teuchos::rcp(new AdjointPreconditionerFactory<Scalar>());
  afac->nonconstInitialize(prec_fac);
  return afac;
}

/** \brief Nonmember constructor function.
 *
 * \relates AdjointPreconditionerFactory
 */
template <class Scalar>
RCP<AdjointPreconditionerFactory<Scalar> > adjointPreconditionerFactory(
    const RCP<const PreconditionerFactoryBase<Scalar> > &prec_fac)
{
  RCP<AdjointPreconditionerFactory<Scalar> > afac =
      Teuchos::rcp(new AdjointPreconditionerFactory<Scalar>());
  afac->initialize(prec_fac);
  return afac;
}

}  // end namespace Thyra

#endif
