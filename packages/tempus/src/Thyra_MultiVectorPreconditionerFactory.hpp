//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Thyra_MultiVectorPreconditionerFactory_hpp
#define Thyra_MultiVectorPreconditionerFactory_hpp

#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Thyra_MultiVectorLinearOp.hpp"
#include "Thyra_MultiVectorPreconditioner.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"

namespace Thyra {

/** \brief Concrete <tt>PreconditionerFactoryBase</tt> subclass that
 * wraps a preconditioner in MultiVectorPreconditioner.
 */
template <class Scalar>
class MultiVectorPreconditionerFactory
  : virtual public PreconditionerFactoryBase<Scalar> {
 public:
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  MultiVectorPreconditionerFactory() {}

  void nonconstInitialize(
      const RCP<PreconditionerFactoryBase<Scalar> > &prec_fac,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecRange,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecDomain)
  {
    validateInitialize(prec_fac, multiVecRange, multiVecDomain);
    prec_fac_       = prec_fac;
    multiVecRange_  = multiVecRange;
    multiVecDomain_ = multiVecDomain;
  }

  void initialize(const RCP<const PreconditionerFactoryBase<Scalar> > &prec_fac,
                  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
                      &multiVecRange,
                  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
                      &multiVecDomain)
  {
    validateInitialize(prec_fac, multiVecRange, multiVecDomain);
    prec_fac_       = prec_fac;
    multiVecRange_  = multiVecRange;
    multiVecDomain_ = multiVecDomain;
  }

  RCP<PreconditionerFactoryBase<Scalar> > getNonconstPreconditionerFactory()
  {
    return prec_fac_.getNonconstObj();
  }

  RCP<const PreconditionerFactoryBase<Scalar> > getPreconditionerFactory() const
  {
    return prec_fac_.getConstObj();
  }

  void uninitialize()
  {
    prec_fac_.uninitialize();
    multiVecRange_  = Teuchos::null;
    multiVecDomain_ = Teuchos::null;
  }

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
    return nonconstMultiVectorPreconditioner(
        prec_fac_.getConstObj()->createPrec(), multiVecRange_, multiVecDomain_);
  }

  void initializePrec(
      const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
      PreconditionerBase<Scalar> *precOp,
      const ESupportSolveUse supportSolveUse = SUPPORT_SOLVE_UNSPECIFIED) const
  {
    using Teuchos::dyn_cast;
    using Teuchos::rcp_dynamic_cast;

    typedef MultiVectorLinearOp<Scalar> MVLO;
    typedef MultiVectorPreconditioner<Scalar> MVP;
    const RCP<const MVLO> mvlo =
        rcp_dynamic_cast<const MVLO>(fwdOpSrc->getOp().assert_not_null());
    MVP &mvp = dyn_cast<MVP>(*precOp);
    prec_fac_.getConstObj()->initializePrec(
        defaultLinearOpSource<Scalar>(mvlo->getLinearOp()),
        mvp.getNonconstPreconditioner().get(), supportSolveUse);
  }

  void uninitializePrec(PreconditionerBase<Scalar> *precOp,
                        RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc = NULL,
                        ESupportSolveUse *supportSolveUse                = NULL) const
  {
    using Teuchos::dyn_cast;

#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(0 == precOp);
#endif
    typedef MultiVectorPreconditioner<Scalar> MVP;
    MVP &mvp = dyn_cast<MVP>(*precOp);
    RCP<const LinearOpSourceBase<Scalar> > inner_fwdOpSrc;
    prec_fac_.getConstObj()->uninitializePrec(
        mvp.getNonconstPreconditioner().get(),
        fwdOpSrc ? &inner_fwdOpSrc : NULL, supportSolveUse);
    if (fwdOpSrc)
      *fwdOpSrc = defaultLinearOpSource<Scalar>(multiVectorLinearOp(
          inner_fwdOpSrc->getOp(), multiVecRange_, multiVecDomain_));
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
  RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > multiVecRange_;
  RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > multiVecDomain_;

  // //////////////////////////////
  // Private member functions

  static void validateInitialize(
      const RCP<const PreconditionerFactoryBase<Scalar> > &prec_fac,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecRange,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecDomain)
  {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(is_null(prec_fac));
    TEUCHOS_TEST_FOR_EXCEPT(is_null(multiVecRange));
    TEUCHOS_TEST_FOR_EXCEPT(is_null(multiVecDomain));
    TEUCHOS_TEST_FOR_EXCEPT(multiVecRange->numBlocks() !=
                            multiVecDomain->numBlocks());
#else
    (void)prec_fac;
    (void)multiVecRange;
    (void)multiVecDomain;
#endif
  }
};

/** \brief Nonmember constructor function.
 *
 * \relates MultiVectorPreconditionerFactory
 */
template <class Scalar>
RCP<MultiVectorPreconditionerFactory<Scalar> >
multiVectorPreconditionerFactory()
{
  return Teuchos::rcp(new MultiVectorPreconditionerFactory<Scalar>());
}

/** \brief Nonmember constructor function.
 *
 * \relates MultiVectorPreconditionerFactory
 */
template <class Scalar>
RCP<MultiVectorPreconditionerFactory<Scalar> >
nonconstMultiVectorPreconditionerFactory(
    const RCP<PreconditionerFactoryBase<Scalar> > &prec_fac,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecRange,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecDomain)
{
  RCP<MultiVectorPreconditionerFactory<Scalar> > mvfac =
      Teuchos::rcp(new MultiVectorPreconditionerFactory<Scalar>());
  mvfac->nonconstInitialize(prec_fac, multiVecRange, multiVecDomain);
  return mvfac;
}

/** \brief Nonmember constructor function.
 *
 * \relates MultiVectorPreconditionerFactory
 */
template <class Scalar>
RCP<MultiVectorPreconditionerFactory<Scalar> > multiVectorPreconditionerFactory(
    const RCP<const PreconditionerFactoryBase<Scalar> > &prec_fac,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecRange,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecDomain)
{
  RCP<MultiVectorPreconditionerFactory<Scalar> > mvfac =
      Teuchos::rcp(new MultiVectorPreconditionerFactory<Scalar>());
  mvfac->initialize(prec_fac, multiVecRange, multiVecDomain);
  return mvfac;
}

}  // end namespace Thyra

#endif
