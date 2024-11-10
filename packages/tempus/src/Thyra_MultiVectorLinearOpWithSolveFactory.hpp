//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Thyra_MultiVectorLinearOpWithSolveFactory_hpp
#define Thyra_MultiVectorLinearOpWithSolveFactory_hpp

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_MultiVectorLinearOp.hpp"
#include "Thyra_MultiVectorPreconditioner.hpp"
#include "Thyra_MultiVectorPreconditionerFactory.hpp"
#include "Thyra_DefaultMultiVectorLinearOpWithSolve.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"

namespace Thyra {

/** \brief Create a LinearOpWithSolveFactory for a flattened-out multi-vector.
 */
template <class Scalar>
class MultiVectorLinearOpWithSolveFactory
  : virtual public LinearOpWithSolveFactoryBase<Scalar> {
 public:
  /** @name Overridden from Constructors/Initializers/Accessors */
  //@{

  /** \brief Construct to uninitialized. */
  MultiVectorLinearOpWithSolveFactory() {}

  /** \brief Initialize given a single non-const LOWSFB object.
   *
   * \param lowsf [in,persisting] The LOWSFB object that will be used to
   * create the LOWSB object for the diagonal blocks.
   * \param multiVecRange [?] Description.
   * \param multiVecDomain [?] Description.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>!is_null(lowsf)</tt>
   * </ul>
   *
   */
  void nonconstInitialize(
      const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecRange,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecDomain);

  /** \brief Initialize given a single const LOWSFB object.
   *
   * \param lowsf [in,persisting] The LOWSFB object that will be used to
   * create the LOWSB object for the diagonal blocks.
   * \param multiVecRange [?] Description.
   * \param multiVecDomain [?] Description.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>!is_null(lowsf)</tt>
   * </ul>
   *
   */
  void initialize(const RCP<const LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
                  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
                      &multiVecRange,
                  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
                      &multiVecDomain);

  RCP<LinearOpWithSolveFactoryBase<Scalar> > getUnderlyingLOWSF();

  RCP<const LinearOpWithSolveFactoryBase<Scalar> > getUnderlyingLOWSF() const;

  //@}

  /** \name Overridden from Teuchos::Describable. */
  //@{

  std::string description() const;

  //@}

  /** @name Overridden from ParameterListAcceptor (simple forwarding functions)
   */
  //@{

  void setParameterList(RCP<ParameterList> const &paramList);
  RCP<ParameterList> getNonconstParameterList();
  RCP<ParameterList> unsetParameterList();
  RCP<const ParameterList> getParameterList() const;
  RCP<const ParameterList> getValidParameters() const;

  //@}

  /** \name Overridden from LinearOpWithSolveFactoyBase */
  //@{

  /** \brief returns false. */
  virtual bool acceptsPreconditionerFactory() const;

  /** \brief Throws exception. */
  virtual void setPreconditionerFactory(
      const RCP<PreconditionerFactoryBase<Scalar> > &precFactory,
      const std::string &precFactoryName);

  /** \brief Returns null . */
  virtual RCP<PreconditionerFactoryBase<Scalar> > getPreconditionerFactory()
      const;

  /** \brief Throws exception. */
  virtual void unsetPreconditionerFactory(
      RCP<PreconditionerFactoryBase<Scalar> > *precFactory,
      std::string *precFactoryName);

  virtual bool isCompatible(const LinearOpSourceBase<Scalar> &fwdOpSrc) const;

  virtual RCP<LinearOpWithSolveBase<Scalar> > createOp() const;

  virtual void initializeOp(
      const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
      LinearOpWithSolveBase<Scalar> *Op,
      const ESupportSolveUse supportSolveUse) const;

  virtual void initializeAndReuseOp(
      const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
      LinearOpWithSolveBase<Scalar> *Op) const;

  virtual void uninitializeOp(
      LinearOpWithSolveBase<Scalar> *Op,
      RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
      RCP<const PreconditionerBase<Scalar> > *prec,
      RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
      ESupportSolveUse *supportSolveUse) const;

  virtual bool supportsPreconditionerInputType(
      const EPreconditionerInputType precOpType) const;

  virtual void initializePreconditionedOp(
      const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
      const RCP<const PreconditionerBase<Scalar> > &prec,
      LinearOpWithSolveBase<Scalar> *Op,
      const ESupportSolveUse supportSolveUse) const;

  virtual void initializeApproxPreconditionedOp(
      const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
      const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
      LinearOpWithSolveBase<Scalar> *Op,
      const ESupportSolveUse supportSolveUse) const;

  //@}

 protected:
  /** \brief Overridden from Teuchos::VerboseObjectBase */
  //@{

  void informUpdatedVerbosityState() const;

  //@}

 private:
  typedef Teuchos::ConstNonconstObjectContainer<
      LinearOpWithSolveFactoryBase<Scalar> >
      LOWSF_t;

  LOWSF_t lowsf_;
  RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > multiVecRange_;
  RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > multiVecDomain_;
};

/** \brief Nonmember constructor.
 *
 * \relates MultiVectorLinearOpWithSolveFactory
 */
template <class Scalar>
RCP<MultiVectorLinearOpWithSolveFactory<Scalar> >
nonconstMultiVectorLinearOpWithSolveFactory(
    const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecRange,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecDomain)
{
  RCP<MultiVectorLinearOpWithSolveFactory<Scalar> > mvlowsf =
      Teuchos::rcp(new MultiVectorLinearOpWithSolveFactory<Scalar>);
  mvlowsf->nonconstInitialize(lowsf, multiVecRange, multiVecDomain);
  return mvlowsf;
}

/** \brief Nonmember constructor.
 *
 * \relates MultiVectorLinearOpWithSolveFactory
 */
template <class Scalar>
RCP<MultiVectorLinearOpWithSolveFactory<Scalar> >
nonconstMultiVectorLinearOpWithSolveFactory(
    const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
    const int num_blocks)
{
  RCP<LinearOpWithSolveBase<Scalar> > op = lowsf->createOp();
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > mv_domain =
      Thyra::multiVectorProductVectorSpace(op->domain(), num_blocks);
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > mv_range =
      Thyra::multiVectorProductVectorSpace(op->range(), num_blocks);
  return nonconstMultiVectorLinearOpWithSolveFactory(lowsf, mv_range,
                                                     mv_domain);
}

/** \brief Nonmember constructor.
 *
 * \relates MultiVectorLinearOpWithSolveFactory
 */
template <class Scalar>
RCP<MultiVectorLinearOpWithSolveFactory<Scalar> >
multiVectorLinearOpWithSolveFactory(
    const RCP<const LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecRange,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecDomain)
{
  RCP<MultiVectorLinearOpWithSolveFactory<Scalar> > mvlowsf =
      Teuchos::rcp(new MultiVectorLinearOpWithSolveFactory<Scalar>);
  mvlowsf->initialize(lowsf, multiVecRange, multiVecDomain);
  return mvlowsf;
}

/** \brief Nonmember constructor.
 *
 * \relates MultiVectorLinearOpWithSolveFactory
 */
template <class Scalar>
RCP<MultiVectorLinearOpWithSolveFactory<Scalar> >
multiVectorLinearOpWithSolveFactory(
    const RCP<const LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
    const int num_blocks)
{
  RCP<LinearOpWithSolveBase<Scalar> > op = lowsf->createOp();
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > mv_domain =
      Thyra::multiVectorProductVectorSpace(op->domain(), num_blocks);
  RCP<const Thyra::DefaultMultiVectorProductVectorSpace<Scalar> > mv_range =
      Thyra::multiVectorProductVectorSpace(op->range(), num_blocks);
  return multiVectorLinearOpWithSolveFactory(lowsf, mv_range, mv_domain);
}

// Overridden from Constructors/Initializers/Accessors

template <class Scalar>
void MultiVectorLinearOpWithSolveFactory<Scalar>::nonconstInitialize(
    const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecRange,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecDomain)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(is_null(lowsf));
#endif
  lowsf_.initialize(lowsf);
  multiVecRange_  = multiVecRange;
  multiVecDomain_ = multiVecDomain;
}

template <class Scalar>
void MultiVectorLinearOpWithSolveFactory<Scalar>::initialize(
    const RCP<const LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecRange,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecDomain)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(is_null(lowsf));
#endif
  lowsf_.initialize(lowsf);
  multiVecRange_  = multiVecRange;
  multiVecDomain_ = multiVecDomain;
}

template <class Scalar>
RCP<LinearOpWithSolveFactoryBase<Scalar> >
MultiVectorLinearOpWithSolveFactory<Scalar>::getUnderlyingLOWSF()
{
  return lowsf_.getNonconstObj();
}

template <class Scalar>
RCP<const LinearOpWithSolveFactoryBase<Scalar> >
MultiVectorLinearOpWithSolveFactory<Scalar>::getUnderlyingLOWSF() const
{
  return lowsf_.getConstObj();
}

// Overridden from Teuchos::Describable

template <class Scalar>
std::string MultiVectorLinearOpWithSolveFactory<Scalar>::description() const
{
  std::ostringstream oss;
  oss << this->Teuchos::Describable::description() << "{"
      << "lowsf=";
  if (!is_null(lowsf_.getConstObj()))
    oss << lowsf_.getConstObj()->description();
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}

// Overridden from ParameterListAcceptor

template <class Scalar>
void MultiVectorLinearOpWithSolveFactory<Scalar>::setParameterList(
    RCP<ParameterList> const &paramList)
{
  lowsf_.getNonconstObj()->setParameterList(paramList);
}

template <class Scalar>
RCP<ParameterList>
MultiVectorLinearOpWithSolveFactory<Scalar>::getNonconstParameterList()
{
  return lowsf_.getNonconstObj()->getNonconstParameterList();
}

template <class Scalar>
RCP<ParameterList>
MultiVectorLinearOpWithSolveFactory<Scalar>::unsetParameterList()
{
  return lowsf_.getNonconstObj()->unsetParameterList();
}

template <class Scalar>
RCP<const ParameterList>
MultiVectorLinearOpWithSolveFactory<Scalar>::getParameterList() const
{
  return lowsf_.getConstObj()->getParameterList();
}

template <class Scalar>
RCP<const ParameterList>
MultiVectorLinearOpWithSolveFactory<Scalar>::getValidParameters() const
{
  return lowsf_.getConstObj()->getValidParameters();
}

// Overridden from LinearOpWithSolveFactoyBase

template <class Scalar>
bool MultiVectorLinearOpWithSolveFactory<Scalar>::acceptsPreconditionerFactory()
    const
{
  return lowsf_.getConstObj()->acceptsPreconditionerFactory();
}

template <class Scalar>
void MultiVectorLinearOpWithSolveFactory<Scalar>::setPreconditionerFactory(
    const RCP<PreconditionerFactoryBase<Scalar> > &precFactory,
    const std::string &precFactoryName)
{
  typedef MultiVectorPreconditionerFactory<Scalar> MVPF;
  RCP<MVPF> mvpf = Teuchos::rcp_dynamic_cast<MVPF>(precFactory);
  lowsf_.getNonconstObj()->setPreconditionerFactory(
      mvpf->getNonconstPreconditionerFactory(), precFactoryName);
}

template <class Scalar>
RCP<PreconditionerFactoryBase<Scalar> >
MultiVectorLinearOpWithSolveFactory<Scalar>::getPreconditionerFactory() const
{
  RCP<PreconditionerFactoryBase<Scalar> > prec_fac =
      lowsf_.getConstObj()->getPreconditionerFactory();
  if (prec_fac == Teuchos::null)
    return Teuchos::null;
  else
    return nonconstMultiVectorPreconditionerFactory(prec_fac, multiVecRange_,
                                                    multiVecDomain_);
}

template <class Scalar>
void MultiVectorLinearOpWithSolveFactory<Scalar>::unsetPreconditionerFactory(
    RCP<PreconditionerFactoryBase<Scalar> > *precFactory,
    std::string *precFactoryName)
{
  RCP<PreconditionerFactoryBase<Scalar> > inner_precFactory;
  lowsf_.getNonconstObj()->unsetPreconditionerFactory(
      precFactory ? &inner_precFactory : NULL, precFactoryName);
  if (precFactory)
    *precFactory = nonconstMultiVectorPreconditionerFactory(
        inner_precFactory, multiVecRange_, multiVecDomain_);
}

template <class Scalar>
bool MultiVectorLinearOpWithSolveFactory<Scalar>::isCompatible(
    const LinearOpSourceBase<Scalar> &fwdOpSrc) const
{
  typedef MultiVectorLinearOp<Scalar> MVLO;
  RCP<const MVLO> mvlo =
      Teuchos::rcp_dynamic_cast<const MVLO>(fwdOpSrc.getOp().assert_not_null());
  RCP<const LinearOpSourceBase<Scalar> > inner_fwdOpSrc =
      defaultLinearOpSource<Scalar>(mvlo->getLinearOp());
  return lowsf_.getConstObj()->isCompatible(*inner_fwdOpSrc);
}

template <class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
MultiVectorLinearOpWithSolveFactory<Scalar>::createOp() const
{
  return nonconstMultiVectorLinearOpWithSolve<Scalar>(
      lowsf_.getConstObj()->createOp(), multiVecRange_, multiVecDomain_);
}

template <class Scalar>
void MultiVectorLinearOpWithSolveFactory<Scalar>::initializeOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(0 == Op);
#endif

  // Set the verbosity settings for the wrapped LOWSF object!
  lowsf_.getConstObj()->setOStream(this->getOStream());
  lowsf_.getConstObj()->setVerbLevel(this->getVerbLevel());

  typedef MultiVectorLinearOp<Scalar> MVLO;
  typedef DefaultMultiVectorLinearOpWithSolve<Scalar> MVLOWS;
  const RCP<const MVLO> mvlo =
      rcp_dynamic_cast<const MVLO>(fwdOpSrc->getOp().assert_not_null());
  MVLOWS &mvlows = dyn_cast<MVLOWS>(*Op);

  lowsf_.getConstObj()->initializeOp(
      defaultLinearOpSource<Scalar>(mvlo->getLinearOp()),
      mvlows.getNonconstLinearOpWithSolve().get(), supportSolveUse);
}

template <class Scalar>
void MultiVectorLinearOpWithSolveFactory<Scalar>::initializeAndReuseOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(0 == Op);
#endif

  // Set the verbosity settings for the wrapped LOWSF object!
  lowsf_.getConstObj()->setOStream(this->getOStream());
  lowsf_.getConstObj()->setVerbLevel(this->getVerbLevel());

  typedef MultiVectorLinearOp<Scalar> MVLO;
  typedef DefaultMultiVectorLinearOpWithSolve<Scalar> MVLOWS;
  const RCP<const MVLO> mvlo =
      rcp_dynamic_cast<const MVLO>(fwdOpSrc->getOp().assert_not_null());
  MVLOWS &mvlows = dyn_cast<MVLOWS>(*Op);

  lowsf_.getConstObj()->initializeAndReuseOp(
      defaultLinearOpSource<Scalar>(mvlo->getLinearOp()),
      mvlows.getNonconstLinearOpWithSolve().get());
}

template <class Scalar>
void MultiVectorLinearOpWithSolveFactory<Scalar>::uninitializeOp(
    LinearOpWithSolveBase<Scalar> *Op,
    RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
    RCP<const PreconditionerBase<Scalar> > *prec,
    RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
    ESupportSolveUse *supportSolveUse) const
{
  using Teuchos::dyn_cast;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(0 == Op);
#endif
  typedef DefaultMultiVectorLinearOpWithSolve<Scalar> MVLOWS;
  MVLOWS &mvlowsOp = dyn_cast<MVLOWS>(*Op);
  RCP<const LinearOpSourceBase<Scalar> > inner_fwdOpSrc;
  RCP<const PreconditionerBase<Scalar> > inner_prec;
  RCP<const LinearOpSourceBase<Scalar> > inner_approxFwdOpSrc;
  lowsf_.getConstObj()->uninitializeOp(
      mvlowsOp.getNonconstLinearOpWithSolve().get(),
      fwdOpSrc ? &inner_fwdOpSrc : NULL, prec ? &inner_prec : NULL,
      approxFwdOpSrc ? &inner_approxFwdOpSrc : NULL, supportSolveUse);
  if (fwdOpSrc)
    *fwdOpSrc = defaultLinearOpSource<Scalar>(multiVectorLinearOp(
        inner_fwdOpSrc->getOp(), multiVecRange_, multiVecDomain_));
  if (prec)
    *prec =
        multiVectorPreconditioner(inner_prec, multiVecRange_, multiVecDomain_);
  if (fwdOpSrc)
    *approxFwdOpSrc = defaultLinearOpSource<Scalar>(multiVectorLinearOp(
        inner_approxFwdOpSrc->getOp(), multiVecRange_, multiVecDomain_));
}

template <class Scalar>
bool MultiVectorLinearOpWithSolveFactory<
    Scalar>::supportsPreconditionerInputType(const EPreconditionerInputType
                                                 precOpType) const
{
  return lowsf_.getConstObj()->supportsPreconditionerInputType(precOpType);
}

template <class Scalar>
void MultiVectorLinearOpWithSolveFactory<Scalar>::initializePreconditionedOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const RCP<const PreconditionerBase<Scalar> > &prec,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(0 == Op);
#endif

  // Set the verbosity settings for the wrapped LOWSF object!
  lowsf_.getConstObj()->setOStream(this->getOStream());
  lowsf_.getConstObj()->setVerbLevel(this->getVerbLevel());

  typedef MultiVectorLinearOp<Scalar> MVLO;
  typedef MultiVectorPreconditioner<Scalar> MVP;
  typedef DefaultMultiVectorLinearOpWithSolve<Scalar> MVLOWS;
  const RCP<const MVLO> mvlo =
      rcp_dynamic_cast<const MVLO>(fwdOpSrc->getOp().assert_not_null());
  const RCP<const MVP> mvp = rcp_dynamic_cast<const MVP>(prec);
  MVLOWS &mvlows           = dyn_cast<MVLOWS>(*Op);

  lowsf_.getConstObj()->initializePreconditionedOp(
      defaultLinearOpSource<Scalar>(mvlo->getLinearOp()),
      mvp->getPreconditioner(), mvlows.getNonconstLinearOpWithSolve().get(),
      supportSolveUse);
}

template <class Scalar>
void MultiVectorLinearOpWithSolveFactory<Scalar>::
    initializeApproxPreconditionedOp(
        const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
        const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
        LinearOpWithSolveBase<Scalar> *Op,
        const ESupportSolveUse supportSolveUse) const
{
  using Teuchos::dyn_cast;
  using Teuchos::rcp_dynamic_cast;

#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(0 == Op);
#endif

  // Set the verbosity settings for the wrapped LOWSF object!
  lowsf_.getConstObj()->setOStream(this->getOStream());
  lowsf_.getConstObj()->setVerbLevel(this->getVerbLevel());

  typedef MultiVectorLinearOp<Scalar> MVLO;
  typedef DefaultMultiVectorLinearOpWithSolve<Scalar> MVLOWS;
  const RCP<const MVLO> mvlo =
      rcp_dynamic_cast<const MVLO>(fwdOpSrc->getOp().assert_not_null());
  const RCP<const MVLO> amvlo =
      rcp_dynamic_cast<const MVLO>(approxFwdOpSrc->getOp().assert_not_null());
  MVLOWS &mvlows = dyn_cast<MVLOWS>(*Op);

  lowsf_.getConstObj()->initializeApproxPreconditionedOp(
      defaultLinearOpSource<Scalar>(mvlo->getLinearOp()),
      defaultLinearOpSource<Scalar>(amvlo->getLinearOp()),
      mvlows.getNonconstLinearOpWithSolve().get(), supportSolveUse);
}

// protected

template <class Scalar>
void MultiVectorLinearOpWithSolveFactory<Scalar>::informUpdatedVerbosityState()
    const
{
  lowsf_.getConstObj()->setVerbLevel(this->getVerbLevel());
  lowsf_.getConstObj()->setOStream(this->getOStream());
}

}  // namespace Thyra

#endif
