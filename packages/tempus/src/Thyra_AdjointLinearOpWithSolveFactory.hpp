//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Thyra_AdjointLinearOpWithSolveFactory_hpp
#define Thyra_AdjointLinearOpWithSolveFactory_hpp

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultAdjointLinearOpWithSolve.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_AdjointPreconditioner.hpp"
#include "Thyra_AdjointPreconditionerFactory.hpp"

namespace Thyra {

/** \brief Create a LinearOpWithSolveFactory for an adjoint linear op.
 */
template <class Scalar>
class AdjointLinearOpWithSolveFactory
  : virtual public LinearOpWithSolveFactoryBase<Scalar> {
 public:
  /** @name Overridden from Constructors/Initializers/Accessors */
  //@{

  /** \brief Construct to uninitialized. */
  AdjointLinearOpWithSolveFactory() {}

  /** \brief Initialize given a single non-const LOWSFB object.
   *
   * \param lowsf [in,persisting] The LOWSFB object that will be used to
   * create the LOWSB object for the original system.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>!is_null(lowsf)</tt>
   * </ul>
   *
   */
  void nonconstInitialize(
      const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf);

  /** \brief Initialize given a single const LOWSFB object.
   *
   * \param lowsf [in,persisting] The LOWSFB object that will be used to
   * create the LOWSB object for the original system.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>!is_null(lowsf)</tt>
   * </ul>
   *
   */
  void initialize(
      const RCP<const LinearOpWithSolveFactoryBase<Scalar> > &lowsf);

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
};

/** \brief Nonmember constructor.
 *
 * \relates AdjointLinearOpWithSolveFactory
 */
template <class Scalar>
RCP<const AdjointLinearOpWithSolveFactory<Scalar> >
adjointLinearOpWithSolveFactory(
    const RCP<const LinearOpWithSolveFactoryBase<Scalar> > &lowsf)
{
  RCP<AdjointLinearOpWithSolveFactory<Scalar> > alowsf =
      Teuchos::rcp(new AdjointLinearOpWithSolveFactory<Scalar>);
  alowsf->initialize(lowsf);
  return alowsf;
}

/** \brief Nonmember constructor.
 *
 * \relates AdjointLinearOpWithSolveFactory
 */
template <class Scalar>
RCP<AdjointLinearOpWithSolveFactory<Scalar> >
nonconstAdjointLinearOpWithSolveFactory(
    const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf)
{
  RCP<AdjointLinearOpWithSolveFactory<Scalar> > alowsf =
      Teuchos::rcp(new AdjointLinearOpWithSolveFactory<Scalar>);
  alowsf->nonconstInitialize(lowsf);
  return alowsf;
}

// Overridden from Constructors/Initializers/Accessors

template <class Scalar>
void AdjointLinearOpWithSolveFactory<Scalar>::nonconstInitialize(
    const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(is_null(lowsf));
#endif
  lowsf_.initialize(lowsf);
}

template <class Scalar>
void AdjointLinearOpWithSolveFactory<Scalar>::initialize(
    const RCP<const LinearOpWithSolveFactoryBase<Scalar> > &lowsf)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(is_null(lowsf));
#endif
  lowsf_.initialize(lowsf);
}

template <class Scalar>
RCP<LinearOpWithSolveFactoryBase<Scalar> >
AdjointLinearOpWithSolveFactory<Scalar>::getUnderlyingLOWSF()
{
  return lowsf_.getNonconstObj();
}

template <class Scalar>
RCP<const LinearOpWithSolveFactoryBase<Scalar> >
AdjointLinearOpWithSolveFactory<Scalar>::getUnderlyingLOWSF() const
{
  return lowsf_.getConstObj();
}

// Overridden from Teuchos::Describable

template <class Scalar>
std::string AdjointLinearOpWithSolveFactory<Scalar>::description() const
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
void AdjointLinearOpWithSolveFactory<Scalar>::setParameterList(
    RCP<ParameterList> const &paramList)
{
  lowsf_.getNonconstObj()->setParameterList(paramList);
}

template <class Scalar>
RCP<ParameterList>
AdjointLinearOpWithSolveFactory<Scalar>::getNonconstParameterList()
{
  return lowsf_.getNonconstObj()->getNonconstParameterList();
}

template <class Scalar>
RCP<ParameterList> AdjointLinearOpWithSolveFactory<Scalar>::unsetParameterList()
{
  return lowsf_.getNonconstObj()->unsetParameterList();
}

template <class Scalar>
RCP<const ParameterList>
AdjointLinearOpWithSolveFactory<Scalar>::getParameterList() const
{
  return lowsf_.getConstObj()->getParameterList();
}

template <class Scalar>
RCP<const ParameterList>
AdjointLinearOpWithSolveFactory<Scalar>::getValidParameters() const
{
  return lowsf_.getConstObj()->getValidParameters();
}

// Overridden from LinearOpWithSolveFactoyBase

template <class Scalar>
bool AdjointLinearOpWithSolveFactory<Scalar>::acceptsPreconditionerFactory()
    const
{
  return lowsf_.getConstObj()->acceptsPreconditionerFactory();
}

template <class Scalar>
void AdjointLinearOpWithSolveFactory<Scalar>::setPreconditionerFactory(
    const RCP<PreconditionerFactoryBase<Scalar> > &precFactory,
    const std::string &precFactoryName)
{
  typedef AdjointPreconditionerFactory<Scalar> APF;
  RCP<APF> apf = Teuchos::rcp_dynamic_cast<APF>(precFactory);
  lowsf_.getNonconstObj()->setPreconditionerFactory(
      apf->getNonconstPreconditionerFactory(), precFactoryName);
}

template <class Scalar>
RCP<PreconditionerFactoryBase<Scalar> >
AdjointLinearOpWithSolveFactory<Scalar>::getPreconditionerFactory() const
{
  RCP<PreconditionerFactoryBase<Scalar> > prec_fac =
      lowsf_.getConstObj()->getPreconditionerFactory();
  if (prec_fac == Teuchos::null)
    return Teuchos::null;
  else
    return nonconstAdjointPreconditionerFactory(prec_fac);
}

template <class Scalar>
void AdjointLinearOpWithSolveFactory<Scalar>::unsetPreconditionerFactory(
    RCP<PreconditionerFactoryBase<Scalar> > *precFactory,
    std::string *precFactoryName)
{
  RCP<PreconditionerFactoryBase<Scalar> > inner_precFactory;
  lowsf_.getNonconstObj()->unsetPreconditionerFactory(
      precFactory ? &inner_precFactory : NULL, precFactoryName);
  if (precFactory)
    *precFactory = nonconstAdjointPreconditionerFactory(inner_precFactory);
}

template <class Scalar>
bool AdjointLinearOpWithSolveFactory<Scalar>::isCompatible(
    const LinearOpSourceBase<Scalar> &fwdOpSrc) const
{
  typedef DefaultScaledAdjointLinearOp<Scalar> ALO;
  RCP<const ALO> alo =
      Teuchos::rcp_dynamic_cast<const ALO>(fwdOpSrc.getOp().assert_not_null());
  RCP<const LinearOpSourceBase<Scalar> > inner_fwdOpSrc =
      defaultLinearOpSource<Scalar>(alo->getOp());
  return lowsf_.getConstObj()->isCompatible(*inner_fwdOpSrc);
}

template <class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
AdjointLinearOpWithSolveFactory<Scalar>::createOp() const
{
  return nonconstAdjointLows<Scalar>(lowsf_.getConstObj()->createOp());
}

template <class Scalar>
void AdjointLinearOpWithSolveFactory<Scalar>::initializeOp(
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

  typedef DefaultScaledAdjointLinearOp<Scalar> ALO;
  typedef DefaultAdjointLinearOpWithSolve<Scalar> ALOWS;
  const RCP<const ALO> alo =
      rcp_dynamic_cast<const ALO>(fwdOpSrc->getOp().assert_not_null());
  ALOWS &alows = dyn_cast<ALOWS>(*Op);

  lowsf_.getConstObj()->initializeOp(
      defaultLinearOpSource<Scalar>(alo->getOrigOp()),
      alows.getNonconstOp().get(), supportSolveUse);
}

template <class Scalar>
void AdjointLinearOpWithSolveFactory<Scalar>::initializeAndReuseOp(
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

  typedef DefaultScaledAdjointLinearOp<Scalar> ALO;
  typedef DefaultAdjointLinearOpWithSolve<Scalar> ALOWS;
  const RCP<const ALO> alo =
      rcp_dynamic_cast<const ALO>(fwdOpSrc->getOp().assert_not_null());
  ALOWS &alows = dyn_cast<ALOWS>(*Op);

  lowsf_.getConstObj()->initializeAndReuseOp(
      defaultLinearOpSource<Scalar>(alo->getOrigOp()),
      alows.getNonconstOp().get());
}

template <class Scalar>
void AdjointLinearOpWithSolveFactory<Scalar>::uninitializeOp(
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
  typedef DefaultAdjointLinearOpWithSolve<Scalar> ALOWS;
  ALOWS &alowsOp = dyn_cast<ALOWS>(*Op);
  RCP<const LinearOpSourceBase<Scalar> > inner_fwdOpSrc;
  RCP<const PreconditionerBase<Scalar> > inner_prec;
  RCP<const LinearOpSourceBase<Scalar> > inner_approxFwdOpSrc;
  lowsf_.getConstObj()->uninitializeOp(
      alowsOp.getNonconstOp().get(), fwdOpSrc ? &inner_fwdOpSrc : NULL,
      prec ? &inner_prec : NULL, approxFwdOpSrc ? &inner_approxFwdOpSrc : NULL,
      supportSolveUse);
  if (fwdOpSrc)
    *fwdOpSrc = defaultLinearOpSource<Scalar>(adjoint(inner_fwdOpSrc->getOp()));
  if (prec) *prec = adjointPreconditioner(inner_prec);
  if (fwdOpSrc)
    *approxFwdOpSrc =
        defaultLinearOpSource<Scalar>(adjoint(inner_approxFwdOpSrc->getOp()));
}

template <class Scalar>
bool AdjointLinearOpWithSolveFactory<Scalar>::supportsPreconditionerInputType(
    const EPreconditionerInputType precOpType) const
{
  return lowsf_.getConstObj()->supportsPreconditionerInputType(precOpType);
}

template <class Scalar>
void AdjointLinearOpWithSolveFactory<Scalar>::initializePreconditionedOp(
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

  typedef DefaultScaledAdjointLinearOp<Scalar> ALO;
  typedef AdjointPreconditioner<Scalar> AP;
  typedef DefaultAdjointLinearOpWithSolve<Scalar> ALOWS;
  const RCP<const ALO> alo =
      rcp_dynamic_cast<const ALO>(fwdOpSrc->getOp().assert_not_null());
  const RCP<const AP> ap = rcp_dynamic_cast<const AP>(prec);
  ALOWS &alows           = dyn_cast<ALOWS>(*Op);

  lowsf_.getConstObj()->initializePreconditionedOp(
      defaultLinearOpSource<Scalar>(alo->getOp()), ap->getPreconditioner(),
      alows.getNonconstOp().get(), supportSolveUse);
}

template <class Scalar>
void AdjointLinearOpWithSolveFactory<Scalar>::initializeApproxPreconditionedOp(
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

  typedef DefaultAdjointLinearOpWithSolve<Scalar> ALO;
  typedef DefaultAdjointLinearOpWithSolve<Scalar> ALOWS;
  const RCP<const ALO> alo =
      rcp_dynamic_cast<const ALO>(fwdOpSrc->getOp().assert_not_null());
  const RCP<const ALO> aalo =
      rcp_dynamic_cast<const ALO>(approxFwdOpSrc->getOp().assert_not_null());
  ALOWS &alows = dyn_cast<ALOWS>(*Op);

  lowsf_.getConstObj()->initializeApproxPreconditionedOp(
      defaultLinearOpSource<Scalar>(alo->getOp()),
      defaultLinearOpSource<Scalar>(aalo->getOp()), alows.getNonconstOp().get(),
      supportSolveUse);
}

// protected

template <class Scalar>
void AdjointLinearOpWithSolveFactory<Scalar>::informUpdatedVerbosityState()
    const
{
  lowsf_.getConstObj()->setVerbLevel(this->getVerbLevel());
  lowsf_.getConstObj()->setOStream(this->getOStream());
}

}  // namespace Thyra

#endif
