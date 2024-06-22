//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Thyra_ReuseLinearOpWithSolveFactory_hpp
#define Thyra_ReuseLinearOpWithSolveFactory_hpp

#include "Thyra_LinearOpWithSolveFactoryBase.hpp"
#include "Thyra_ReusePreconditionerFactory.hpp"

namespace Thyra {

/** \brief A LinearOpWithSolveFactory that is designed to reuse an already
 * created/initialized preconditioner.
 */
template <class Scalar>
class ReuseLinearOpWithSolveFactory
  : virtual public LinearOpWithSolveFactoryBase<Scalar> {
 public:
  /** @name Overridden from Constructors/Initializers/Accessors */
  //@{

  /** \brief Construct to uninitialized. */
  ReuseLinearOpWithSolveFactory() {}

  /** \brief Initialize given a single non-const LOWSFB object.
   *
   * \param lowsf [in,persisting] The LOWSFB object that will be used to
   * create the LOWSB object for the diagonal blocks.
   * \param prec [?] Description.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>!is_null(lowsf)</tt>
   * </ul>
   *
   */
  void nonconstInitialize(
      const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
      const RCP<PreconditionerBase<Scalar> > &prec);

  /** \brief Initialize given a single const LOWSFB object.
   *
   * \param lowsf [in,persisting] The LOWSFB object that will be used to
   * create the LOWSB object for the diagonal blocks.
   * \param prec [?] Description.
   *
   * <b>Preconditions:</b><ul>
   * <li><tt>!is_null(lowsf)</tt>
   * </ul>
   *
   */
  void initialize(const RCP<const LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
                  const RCP<PreconditionerBase<Scalar> > &prec);

  RCP<LinearOpWithSolveFactoryBase<Scalar> > getUnderlyingLOWSF();

  RCP<const LinearOpWithSolveFactoryBase<Scalar> > getUnderlyingLOWSF() const;

  RCP<PreconditionerBase<Scalar> > getUnderlyingPreconditioner();

  RCP<const PreconditionerBase<Scalar> > getUnderlyingPreconditioner() const;

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
  RCP<PreconditionerBase<Scalar> > prec_;
};

/** \brief Nonmember constructor.
 *
 * \relates ReuseLinearOpWithSolveFactory
 */
template <class Scalar>
RCP<ReuseLinearOpWithSolveFactory<Scalar> >
nonconstReuseLinearOpWithSolveFactory(
    const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
    const RCP<PreconditionerBase<Scalar> > &prec)
{
  RCP<ReuseLinearOpWithSolveFactory<Scalar> > rlowsf =
      Teuchos::rcp(new ReuseLinearOpWithSolveFactory<Scalar>);
  rlowsf->nonconstInitialize(lowsf, prec);
  return rlowsf;
}

/** \brief Nonmember constructor.
 *
 * \relates ReuseLinearOpWithSolveFactory
 */
template <class Scalar>
RCP<ReuseLinearOpWithSolveFactory<Scalar> > reuseLinearOpWithSolveFactory(
    const RCP<const LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
    const RCP<PreconditionerBase<Scalar> > &prec)
{
  RCP<ReuseLinearOpWithSolveFactory<Scalar> > rlowsf =
      Teuchos::rcp(new ReuseLinearOpWithSolveFactory<Scalar>);
  rlowsf->initialize(lowsf, prec);
  return rlowsf;
}

// Overridden from Constructors/Initializers/Accessors

template <class Scalar>
void ReuseLinearOpWithSolveFactory<Scalar>::nonconstInitialize(
    const RCP<LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
    const RCP<PreconditionerBase<Scalar> > &prec)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(is_null(lowsf));
  TEUCHOS_TEST_FOR_EXCEPT(is_null(prec));
#endif
  lowsf_.initialize(lowsf);
  prec_ = prec;
}

template <class Scalar>
void ReuseLinearOpWithSolveFactory<Scalar>::initialize(
    const RCP<const LinearOpWithSolveFactoryBase<Scalar> > &lowsf,
    const RCP<PreconditionerBase<Scalar> > &prec)
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_TEST_FOR_EXCEPT(is_null(lowsf));
  TEUCHOS_TEST_FOR_EXCEPT(is_null(prec));
#endif
  lowsf_.initialize(lowsf);
  prec_ = prec;
}

template <class Scalar>
RCP<LinearOpWithSolveFactoryBase<Scalar> >
ReuseLinearOpWithSolveFactory<Scalar>::getUnderlyingLOWSF()
{
  return lowsf_.getNonconstObj();
}

template <class Scalar>
RCP<const LinearOpWithSolveFactoryBase<Scalar> >
ReuseLinearOpWithSolveFactory<Scalar>::getUnderlyingLOWSF() const
{
  return lowsf_.getConstObj();
}

template <class Scalar>
RCP<PreconditionerBase<Scalar> >
ReuseLinearOpWithSolveFactory<Scalar>::getUnderlyingPreconditioner()
{
  return prec_;
}

template <class Scalar>
RCP<const PreconditionerBase<Scalar> >
ReuseLinearOpWithSolveFactory<Scalar>::getUnderlyingPreconditioner() const
{
  return prec_;
}

// Overridden from Teuchos::Describable

template <class Scalar>
std::string ReuseLinearOpWithSolveFactory<Scalar>::description() const
{
  std::ostringstream oss;
  oss << this->Teuchos::Describable::description() << "{"
      << "lowsf=";
  if (!is_null(lowsf_.getConstObj()))
    oss << lowsf_.getConstObj()->description();
  else
    oss << "NULL";
  oss << std::endl
      << "prec=";
  if (!is_null(prec_))
    oss << prec_->description();
  else
    oss << "NULL";
  oss << "}";
  return oss.str();
}

// Overridden from ParameterListAcceptor

template <class Scalar>
void ReuseLinearOpWithSolveFactory<Scalar>::setParameterList(
    RCP<ParameterList> const &paramList)
{
  lowsf_.getNonconstObj()->setParameterList(paramList);
}

template <class Scalar>
RCP<ParameterList>
ReuseLinearOpWithSolveFactory<Scalar>::getNonconstParameterList()
{
  return lowsf_.getNonconstObj()->getNonconstParameterList();
}

template <class Scalar>
RCP<ParameterList> ReuseLinearOpWithSolveFactory<Scalar>::unsetParameterList()
{
  return lowsf_.getNonconstObj()->unsetParameterList();
}

template <class Scalar>
RCP<const ParameterList>
ReuseLinearOpWithSolveFactory<Scalar>::getParameterList() const
{
  return lowsf_.getConstObj()->getParameterList();
}

template <class Scalar>
RCP<const ParameterList>
ReuseLinearOpWithSolveFactory<Scalar>::getValidParameters() const
{
  return lowsf_.getConstObj()->getValidParameters();
}

// Overridden from LinearOpWithSolveFactoyBase

template <class Scalar>
bool ReuseLinearOpWithSolveFactory<Scalar>::acceptsPreconditionerFactory() const
{
  return false;
}

template <class Scalar>
void ReuseLinearOpWithSolveFactory<Scalar>::setPreconditionerFactory(
    const RCP<PreconditionerFactoryBase<Scalar> > & /* precFactory */,
    const std::string & /* precFactoryName */
)
{
}

template <class Scalar>
RCP<PreconditionerFactoryBase<Scalar> >
ReuseLinearOpWithSolveFactory<Scalar>::getPreconditionerFactory() const
{
  return Thyra::reusePreconditionerFactory<Scalar>(prec_);
}

template <class Scalar>
void ReuseLinearOpWithSolveFactory<Scalar>::unsetPreconditionerFactory(
    RCP<PreconditionerFactoryBase<Scalar> > * /* precFactory */,
    std::string * /* precFactoryName */
)
{
}

template <class Scalar>
bool ReuseLinearOpWithSolveFactory<Scalar>::isCompatible(
    const LinearOpSourceBase<Scalar> &fwdOpSrc) const
{
  return lowsf_.getConstObj()->isCompatible(fwdOpSrc);
}

template <class Scalar>
RCP<LinearOpWithSolveBase<Scalar> >
ReuseLinearOpWithSolveFactory<Scalar>::createOp() const
{
  return lowsf_.getConstObj()->createOp();
}

template <class Scalar>
void ReuseLinearOpWithSolveFactory<Scalar>::initializeOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse) const
{
  lowsf_.getConstObj()->initializeOp(fwdOpSrc, Op, supportSolveUse);
}

template <class Scalar>
void ReuseLinearOpWithSolveFactory<Scalar>::initializeAndReuseOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op) const
{
  lowsf_.getConstObj()->initializeAndReuseOp(fwdOpSrc, Op);
}

template <class Scalar>
void ReuseLinearOpWithSolveFactory<Scalar>::uninitializeOp(
    LinearOpWithSolveBase<Scalar> *Op,
    RCP<const LinearOpSourceBase<Scalar> > *fwdOpSrc,
    RCP<const PreconditionerBase<Scalar> > *prec,
    RCP<const LinearOpSourceBase<Scalar> > *approxFwdOpSrc,
    ESupportSolveUse *supportSolveUse) const
{
  lowsf_.getConstObj()->uninitializeOp(Op, fwdOpSrc, prec, approxFwdOpSrc,
                                       supportSolveUse);
}

template <class Scalar>
bool ReuseLinearOpWithSolveFactory<Scalar>::supportsPreconditionerInputType(
    const EPreconditionerInputType precOpType) const
{
  return lowsf_.getConstObj()->supportsPreconditionerInputType(precOpType);
}

template <class Scalar>
void ReuseLinearOpWithSolveFactory<Scalar>::initializePreconditionedOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const RCP<const PreconditionerBase<Scalar> > &prec,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse) const
{
  lowsf_.getConstObj()->initializePreconditionedOp(fwdOpSrc, prec, Op,
                                                   supportSolveUse);
}

template <class Scalar>
void ReuseLinearOpWithSolveFactory<Scalar>::initializeApproxPreconditionedOp(
    const RCP<const LinearOpSourceBase<Scalar> > &fwdOpSrc,
    const RCP<const LinearOpSourceBase<Scalar> > &approxFwdOpSrc,
    LinearOpWithSolveBase<Scalar> *Op,
    const ESupportSolveUse supportSolveUse) const
{
  lowsf_.getConstObj()->initializeApproxPreconditionedOp(
      fwdOpSrc, approxFwdOpSrc, Op, supportSolveUse);
}

// protected

template <class Scalar>
void ReuseLinearOpWithSolveFactory<Scalar>::informUpdatedVerbosityState() const
{
  lowsf_.getConstObj()->setVerbLevel(this->getVerbLevel());
  lowsf_.getConstObj()->setOStream(this->getOStream());
}

}  // namespace Thyra

#endif
