//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Thyra_MultiVectorPreconditioner_hpp
#define Thyra_MultiVectorPreconditioner_hpp

#include "Thyra_PreconditionerBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Thyra_MultiVectorLinearOp.hpp"
#include "Thyra_DefaultMultiVectorProductVectorSpace.hpp"

namespace Thyra {

/** \brief Concrete <tt>PreconditionerBase</tt> subclass that
 * wraps a preconditioner operator in MultiVectorLinearOp.
 */
template <class Scalar>
class MultiVectorPreconditioner : virtual public PreconditionerBase<Scalar> {
 public:
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  MultiVectorPreconditioner() {}

  void nonconstInitialize(
      const RCP<PreconditionerBase<Scalar> > &prec,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecRange,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecDomain)
  {
    validateInitialize(prec, multiVecRange, multiVecDomain);
    prec_           = prec;
    multiVecRange_  = multiVecRange;
    multiVecDomain_ = multiVecDomain;
  }

  void initialize(const RCP<const PreconditionerBase<Scalar> > &prec,
                  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
                      &multiVecRange,
                  const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
                      &multiVecDomain)
  {
    validateInitialize(prec, multiVecRange, multiVecDomain);
    prec_           = prec;
    multiVecRange_  = multiVecRange;
    multiVecDomain_ = multiVecDomain;
  }

  RCP<PreconditionerBase<Scalar> > getNonconstPreconditioner()
  {
    return prec_.getNonconstObj();
  }

  RCP<const PreconditionerBase<Scalar> > getPreconditioner() const
  {
    return prec_.getConstObj();
  }

  void uninitialize()
  {
    prec_.uninitialize();
    multiVecRange_  = Teuchos::null;
    multiVecDomain_ = Teuchos::null;
  }

  //@}

  /** @name Overridden from PreconditionerBase */
  //@{

  bool isLeftPrecOpConst() const
  {
    return prec_.getConstObj()->isLeftPrecOpConst();
  }

  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstLeftPrecOp()
  {
    return nonconstMultiVectorLinearOp(
        prec_.getNonconstObj()->getNonconstLeftPrecOp(), multiVecRange_,
        multiVecDomain_);
  }

  Teuchos::RCP<const LinearOpBase<Scalar> > getLeftPrecOp() const
  {
    return multiVectorLinearOp(prec_.getConstObj()->getLeftPrecOp(),
                               multiVecRange_, multiVecDomain_);
  }

  bool isRightPrecOpConst() const
  {
    return prec_.getConstObj()->isRightPrecOpConst();
  }

  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstRightPrecOp()
  {
    return nonconstMultiVectorLinearOp(
        prec_.getNonconstObj()->getNonconstRightPrecOp(), multiVecRange_,
        multiVecDomain_);
  }

  Teuchos::RCP<const LinearOpBase<Scalar> > getRightPrecOp() const
  {
    return multiVectorLinearOp(prec_.getConstObj()->getRightPrecOp(),
                               multiVecRange_, multiVecDomain_);
  }

  bool isUnspecifiedPrecOpConst() const
  {
    return prec_.getConstObj()->isUnspecifiedPrecOpConst();
  }

  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstUnspecifiedPrecOp()
  {
    return nonconstMultiVectorLinearOp(
        prec_.getNonconstObj()->getNonconstUnspecifiedPrecOp(), multiVecRange_,
        multiVecDomain_);
  }

  Teuchos::RCP<const LinearOpBase<Scalar> > getUnspecifiedPrecOp() const
  {
    return multiVectorLinearOp(prec_.getNonconstObj()->getUnspecifiedPrecOp(),
                               multiVecRange_, multiVecDomain_);
  }

  //@}

 private:
  // //////////////////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<PreconditionerBase<Scalar> >
      CNPB;

  // //////////////////////////////
  // Private data members

  CNPB prec_;
  RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > multiVecRange_;
  RCP<const DefaultMultiVectorProductVectorSpace<Scalar> > multiVecDomain_;

  // //////////////////////////////
  // Private member functions

  static void validateInitialize(
      const RCP<const PreconditionerBase<Scalar> > &prec,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecRange,
      const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
          &multiVecDomain)
  {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(is_null(prec));
    TEUCHOS_TEST_FOR_EXCEPT(is_null(multiVecRange));
    TEUCHOS_TEST_FOR_EXCEPT(is_null(multiVecDomain));
    TEUCHOS_TEST_FOR_EXCEPT(multiVecRange->numBlocks() !=
                            multiVecDomain->numBlocks());
#else
    (void)prec;
    (void)multiVecRange;
    (void)multiVecDomain;
#endif
  }
};

/** \brief Nonmember constructor function.
 *
 * \relates MultiVectorPreconditioner
 */
template <class Scalar>
RCP<MultiVectorPreconditioner<Scalar> > multiVectorPreconditioner()
{
  return Teuchos::rcp(new MultiVectorPreconditioner<Scalar>());
}

/** \brief Nonmember constructor function.
 *
 * \relates MultiVectorPreconditioner
 */
template <class Scalar>
RCP<MultiVectorPreconditioner<Scalar> > nonconstMultiVectorPreconditioner(
    const RCP<PreconditionerBase<Scalar> > &prec,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecRange,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecDomain)
{
  RCP<MultiVectorPreconditioner<Scalar> > mvprec =
      Teuchos::rcp(new MultiVectorPreconditioner<Scalar>());
  mvprec->nonconstInitialize(prec, multiVecRange, multiVecDomain);
  return mvprec;
}

/** \brief Nonmember constructor function.
 *
 * \relates MultiVectorPreconditioner
 */
template <class Scalar>
RCP<MultiVectorPreconditioner<Scalar> > multiVectorPreconditioner(
    const RCP<const PreconditionerBase<Scalar> > &prec,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecRange,
    const RCP<const DefaultMultiVectorProductVectorSpace<Scalar> >
        &multiVecDomain)
{
  RCP<MultiVectorPreconditioner<Scalar> > mvprec =
      Teuchos::rcp(new MultiVectorPreconditioner<Scalar>());
  mvprec->initialize(prec, multiVecRange, multiVecDomain);
  return mvprec;
}

}  // end namespace Thyra

#endif
