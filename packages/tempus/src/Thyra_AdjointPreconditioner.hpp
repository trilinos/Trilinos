//@HEADER
// *****************************************************************************
//          Tempus: Time Integration and Sensitivity Analysis Package
//
// Copyright 2017 NTESS and the Tempus contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
//@HEADER

#ifndef Thyra_AdjointPreconditioner_hpp
#define Thyra_AdjointPreconditioner_hpp

#include "Thyra_PreconditionerBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"

namespace Thyra {

/** \brief Concrete <tt>PreconditionerBase</tt> subclass that
 * wraps a preconditioner operator in MultiVectorLinearOp.
 */
template <class Scalar>
class AdjointPreconditioner : virtual public PreconditionerBase<Scalar> {
 public:
  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  AdjointPreconditioner() {}

  void nonconstInitialize(const RCP<PreconditionerBase<Scalar> > &prec)
  {
    validateInitialize(prec);
    prec_ = prec;
  }

  void initialize(const RCP<const PreconditionerBase<Scalar> > &prec)
  {
    validateInitialize(prec);
    prec_ = prec;
  }

  RCP<PreconditionerBase<Scalar> > getNonconstPreconditioner()
  {
    return prec_.getNonconstObj();
  }

  RCP<const PreconditionerBase<Scalar> > getPreconditioner() const
  {
    return prec_.getConstObj();
  }

  void uninitialize() { prec_.uninitialize(); }

  //@}

  /** @name Overridden from PreconditionerBase */
  //@{

  bool isLeftPrecOpConst() const
  {
    return prec_.getConstObj()->isLeftPrecOpConst();
  }

  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstLeftPrecOp()
  {
    return nonconstAdjoint(prec_.getNonconstObj()->getNonconstLeftPrecOp());
  }

  Teuchos::RCP<const LinearOpBase<Scalar> > getLeftPrecOp() const
  {
    return adjoint(prec_.getConstObj()->getLeftPrecOp());
  }

  bool isRightPrecOpConst() const
  {
    return prec_.getConstObj()->isRightPrecOpConst();
  }

  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstRightPrecOp()
  {
    return nonconstAdjoint(prec_.getNonconstObj()->getNonconstRightPrecOp());
  }

  Teuchos::RCP<const LinearOpBase<Scalar> > getRightPrecOp() const
  {
    return adjoint(prec_.getConstObj()->getRightPrecOp());
  }

  bool isUnspecifiedPrecOpConst() const
  {
    return prec_.getConstObj()->isUnspecifiedPrecOpConst();
  }

  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstUnspecifiedPrecOp()
  {
    return nonconstAdjoint(
        prec_.getNonconstObj()->getNonconstUnspecifiedPrecOp());
  }

  Teuchos::RCP<const LinearOpBase<Scalar> > getUnspecifiedPrecOp() const
  {
    return adjoint(prec_.getNonconstObj()->getUnspecifiedPrecOp());
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

  // //////////////////////////////
  // Private member functions

  static void validateInitialize(
      const RCP<const PreconditionerBase<Scalar> > &prec)
  {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(is_null(prec));
#else
    (void)prec;
#endif
  }
};

/** \brief Nonmember constructor function.
 *
 * \relates AdjointPreconditioner
 */
template <class Scalar>
RCP<AdjointPreconditioner<Scalar> > adjointPreconditioner()
{
  return Teuchos::rcp(new AdjointPreconditioner<Scalar>());
}

/** \brief Nonmember constructor function.
 *
 * \relates AdjointPreconditioner
 */
template <class Scalar>
RCP<AdjointPreconditioner<Scalar> > nonconstAdjointPreconditioner(
    const RCP<PreconditionerBase<Scalar> > &prec)
{
  RCP<AdjointPreconditioner<Scalar> > aprec =
      Teuchos::rcp(new AdjointPreconditioner<Scalar>());
  aprec->nonconstInitialize(prec);
  return aprec;
}

/** \brief Nonmember constructor function.
 *
 * \relates AdjointPreconditioner
 */
template <class Scalar>
RCP<AdjointPreconditioner<Scalar> > adjointPreconditioner(
    const RCP<const PreconditionerBase<Scalar> > &prec)
{
  RCP<AdjointPreconditioner<Scalar> > aprec =
      Teuchos::rcp(new AdjointPreconditioner<Scalar>());
  aprec->initialize(prec);
  return aprec;
}

}  // end namespace Thyra

#endif
