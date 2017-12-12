// @HEADER
// ****************************************************************************
//                Tempus: Copyright (2017) Sandia Corporation
//
// Distributed under BSD 3-clause license (See accompanying file Copyright.txt)
// ****************************************************************************
// @HEADER

#ifndef Thyra_AdjointPreconditioner_hpp
#define Thyra_AdjointPreconditioner_hpp

#include "Thyra_PreconditionerBase.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"

namespace Thyra {

/** \brief Concrete <tt>PreconditionerBase</tt> subclass that
 * wraps a preconditioner operator in MultiVectorLinearOp.
 */
template<class Scalar>
class AdjointPreconditioner : virtual public PreconditionerBase<Scalar>
{
public:

  /** @name Constructors/initializers/accessors */
  //@{

  /** \brief Construct to uninitialized. */
  AdjointPreconditioner() {}

  /** \brief . */
  void nonconstInitialize(
    const RCP<PreconditionerBase<Scalar> > &prec) {
    validateInitialize(prec);
    prec_ = prec;
  }

  /** \brief . */
  void initialize(
    const RCP<const PreconditionerBase<Scalar> > &prec) {
    validateInitialize(prec);
    prec_ = prec;
  }

  /** \brief . */
  RCP<PreconditionerBase<Scalar> >
  getNonconstPreconditioner() { return prec_.getNonconstObj(); }

  /** \brief . */
  RCP<const PreconditionerBase<Scalar> >
  getPreconditioner() const { return prec_.getConstObj(); }

  /** \brief . */
  void uninitialize() {
    prec_.uninitialize();
  }

  //@}

  /** @name Overridden from PreconditionerBase */
  //@{

  /** \brief . */
  bool isLeftPrecOpConst() const
  { return prec_.getConstObj()->isLeftPrecOpConst(); }

  /** \brief . */
  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstLeftPrecOp()
  { return nonconstAdjoint(prec_.getNonconstObj()->getNonconstLeftPrecOp()); }

  /** \brief . */
  Teuchos::RCP<const LinearOpBase<Scalar> > getLeftPrecOp() const
  { return adjoint(prec_.getConstObj()->getLeftPrecOp()); }

  /** \brief . */
  bool isRightPrecOpConst() const
  { return prec_.getConstObj()->isRightPrecOpConst(); }

  /** \brief . */
  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstRightPrecOp()
  { return nonconstAdjoint(prec_.getNonconstObj()->getNonconstRightPrecOp()); }

  /** \brief . */
  Teuchos::RCP<const LinearOpBase<Scalar> > getRightPrecOp() const
  { return adjoint(prec_.getConstObj()->getRightPrecOp()); }

  /** \brief . */
  bool isUnspecifiedPrecOpConst() const
  { return prec_.getConstObj()->isUnspecifiedPrecOpConst(); }

  /** \brief . */
  Teuchos::RCP<LinearOpBase<Scalar> > getNonconstUnspecifiedPrecOp()
  { return nonconstAdjoint(
      prec_.getNonconstObj()->getNonconstUnspecifiedPrecOp()); }

  /** \brief . */
  Teuchos::RCP<const LinearOpBase<Scalar> > getUnspecifiedPrecOp() const
  { return adjoint(prec_.getNonconstObj()->getUnspecifiedPrecOp()); }

  //@}

private:

  // //////////////////////////////
  // Private types

  typedef Teuchos::ConstNonconstObjectContainer<PreconditionerBase<Scalar> > CNPB;

  // //////////////////////////////
  // Private data members

  CNPB prec_;

  // //////////////////////////////
  // Private member functions

  static void validateInitialize(
    const RCP<const PreconditionerBase<Scalar> > &prec) {
#ifdef TEUCHOS_DEBUG
    TEUCHOS_TEST_FOR_EXCEPT(is_null(prec));
#endif
  }

};

/** \brief Nonmember constructor function.
 *
 * \relates AdjointPreconditioner
 */
template<class Scalar>
RCP<AdjointPreconditioner<Scalar> >
adjointPreconditioner()
{
  return Teuchos::rcp(new AdjointPreconditioner<Scalar>());
}

/** \brief Nonmember constructor function.
 *
 * \relates AdjointPreconditioner
 */
template<class Scalar>
RCP<AdjointPreconditioner<Scalar> >
nonconstAdjointPreconditioner(
  const RCP<PreconditionerBase<Scalar> > &prec
  )
{
  RCP<AdjointPreconditioner<Scalar> >
    aprec = Teuchos::rcp(new AdjointPreconditioner<Scalar>());
  aprec->nonconstInitialize(prec);
  return aprec;
}

/** \brief Nonmember constructor function.
 *
 * \relates AdjointPreconditioner
 */
template<class Scalar>
RCP<AdjointPreconditioner<Scalar> >
adjointPreconditioner(
  const RCP<const PreconditionerBase<Scalar> > &prec
  )
{
  RCP<AdjointPreconditioner<Scalar> >
    aprec = Teuchos::rcp(new AdjointPreconditioner<Scalar>());
  aprec->initialize(prec);
  return aprec;
}

}       // end namespace Thyra

#endif
