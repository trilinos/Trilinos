#ifndef __Teko_PreconditionerLinearOp_hpp__
#define __Teko_PreconditionerLinearOp_hpp__

#include "Teko_PreconditionerLinearOpDecl.hpp"

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_PreconditionerBase.hpp"

namespace Teko {

template <typename ScalarT>
PreconditionerLinearOp<ScalarT>::PreconditionerLinearOp()
{ }

template <typename ScalarT>
PreconditionerLinearOp<ScalarT>::PreconditionerLinearOp(const Teuchos::RCP<Thyra::PreconditionerBase<ScalarT> > & prec)
{
   preconditioner_.initialize(prec);
}

template <typename ScalarT>
PreconditionerLinearOp<ScalarT>::PreconditionerLinearOp(const Teuchos::RCP<const Thyra::PreconditionerBase<ScalarT> > & prec)
{
   preconditioner_.initialize(prec);
}

//! build a linear operator using this preconditioner, this initialization permits changes
template <typename ScalarT>
void PreconditionerLinearOp<ScalarT>::initialize(const Teuchos::RCP<Thyra::PreconditionerBase<ScalarT> > & prec)
{
   uninitialize();
   preconditioner_.initialize(prec);
}

//! build a linear operator using this preconditioner, this initialization refuses changes
template <typename ScalarT>
void PreconditionerLinearOp<ScalarT>::initialize(const Teuchos::RCP<const Thyra::PreconditionerBase<ScalarT> > & prec)
{
   uninitialize();
   preconditioner_.initialize(prec);
}

//! Disassociate this object with the currently owned preconditioner
template <typename ScalarT>
void PreconditionerLinearOp<ScalarT>::uninitialize()
{
   preconditioner_.uninitialize();
}

/** @brief Range space of this operator */
template <typename ScalarT>
Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > PreconditionerLinearOp<ScalarT>::range() const
{
   return getOperator()->range();
}

/** @brief Domain space of this operator */
template <typename ScalarT>
Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > PreconditionerLinearOp<ScalarT>::domain() const
{
   return getOperator()->domain();
}

template <typename ScalarT>
bool PreconditionerLinearOp<ScalarT>::opSupportedImpl(
  const Thyra::EOpTransp M_trans) const
{
  return getOperator()->opSupported(M_trans);
}

template <typename ScalarT>
void PreconditionerLinearOp<ScalarT>::applyImpl(
  const Thyra::EOpTransp M_trans,
  const Thyra::MultiVectorBase<ScalarT> & x,
  const Teuchos::Ptr<Thyra::MultiVectorBase<ScalarT> > & y,
  const ScalarT alpha,
  const ScalarT beta
  ) const
{
   getOperator()->apply(M_trans, x, y, alpha, beta);
}


//! Get a nonconstant <code>PreconditionerBase</code> object
template <typename ScalarT>
Teuchos::RCP<Thyra::PreconditionerBase<ScalarT> > PreconditionerLinearOp<ScalarT>::getNonconstPreconditioner()
{
   return preconditioner_.getNonconstObj();
}

//! Get a constant <code>PreconditionerBase</code> object
template <typename ScalarT>
Teuchos::RCP<const Thyra::PreconditionerBase<ScalarT> > PreconditionerLinearOp<ScalarT>::getPreconditioner() const
{
   return preconditioner_.getConstObj();
}

//! get operator associated with the preconditioner
template <typename ScalarT>
Teuchos::ConstNonconstObjectContainer<Thyra::LinearOpBase<ScalarT> > PreconditionerLinearOp<ScalarT>::getOperator() const
{
   Teuchos::ConstNonconstObjectContainer<Thyra::LinearOpBase<ScalarT> > oper;
   oper.initialize(preconditioner_.getConstObj()->getUnspecifiedPrecOp());

   return oper;
}

//! get operator associated with the preconditioner
template <typename ScalarT>
Teuchos::ConstNonconstObjectContainer<Thyra::LinearOpBase<ScalarT> > PreconditionerLinearOp<ScalarT>::getOperator()
{
   Teuchos::ConstNonconstObjectContainer<Thyra::LinearOpBase<ScalarT> > oper;
   oper.initialize(preconditioner_.getNonconstObj()->getNonconstUnspecifiedPrecOp());

   return oper;
}

} // end namespace Teko

#endif
