#ifndef __Teko_PreconditionerLinearOpDecl_hpp__
#define __Teko_PreconditionerLinearOpDecl_hpp__

#include "Thyra_LinearOpBase.hpp"
#include "Thyra_PreconditionerBase.hpp"
#include "Thyra_VectorSpaceBase.hpp"

#include "Teuchos_ConstNonconstObjectContainer.hpp"

namespace Teko {

/** \brief Class that wraps a <code>PreconditionerBase</code> object it makes it behave
  *        like a linear operator.
  *        
  * Class that wraps a <code>PreconditionerBase</code> object it makes it behave
  * like a linear operator.
  */
template <typename ScalarT>
class PreconditionerLinearOp : public Thyra::LinearOpBase<ScalarT> {
public:
   PreconditionerLinearOp();
   PreconditionerLinearOp(const Teuchos::RCP<Thyra::PreconditionerBase<ScalarT> > & prec);
   PreconditionerLinearOp(const Teuchos::RCP<const Thyra::PreconditionerBase<ScalarT> > & prec);

   //! build a linear operator using this preconditioner, this initialization permits changes
   void initialize(const Teuchos::RCP<Thyra::PreconditionerBase<ScalarT> > & prec);

   //! build a linear operator using this preconditioner, this initialization refuses changes
   void initialize(const Teuchos::RCP<const Thyra::PreconditionerBase<ScalarT> > & prec);

   //! Disassociate this object with the currently owned preconditioner
   void uninitialize();

   /** @brief Range space of this operator */
   virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > range() const;

   /** @brief Domain space of this operator */
   virtual Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > domain() const;

   virtual bool opSupportedImpl(const Thyra::EOpTransp M_trans) const;

   //! @brief Apply operation
   virtual void applyImpl(
     const Thyra::EOpTransp M_trans,
     const Thyra::MultiVectorBase<ScalarT> & x,
     const Teuchos::Ptr<Thyra::MultiVectorBase<ScalarT> > & y,
     const ScalarT alpha,
     const ScalarT beta
     ) const;

   //! Get a nonconstant <code>PreconditionerBase</code> object
   virtual Teuchos::RCP<Thyra::PreconditionerBase<ScalarT> > getNonconstPreconditioner();

   //! Get a constant <code>PreconditionerBase</code> object
   virtual Teuchos::RCP<const Thyra::PreconditionerBase<ScalarT> > getPreconditioner() const;

protected:
   //! get operator associated with the preconditioner
   Teuchos::ConstNonconstObjectContainer<Thyra::LinearOpBase<ScalarT> > getOperator() const;

   //! get operator associated with the preconditioner
   Teuchos::ConstNonconstObjectContainer<Thyra::LinearOpBase<ScalarT> > getOperator();

   Teuchos::ConstNonconstObjectContainer<Thyra::PreconditionerBase<ScalarT> > preconditioner_;
};

} // end namespace Teko

#endif
