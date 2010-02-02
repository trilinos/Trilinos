#ifndef __Teko_Preconditioner_hpp__
#define __Teko_Preconditioner_hpp__

// Thyra includes
#include "Thyra_DefaultPreconditioner.hpp"

// Teko includes
#include "Teko_Utilities.hpp"
#include "Teko_PreconditionerState.hpp"

namespace Teko {

using Thyra::LinearOpBase;
using Thyra::DefaultPreconditioner;
using Teuchos::RCP;

/** \brief An extension of the <code>Thyra::DefaultPreconditioner</code>
  *        class with some specializations useful for use within Teko.
  *
  * An extension of the <code>Thyra::DefaultPreconditioner</code>
  * class with some specializations useful for use within Teko. This includes
  * having facilities to store the source vector and the state object.
  */
class Preconditioner : public DefaultPreconditioner<double> {
public:
   //! \name Constructors based from Thyra::DefaultPreconditioner
   //@{
   Preconditioner()
      : DefaultPreconditioner<double>() {}
   Preconditioner(const RCP<LinearOpBase<double> > & leftPrecOp,
                  const RCP<LinearOpBase<double> > & rightPrecOp)
      : DefaultPreconditioner<double>(leftPrecOp,rightPrecOp) {}
   Preconditioner(const RCP<const LinearOpBase<double> > & leftPrecOp,
                  const RCP<const LinearOpBase<double> > & rightPrecOp)
      : DefaultPreconditioner<double>(leftPrecOp,rightPrecOp) {}
   Preconditioner(const RCP<LinearOpBase<double> > & unspecifiedPrecOp)
      : DefaultPreconditioner<double>(unspecifiedPrecOp) {}
   Preconditioner(const RCP<const LinearOpBase<double> > & unspecifiedPrecOp)
      : DefaultPreconditioner<double>(unspecifiedPrecOp) {}
   //@}

   /** Set the vector associated with this operator (think nonlinear system)
     *
     * \param[in] srcVec The source vector associated with this preconditioner.
     */
   virtual void setSourceVector(const RCP<Thyra::MultiVectorBase<double> > & srcVec)
   { if(srcVec!=Teuchos::null) state_->setSourceVector(srcVec);
     else                      state_->setSourceVector(Teuchos::null); }

   /** Set the state object associated with this preconditioner
     *
     * \param[in] state The state object to use
     */
   virtual void setStateObject(const RCP<PreconditionerState> & state)
   { state_ = state; }

   /** Get the state object associated with this preconditioner
     *
     * \returns The state object used by this preconditioner
     */
   virtual const RCP<PreconditionerState> getStateObject()
   { return state_; }

   /** Get the state object associated with this preconditioner
     *
     * \returns The state object used by this preconditioner
     */
   virtual const RCP<const PreconditionerState> getStateObject() const
   { return state_; }

protected:
   //! User defined state for this preconditioner
   RCP<PreconditionerState> state_;
};

} // end namespace Teko

#endif
