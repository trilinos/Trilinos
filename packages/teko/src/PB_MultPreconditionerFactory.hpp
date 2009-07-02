#ifndef __PB_MultiPreconditionerFactory_hpp__
#define __PB_MultiPreconditionerFactory_hpp__

#include "PB_BlockPreconditionerFactory.hpp"
#include "PB_Utilities.hpp"
#include "PB_BlockImplicitLinearOp.hpp"

namespace PB {

/** Preconditioning factories must supply a 'State' class which
  * is where data specific to the preconditioner construction 
  * is stored. The constructor will be invoked elsewhere.
  */
class MultPrecondState : public PB::BlockPreconditionerState {
public:
   MultPrecondState() {}

   Teuchos::RCP<BlockPreconditionerState> StateOne_;
   Teuchos::RCP<BlockPreconditionerState> StateTwo_;
};

/** Declaration of a PB::BlockImplicitLinearOp.
  * BlockImplicitLinearOp's are useful
  * when there is no simple or cheap matrix representation 
  * of something like f(r). This particular class 
  * corresponds to
  *         f(r) = (M2 + M1 - M2 * A * M1) r
  * 
  * which is an application of a multiplicative preconditioner.
  * Notice that when M1 = inv(A) or when M2 = inv(A), we get
  * f(r) = inv(A)*r. 
  *
  * While the above matrix represenation could be used
  * instead of writing an implicit function, it requires
  * an additional matrix-vector product than an efficient
  * implementation. It should be noted (see comments below)
  * that there is an efficient matrix represenation of this
  * f(r). Namely,
  *
  *             f(r) =  [M2 I] [I  -A] [ I]
  *                            [0   I] [ M1]
  *
  * so we didn't really need to create an implicit operator.
  */
class MultPrecsLinearOp : public PB::BlockImplicitLinearOp {
public:

   //! Constructor
   MultPrecsLinearOp(const PB::LinearOp &A, const PB::LinearOp &M1, 
            const PB::LinearOp &M2): A_(A), M1_(M1), M2_(M2) { }

   virtual PB::VectorSpace  range() const { return M1_->range(); }
   virtual PB::VectorSpace domain() const { return M1_->domain();}
   virtual void implicitApply(const PB::BlockedMultiVector & r, PB::BlockedMultiVector & y,
             const double alpha = 1.0, const double beta = 0.0) const;

protected:
   PB::LinearOp A_, M1_, M2_;

private:
   // hide me!
   MultPrecsLinearOp();
   MultPrecsLinearOp(const MultPrecsLinearOp &);
};

/** Declaration of preconditioner factory that creates
  * a preconditioner which is the multiplicative combination
  * of two other preconditioners.
  */
class MultPreconditionerFactory
   : public PB::BlockPreconditionerFactory {
public:
   //! Constructor
   MultPreconditionerFactory(const Teuchos::RCP<const PB::BlockPreconditionerFactory> & FirstFactory,
                             const Teuchos::RCP<const PB::BlockPreconditionerFactory> & SecondFactory);

   MultPreconditionerFactory();

   //! Function inherited from PB::BlockPreconditionerFactory
   PB::LinearOp buildPreconditionerOperator(PB::BlockedLinearOp & blo,
                                            PB::BlockPreconditionerState & state) const;
    
   //! Build the MultPrecondState object
   virtual Teuchos::RCP<PB::BlockPreconditionerState> buildPreconditionerState() const;

protected:
   // class members
   Teuchos::RCP<const PB::BlockPreconditionerFactory> FirstFactory_;
   Teuchos::RCP<const PB::BlockPreconditionerFactory> SecondFactory_;
   
   //! Initialize from a parameter list
   virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);
};

} // end namespace PB

#endif
