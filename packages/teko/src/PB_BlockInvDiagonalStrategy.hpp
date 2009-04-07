#ifndef __PB_BlockInvDiagonalStrategy_hpp__
#define __PB_BlockInvDiagonalStrategy_hpp__

#include <vector>

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

// PB includes
#include "PB_Utilities.hpp"
#include "PB_InverseFactory.hpp"
#include "PB_BlockPreconditionerFactory.hpp"

namespace PB {

/** this should be paired with a BlockJacobiPreconditionerFactory
  * or BlockGSPreconditionerFactory.  The idea is that this object
  * provides an approximate inverse operator for each of the diagonal
  * blocks.  Then, the [whatever]PreconditionerFactory can easily
  * construct an approximate inverse. The system under consideration
  * is
  * 
  *    A = [ D0  U01 U02 ...]
  *        [ L10  D1 U12 ...]
  *        [ L20 L21  D2 ...]
  *        [  .   .   .  ...]
  *
  * where inverses of D0,D1,D2...DN are needed.
  */
class BlockInvDiagonalStrategy {
public:
   //! returns an (approximate) inverse of the diagonal blocks of A
   virtual void getInvD(const BlockedLinearOp & A,BlockPreconditionerState & state,
                        std::vector<LinearOp> & invDiag) const = 0;

};

/** This is a simple strategy for a [whatever]PreconditionerFactory
  * it simply returns statically set RCP pointers to the passed in
  * inv(D0) and inv(D1) operators. Not this will _not_ permit 
  * efficient implementations when the preconditioner has to be rebuilt
  * or reused often.
  */
class StaticInvDiagStrategy : public BlockInvDiagonalStrategy {
public:
   StaticInvDiagStrategy(const LinearOp & invD0,
                          const LinearOp & invD1)
   { invDiag_.push_back(invD0); invDiag_.push_back(invD1); }

   StaticInvDiagStrategy(const std::vector<LinearOp> & invD)
      : invDiag_(invD)
   { }

   /** returns an (approximate) inverse of the diagonal blocks of A
     * where A is closely related to the original source for invD0 and invD1
     */
   virtual void getInvD(const BlockedLinearOp & A, BlockPreconditionerState & state,
                        std::vector<LinearOp> & invDiag) const
   { invDiag.clear(); invDiag = invDiag_; }

protected:
   // stored inverse operators
   std::vector<Teuchos::RCP<const Thyra::LinearOpBase<double> > > invDiag_;
};

class InvFactoryDiagStrategy : public BlockInvDiagonalStrategy {
public:
   InvFactoryDiagStrategy(const Teuchos::RCP<const InverseFactory> & factory);

   /** returns an (approximate) inverse of the diagonal blocks of A
     * where A is closely related to the original source for invD0 and invD1
     */
   virtual void getInvD(const BlockedLinearOp & A, BlockPreconditionerState & state,
                        std::vector<LinearOp> & invDiag) const;

protected:
   // stored inverse operators
   Teuchos::RCP<const InverseFactory> invDiagFact_;

private:
   InvFactoryDiagStrategy();
   InvFactoryDiagStrategy(const InvFactoryDiagStrategy &);
};

} // end namespace PB

#endif
