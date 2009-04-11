#ifndef __PB_LSCPreconditionerStrategy_hpp__
#define __PB_LSCPreconditionerStrategy_hpp__

#include "Teuchos_RCP.hpp"

#include "Thyra_LinearOpBase.hpp"

#include "PB_Utilities.hpp"
#include "PB_InverseFactory.hpp"
#include "PB_BlockPreconditionerFactory.hpp"

namespace PB {
namespace NS {

class LSCPrecondState; // forward declration

// simple strategy for driving LSCPreconditionerFactory
class LSCStrategy {
public:
   virtual void buildState(BlockedLinearOp & A,BlockPreconditionerState & state) const = 0;

   virtual LinearOp getInvBQBt(const BlockedLinearOp & A,BlockPreconditionerState & state) const = 0;

   virtual LinearOp getInvF(const BlockedLinearOp & A,BlockPreconditionerState & state) const = 0;

   virtual LinearOp getInvD(const BlockedLinearOp & A,BlockPreconditionerState & state) const = 0;

   virtual LinearOp getInvMass(const BlockedLinearOp & A,BlockPreconditionerState & state) const = 0;
};

// constant, not very flexible strategy for driving LSCPreconditioenrFactory
class StaticLSCStrategy : public LSCStrategy {
public:

   // Staiblized constructor
   StaticLSCStrategy(const LinearOp & invF,
                     const LinearOp & invBQBtmC,
                     const LinearOp & invD,
                     const LinearOp & invMass);
 
   // Stable constructor
   StaticLSCStrategy(const LinearOp & invF,
                     const LinearOp & invBQBtmC,
                     const LinearOp & invMass);

   virtual void buildState(BlockedLinearOp & A,BlockPreconditionerState & state) const {}

   virtual LinearOp getInvF(const BlockedLinearOp & A,BlockPreconditionerState & state) const
   { return invF_; }

   virtual LinearOp getInvBQBt(const BlockedLinearOp & A,BlockPreconditionerState & state) const
   { return invBQBtmC_; }

   virtual LinearOp getInvD(const BlockedLinearOp & A,BlockPreconditionerState & state) const
   { return invD_; }

   virtual LinearOp getInvMass(const BlockedLinearOp & A,BlockPreconditionerState & state) const
   { return invMass_; }

protected:
   // protected memebers
   LinearOp invF_;
   LinearOp invBQBtmC_;
   LinearOp invD_;
   LinearOp invMass_;
};

/** \brief A strategy that takes a single inverse factory and
  *        uses that for all inverses. If no mass matrix is 
  *        passed in the diagonal of the 1,1 block is used.
  *
  * A strategy that takes a single inverse factory and uses that
  * for all inverses. Optionally the mass matrix can be passed
  * in, if it is the diagonal is extracted and that is used to
  * form the inverse approximation.
  */
class InvLSCStrategy : public LSCStrategy {
public:
   //! \name Constructors
   //@{
   InvLSCStrategy(const Teuchos::RCP<const InverseFactory> & factory);
   InvLSCStrategy(const Teuchos::RCP<const InverseFactory> & factory,LinearOp & mass);
   //@}

   //! Functions inherited from LSCStrategy
   //@{
   virtual void buildState(BlockedLinearOp & A,BlockPreconditionerState & state) const;

   virtual LinearOp getInvBQBt(const BlockedLinearOp & A,BlockPreconditionerState & state) const;

   virtual LinearOp getInvF(const BlockedLinearOp & A,BlockPreconditionerState & state) const;

   virtual LinearOp getInvD(const BlockedLinearOp & A,BlockPreconditionerState & state) const;

   virtual LinearOp getInvMass(const BlockedLinearOp & A,BlockPreconditionerState & state) const;
   //@}

   //! Initialize the state object using this blocked linear operator
   void initializeState(const BlockedLinearOp & A,LSCPrecondState * state) const;

   //! Initialize the state object using this blocked linear operator
   void reinitializeState(const BlockedLinearOp & A,LSCPrecondState * state) const;

protected:
   LinearOp massMatrix_;
   Teuchos::RCP<const InverseFactory> invFactory_;
};

} // end namespace NS
} // end namespace PB

#endif
