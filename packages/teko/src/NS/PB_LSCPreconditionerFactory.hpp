#ifndef __PB_LSCPreconditionerFactory_hpp__
#define __PB_LSCPreconditionerFactory_hpp__

#include "PB_BlockPreconditionerFactory.hpp"
#include "PB_LSCStrategy.hpp"

namespace PB {
namespace NS { // Navier-Stokes specialization

/** \brief Preconditioner state for the LSC factory. 
  *
  * Preconditioner state for the LSC factory. This is based
  * on the notation and concepts found in
  *
  * Elman, Howle, Shadid, Silvester, and Tuminaro, "Least Squares Preconditioners
  * for Stabilized Discretizations of the Navier-Stokes Euqations," SISC-2007.
  */
class LSCPrecondState : public BlockPreconditionerState {
public:
   LSCPrecondState() : initialized_(false) {}

   //! Is this state constructed?
   bool isInitialized() { return initialized_; }

   //! Is this state constructed?
   void setInitialized(bool i) { initialized_ = i; }

   //! Inverse mass operator (\f$Q_u^{-1}\f$)
   LinearOp invMass_;

   /** \f$B Q_u^{-1} B^T\f$
     * \f$D = diag(B \; diag(F)^{-1} B^T + C)\f$.
     */
   LinearOp BQBt_;

   /** \f$B Q_u^{-1} B^T-\gamma C\f$
     */
   LinearOp BQBtmC_;
   InverseLinearOp invBQBtmC_;

   InverseLinearOp invF_;

   //! \f$\alpha D^{-1}\f$ where
   LinearOp aiD_;

   //! \f$\gamma = \rho(Q_u^{-1} F / 3)\f$
   double gamma_;

   /** \f$\alpha = 1/\rho(B \; diag(F)^{-1} B^T D^{-1})\f$ where
     * \f$D = diag(B \; diag(F)^{-1} B^T + C)\f$.
     */
   double alpha_;

protected:
   bool initialized_;
};

class LSCPreconditionerFactory 
   : public BlockPreconditionerFactory {
public:
   //! \name Constructors
   //@{

   //! Staiblized constructor
   LSCPreconditionerFactory(const LinearOp & invF,const LinearOp & invBQBtmC,
                            const LinearOp & invD,const LinearOp & invMass);
 
   //! Stable constructor
   LSCPreconditionerFactory(const LinearOp & invF,
                            const LinearOp & invBQBtmC,
                            const LinearOp & invMass);

   //! fully generic constructor
   LSCPreconditionerFactory(const Teuchos::RCP<const LSCStrategy> & strategy);

   //! Default constructor
   LSCPreconditionerFactory();
   //@}

   //! for PreconditionerFactoryBase
   virtual LinearOp buildPreconditionerOperator(BlockedLinearOp & blo,BlockPreconditionerState & state) const;

   //! Build the LSCPrecondState object
   virtual RCP<BlockPreconditionerState> buildPreconditionerState() const
   { return rcp(new LSCPrecondState()); }

protected:
   // Gimmie object
   Teuchos::RCP<const LSCStrategy> invOpsStrategy_;

   //! Initialize from a parameter list
   virtual void initializeFromParameterList(const Teuchos::ParameterList & pl);
};

} // end namespace NS
} // end namespace PB

#endif
