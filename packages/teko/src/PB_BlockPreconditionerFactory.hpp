#ifndef __PB_BlockPreconditionerFactory_hpp__
#define __PB_BlockPreconditionerFactory_hpp__

#include "Teuchos_ParameterListAcceptor.hpp"

// Thyra includes
#include "Thyra_SolveSupportTypes.hpp"
#include "Thyra_LinearOpSourceBase.hpp"
#include "Thyra_PreconditionerFactoryBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

// PB includes
#include "PB_Utilities.hpp"

namespace PB {

using Teuchos::RCP;
using Teuchos::ParameterList;

inline BlockedLinearOp createNewBlockedOp()
{
   return rcp(new Thyra::DefaultBlockedLinearOp<double>());
}

class BlockPreconditionerFactory 
   : public virtual Thyra::PreconditionerFactoryBase<double> {
   public:
      // for BlockedPreconditionerFactory
      ///////////////////////////////////////////////////////////////////////

      virtual LinearOp buildPreconditionerOperator(BlockedLinearOp & blo) const = 0;

 
      //! @name Ignore me methods.
      //@{
 
      virtual LinearOp buildPreconditionerOperator(LinearOp & blo) const;

      //! is this operator compatiable with the preconditioner factory?
      bool isCompatible(const Thyra::LinearOpSourceBase<double> &fwdOpSrc) const;

      //! create an instance of the preconditioner
      RCP<Thyra::PreconditionerBase<double> > createPrec() const;

      //! initialize a newly created preconditioner object
      void initializePrec(const RCP<const Thyra::LinearOpSourceBase<double> > & fwdOpSrc,
                          Thyra::PreconditionerBase<double> * precOp,
                          const Thyra::ESupportSolveUse supportSolveUse) const;

      //! wipe clean a already initialized preconditioner object
      void uninitializePrec(Thyra::PreconditionerBase<double> * prec, 
                            RCP<const Thyra::LinearOpSourceBase<double> > * fwdOpSrc,
                            Thyra::ESupportSolveUse *supportSolveUse) const;

      // for ParameterListAcceptor
      ///////////////////////////////////////////////////////////////////////

      //! Set parameters from a parameter list and return with default values.
      void setParameterList(const RCP<ParameterList> & paramList); 

      //! Get the parameter list that was set using setParameterList().
      RCP< ParameterList > getNonconstParameterList();

      //! Unset the parameter list that was set using setParameterList(). 
      RCP< ParameterList > unsetParameterList();
      //@}

   protected:
      //! for ParameterListAcceptor
      mutable RCP<ParameterList>  validPL_;
      RCP<ParameterList>          paramList_;
};

} // end namespace PB

#endif
