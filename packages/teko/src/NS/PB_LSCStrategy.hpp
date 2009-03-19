#ifndef __PB_LSCPreconditionerStrategy_hpp__
#define __PB_LSCPreconditionerStrategy_hpp__

#include "Teuchos_RCP.hpp"

#include "Thyra_LinearOpBase.hpp"

namespace PB {
namespace NS {

// simple strategy for driving LSCPreconditionerFactory
class LSCStrategy {
public:
   virtual Teuchos::RCP<const Thyra::LinearOpBase<double> > 
   getInvBQBt(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const = 0;

   virtual Teuchos::RCP<const Thyra::LinearOpBase<double> > 
   getInvF(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const = 0;

   virtual Teuchos::RCP<const Thyra::LinearOpBase<double> > 
   getInvD(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const = 0;

   virtual Teuchos::RCP<const Thyra::LinearOpBase<double> > 
   getInvMass(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const = 0;
};

// constant, not very flexible strategy for driving LSCPreconditioenrFactory
class StaticLSCStrategy : public LSCStrategy {
public:

   // Staiblized constructor
   StaticLSCStrategy(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invF,
                     const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invBQBtmC,
                     const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invD,
                     const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invMass);
 
   // Stable constructor
   StaticLSCStrategy(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invF,
                     const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invBQBtmC,
                     const Teuchos::RCP<const Thyra::LinearOpBase<double> > & invMass);

   virtual Teuchos::RCP<const Thyra::LinearOpBase<double> > 
   getInvF(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const
   { return invF_; }

   virtual Teuchos::RCP<const Thyra::LinearOpBase<double> > 
   getInvBQBt(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const
   { return invBQBtmC_; }

   virtual Teuchos::RCP<const Thyra::LinearOpBase<double> > 
   getInvD(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const
   { return invD_; }

   virtual Teuchos::RCP<const Thyra::LinearOpBase<double> > 
   getInvMass(const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A) const
   { return invMass_; }

protected:
   // protected memebers
   Teuchos::RCP<const Thyra::LinearOpBase<double> > invF_;
   Teuchos::RCP<const Thyra::LinearOpBase<double> > invBQBtmC_;
   Teuchos::RCP<const Thyra::LinearOpBase<double> > invD_;
   Teuchos::RCP<const Thyra::LinearOpBase<double> > invMass_;
};

} // end namespace NS
} // end namespace PB

#endif
