#ifndef __tBlockJacobiPreconditionerFactory_hpp__
#define __tBlockJacobiPreconditionerFactory_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

#include "Epetra_SerialComm.h"

#include <string>

#include "Test_Utils.hpp"

namespace PB {
namespace Test {

class tBlockJacobiPreconditionerFactory : public UnitTest {
public:
   virtual ~tBlockJacobiPreconditionerFactory() {}

   virtual void initializeTest();
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun);
   virtual bool isParallel() const { return true; }

   bool test_createPrec(int verbosity,std::ostream & os);
   bool test_initializePrec(int verbosity,std::ostream & os);
   bool test_uninitializePrec(int verbosity,std::ostream & os);
   bool test_isCompatible(int verbosity,std::ostream & os);

protected:
   double tolerance_;

   Teuchos::RCP<const Thyra::LinearOpBase<double> > F_,C_,B_,Bt_;
   Teuchos::RCP<const Thyra::LinearOpBase<double> > invF_,invC_;
};

} // end namespace Tests
} // end namespace PB

#endif
