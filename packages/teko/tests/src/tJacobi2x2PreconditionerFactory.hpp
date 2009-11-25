#ifndef __tJacobi2x2PreconditionerFactory_hpp__
#define __tJacobi2x2PreconditionerFactory_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"

#include "Epetra_SerialComm.h"

#include <string>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tJacobi2x2PreconditionerFactory : public UnitTest {
public:
   virtual ~tJacobi2x2PreconditionerFactory() {}

   virtual void initializeTest();
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun);
   virtual bool isParallel() const { return false; }

   bool test_createPrec(int verbosity,std::ostream & os);
   bool test_initializePrec(int verbosity,std::ostream & os);
   bool test_uninitializePrec(int verbosity,std::ostream & os);
   bool test_isCompatable(int verbosity,std::ostream & os);

   // non-member tests
   bool test_result(int verbosity,std::ostream & os);
   bool test_identity(int verbosity,std::ostream & os);
   bool test_diagonal(int verbosity,std::ostream & os);

protected:
   // some simple matrix subblocks
   Teuchos::RCP<const Thyra::LinearOpBase<double> > A_; 
   Teuchos::RCP<const Thyra::LinearOpBase<double> > F_; 
   Teuchos::RCP<const Thyra::LinearOpBase<double> > D_; 
   Teuchos::RCP<const Thyra::LinearOpBase<double> > G_; 
   Teuchos::RCP<const Thyra::LinearOpBase<double> > C_; 
   Teuchos::RCP<const Thyra::LinearOpBase<double> > invF_; 
   Teuchos::RCP<const Thyra::LinearOpBase<double> > invC_; 
   Teuchos::RCP<Epetra_SerialComm> comm;

   double tolerance_;
};

} // end namespace Tests
} // end namespace Teko

#endif
