#ifndef __tBlockLowerTriInverseOp_hpp__
#define __tBlockLowerTriInverseOp_hpp__

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

#include "Epetra_SerialComm.h"

#include <string>
#include <vector>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tBlockLowerTriInverseOp : public UnitTest {
public:
   virtual ~tBlockLowerTriInverseOp() {}

   virtual void initializeTest();
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun);
   virtual bool isParallel() const { return false; }

   bool test_apply(int verbosity,std::ostream & os);
   bool test_alphabeta(int verbosity,std::ostream & os);

protected:
   double tolerance_;
   Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > A_;
   Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > invA_;

   std::vector<Teuchos::RCP<const Thyra::LinearOpBase<double> > > invDiag_;
};

} // end namespace Tests
} // end namespace Teko

#endif
