#ifndef __tLSCHIntegrationTest_hpp__
#define __tLSCHIntegrationTest_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include <string>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tLSCHIntegrationTest : public UnitTest {
public:
   virtual ~tLSCHIntegrationTest() {}

   virtual void initializeTest();
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun);
   virtual bool isParallel() const { return false; }

   bool test_hScaling(int verbosity,std::ostream & os);

protected:
   double tolerance_;
};

} // end namespace Tests
} // end namespace Teko

#endif
