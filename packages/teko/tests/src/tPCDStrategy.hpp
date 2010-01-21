#ifndef __tPCDStrategy_hpp__
#define __tPCDStrategy_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

// Teko includes
#include "Teko_Utilities.hpp"

#include <string>

#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

class tPCDStrategy : public UnitTest {
public:
   virtual ~tPCDStrategy() {}

   virtual void initializeTest();
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun);
   virtual bool isParallel() const { return false; }

   bool test_PCDStrategy(int verbosity,std::ostream & os);

protected:
   double tolerance_;
};

} // end namespace Tests
} // end namespace Teko

#endif
