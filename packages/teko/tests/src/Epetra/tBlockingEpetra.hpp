#ifndef __tBlockingEpetra_hpp__
#define __tBlockingEpetra_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

#include <string>

#include "Test_Utils.hpp"

namespace PB {
namespace Test {

class tBlockingEpetra : public UnitTest {
public:
   virtual ~tBlockingEpetra() {}

   virtual void initializeTest();
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun);
   virtual bool isParallel() const { return false; }

   bool test_buildMaps(int verbosity,std::ostream & os);
   bool test_one2many(int verbosity,std::ostream & os);
   bool test_many2one(int verbosity,std::ostream & os);

protected:
   double tolerance_;
};

} // end namespace Tests
} // end namespace PB

#endif
