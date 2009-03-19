#ifndef __tInterlacedEpetra_hpp__
#define __tInterlacedEpetra_hpp__

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

class tInterlacedEpetra : public UnitTest {
public:
   virtual ~tInterlacedEpetra() {}

   virtual void initializeTest();
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun);
   virtual bool isParallel() const { return false; }

   bool test_buildSubMaps_num(int verbosity,std::ostream & os);
   bool test_buildSubMaps_vec(int verbosity,std::ostream & os);

protected:
};

} // end namespace Tests
} // end namespace PB

#endif
