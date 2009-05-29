#ifndef __tParallelInverse_hpp__
#define __tParallelInverse_hpp__

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

// PB includes
#include "PB_Utilities.hpp"

#include <string>

#include "Test_Utils.hpp"

namespace PB {
namespace Test {

class tParallelInverse : public UnitTest {
public:
   virtual ~tParallelInverse() {}

   virtual void initializeTest();
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun);
   virtual bool isParallel() const { return true; }

   bool test_inverse(int verbosity,std::ostream & os);
   bool test_stridedInverse(int verbosity,std::ostream & os);

protected:
   void loadMatrix();
   void loadStridedMatrix();

   double tolerance_;

   // matrix to invert
   PB::LinearOp F_;
};

} // end namespace Tests
} // end namespace PB

#endif
