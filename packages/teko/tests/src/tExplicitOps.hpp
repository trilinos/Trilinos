#ifndef __tExplicitOps_hpp__
#define __tExplicitOps_hpp__

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

class tExplicitOps : public UnitTest {
public:
   virtual ~tExplicitOps() {}

   virtual void initializeTest();
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun);
   virtual bool isParallel() const { return false; }

   bool test_mult_diagScaleMatProd(int verbosity,std::ostream & os);
   bool test_mult_diagScaling(int verbosity,std::ostream & os);
   bool test_mult_modScaleMatProd(int verbosity,std::ostream & os);

protected:
   double tolerance_;

   // matrix to invert
   PB::ModifiableLinearOp F_;
   PB::ModifiableLinearOp G_;
   PB::LinearOp D_;
};

} // end namespace Tests
} // end namespace PB

#endif
