#ifndef __tEpetraLSCIntegrationTest_hpp__
#define __tEpetraLSCIntegrationTest_hpp__

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

class tEpetraLSCIntegrationTest : public UnitTest {
public:
   virtual ~tEpetraLSCIntegrationTest() {}

   virtual void initializeTest();
   virtual int runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun);
   virtual bool isParallel() const { return true; }

   bool test_withmassStable(int verbosity,std::ostream & os);
   bool test_nomassStable(int verbosity,std::ostream & os);

protected:
   void loadStableSystem();
   const Teuchos::RCP<const Epetra_Operator> mlSolve(const Teuchos::RCP<const Epetra_RowMatrix> mat,int cycles);

   double tolerance_;

   Teuchos::RCP<const Epetra_Map> velMap_;  // map of velocity space
   Teuchos::RCP<const Epetra_Map> prsMap_;  // map of pressure space
   Teuchos::RCP<const Epetra_Map> fullMap_; // map of pressure space

   // stable discretizations matrices
   Teuchos::RCP<const Epetra_CrsMatrix> sF_;
   Teuchos::RCP<const Epetra_CrsMatrix> sB_;
   Teuchos::RCP<const Epetra_CrsMatrix> sBt_;
   Teuchos::RCP<const Epetra_CrsMatrix> sQu_;
   Teuchos::RCP<Epetra_Operator> sA_;

   // stable rhs and IFISS solution
   Teuchos::RCP<Epetra_Vector> rhs_;
   Teuchos::RCP<const Epetra_Vector> sExact_;
};

} // end namespace Tests
} // end namespace PB

#endif
