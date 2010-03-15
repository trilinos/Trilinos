#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <Teuchos_RCP.hpp>

#include <string>
#include <iostream>

#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
#else
   #include "Epetra_SerialComm.h"
#endif

#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_InverseFactory.hpp"
#include "Teko_LU2x2InverseOp.hpp"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

// Test-rig

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::epetraLinearOp;

const RCP<const Thyra::LinearOpBase<double> > build2x2(const Epetra_Comm & comm,double a,double b,double c,double d)
{
   RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,comm));

   int indicies[2];
   double row0[2];
   double row1[2];

   indicies[0] = 0;
   indicies[1] = 1;

   // build a CrsMatrix
   RCP<Epetra_CrsMatrix> blk  = rcp(new Epetra_CrsMatrix(Copy,*map,2));
   row0[0] = a; row0[1] = b;  // do a transpose here!
   row1[0] = c; row1[1] = d;
   blk->InsertGlobalValues(0,2,&row0[0],&indicies[0]);
   blk->InsertGlobalValues(1,2,&row1[0],&indicies[0]);
   blk->FillComplete();

   return Thyra::epetraLinearOp(blk);
}

TEUCHOS_UNIT_TEST(tLU2x2InverseOp, exact_test)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   Teko::LinearOp  A_00 = build2x2(Comm,1,2,3,4);
   Teko::LinearOp  A_01 = build2x2(Comm,5,6,7,8);
   Teko::LinearOp  A_10 = build2x2(Comm,9,10,11,12);
   Teko::LinearOp  A_11 = build2x2(Comm,-13,-14,-15,-16);
   Teko::BlockedLinearOp A = Teko::toBlockedLinearOp(Thyra::block2x2(A_00,A_01,A_10,A_11));

   Teko::LinearOp  S = build2x2(Comm,26.000000000000000, 28.000000000000004, 30.000000000000000, 32.000000000000000);

   Teko::LinearOp iA_00 = build2x2(Comm,-1.000000000000000,  0.500000000000000,  0.749999999999998, -0.249999999999998);
   Teko::LinearOp iA_01 = build2x2(Comm,-3.000000000000004,  2.500000000000003,  2.750000000000008, -2.250000000000007);
   Teko::LinearOp iA_10 = build2x2(Comm,-1.999999999999999,  1.499999999999999,  1.749999999999999, -1.249999999999999);
   Teko::LinearOp iA_11 = build2x2(Comm, 4.000000000000001, -3.500000000000001, -3.750000000000001,  3.250000000000001);
   Teko::LinearOp iA = Thyra::block2x2(iA_00,iA_01,iA_10,iA_11);

   Thyra::LinearOpTester<double> tester;
   tester.dump_all(true);
   tester.show_all_tests(true);

   {
      RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromStratimikos();
      RCP<Teko::InverseFactory> invFact = invLib->getInverseFactory("Amesos");

      Teko::LinearOp invA_00 = Teko::buildInverse(*invFact,A_00);
      Teko::LinearOp invS    = Teko::buildInverse(*invFact,S);

      Teko::LinearOp invA = Teko::createLU2x2InverseOp(A,invA_00,invA_00,invS,"Approximation");
      
      const bool result = tester.compare( *invA, *iA, &out);
      if (!result) {
         out << "Apply: FAILURE" << std::endl;
         success = false;
      }
      else
         out << "Apply: SUCCESS" << std::endl;
   }
}
