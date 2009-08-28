#include "tLSCHIntegrationTest.hpp"

#include <string>

// Epetra includes
#include "Epetra_Export.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// ML includes
#include "ml_include.h"
#include "ml_MultiLevelPreconditioner.h"

// PB-Package includes
#include "Epetra/PB_EpetraHelpers.hpp"
#include "Epetra/PB_EpetraBlockPreconditioner.hpp"
#include "NS/PB_LSCPreconditionerFactory.hpp"
#include "NS/PB_InvLSCStrategy.hpp"
#include "PB_Utilities.hpp"

#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

// Aztec includes
#include "AztecOO.h"
#include "AztecOO_Operator.h"

// Stratimikos includes
#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

// Test-rig
#include "Test_Utils.hpp"

#include "PB_InverseFactory.hpp"
#include "PB_Utilities.hpp"

namespace PB {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Thyra::epetraLinearOp;

void tLSCHIntegrationTest::initializeTest()
{
   tolerance_ = 1.0e-7;
}

int tLSCHIntegrationTest::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tLSCHIntegrationTest";

   status = test_hScaling(verbosity,failstrm);
   PB_TEST_MSG(stdstrm,1,"   \"hScaling\" ... PASSED","   \"hScaling\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      PB_TEST_MSG(failstrm,0,"tLSCHIntegrationTest...PASSED","tLSCHIntegrationTest...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      PB_TEST_MSG(failstrm,0,"...PASSED","tLSCHIntegrationTest...FAILED");
   }

   return failcount;
}

bool tLSCHIntegrationTest::test_hScaling(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;
   
   RCP<const Epetra_Comm> comm = GetComm();

   // build some operators
   PB::LinearOp F = PB::Test::build2x2(*comm,1,2,2,1);
   PB::LinearOp G = PB::Test::build2x2(*comm,1,-1,-3,1);
   PB::LinearOp D = PB::Test::build2x2(*comm,1,-3,-1,1);

   double diag[2];
 
   diag[0] = 1.0/3.0; diag[1] = 1.0/2.0;
   PB::LinearOp M = PB::Test::DiagMatrix(2,diag,"M");

   diag[0] = 5.0; diag[1] = 9.0;
   PB::LinearOp H = PB::Test::DiagMatrix(2,diag,"H");

   PB::LinearOp A = Thyra::block2x2<double>(F,G,D,Teuchos::null);

   PB::LinearOp exact;
   {
      // build some operators
      PB::LinearOp D0 = PB::Test::build2x2(*comm,-1.0/3.0,2.0/3.0,2.0/3.0,-1.0/3.0);
      PB::LinearOp D1 = PB::Test::build2x2(*comm,-1.5,-3.0,-3.0,-5.5);
      PB::LinearOp U  = PB::Test::build2x2(*comm,-0.5,-1.5,-0.5,-0.5);
      
      exact = Thyra::block2x2<double>(D0,U,Teuchos::null,D1);
   }

   RCP<PB::InverseLibrary> invLib = PB::InverseLibrary::buildFromStratimikos();
   RCP<const PB::InverseFactory> invFact = invLib->getInverseFactory("Amesos");
   RCP<PB::NS::InvLSCStrategy> strategy = rcp(new PB::NS::InvLSCStrategy(invFact,M));
   strategy->setHScaling(PB::getDiagonal(H));
   strategy->setUseFullLDU(false);

   RCP<PB::BlockPreconditionerFactory> precFact = rcp(new PB::NS::LSCPreconditionerFactory(strategy));
   RCP<PB::BlockPreconditionerState> bps = precFact->buildPreconditionerState();
   PB::LinearOp prec = precFact->buildPreconditionerOperator(A,*bps);

   PB::BlockedLinearOp bA = PB::toBlockedLinearOp(A);
   std::stringstream ss;
   ss << "invF = " << Teuchos::describe(*strategy->getInvF(bA,*bps),Teuchos::VERB_EXTREME) << std::endl;
   ss << "invBQBt = " << Teuchos::describe(*strategy->getInvBQBt(bA,*bps),Teuchos::VERB_EXTREME) << std::endl;
   ss << "invF = " << Teuchos::describe(*strategy->getInvBHBt(bA,*bps),Teuchos::VERB_EXTREME) << std::endl;
   ss << "invMass = " << Teuchos::describe(*strategy->getInvMass(bA,*bps),Teuchos::VERB_EXTREME) << std::endl;
   ss << "HScaling = " << Teuchos::describe(*strategy->getHScaling(bA,*bps),Teuchos::VERB_EXTREME) << std::endl;
   ss << "prec = " << Teuchos::describe(*prec,Teuchos::VERB_EXTREME) << std::endl;

   // construct a couple of vectors
   const RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,*GetComm()));
   Epetra_Vector ea(*map),eb(*map);
   const RCP<const Thyra::VectorBase<double> > x = BlockVector(ea,eb,prec->domain());
   const RCP<Thyra::VectorBase<double> > y = Thyra::createMember(prec->range()); 
   ea[0] = 1.0; ea[1] = 0.0; eb[0] = 0.0; eb[1] = 0.0;
   Thyra::apply(*prec,Thyra::NONCONJ_ELE,*x,&*y);
   ss << "prec = " << Teuchos::describe(*y,Teuchos::VERB_EXTREME) << std::endl;
   ea[0] = 0.0; ea[1] = 1.0; eb[0] = 0.0; eb[1] = 0.0;
   Thyra::apply(*prec,Thyra::NONCONJ_ELE,*x,&*y);
   ss << "prec = " << Teuchos::describe(*y,Teuchos::VERB_EXTREME) << std::endl;
   ea[0] = 0.0; ea[1] = 0.0; eb[0] = 1.0; eb[1] = 0.0;
   Thyra::apply(*prec,Thyra::NONCONJ_ELE,*x,&*y);
   ss << "prec = " << Teuchos::describe(*y,Teuchos::VERB_EXTREME) << std::endl;
   ea[0] = 0.0; ea[1] = 0.0; eb[0] = 0.0; eb[1] = 1.0;
   Thyra::apply(*prec,Thyra::NONCONJ_ELE,*x,&*y);
   ss << "prec = " << Teuchos::describe(*y,Teuchos::VERB_EXTREME) << std::endl;
  
   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);
   Teuchos::FancyOStream fos(Teuchos::rcpFromRef(ss),"|||");
   const bool result = tester.compare( *prec, *exact, &fos );
   TEST_ASSERT(result,
          std::endl << "   tLSCHIntegration::test_hScaling "
                    << ": Comparing precondtioner to exactly computed version");

   const bool result2 = tester.compare( *H, *strategy->getHScaling(bA,*bps), &fos );
   TEST_ASSERT(result2,
          std::endl << "   tLSCHIntegration::test_hScaling "
                    << ": Comparing scaling of H operator");

   const bool result3 = tester.compare( *build2x2(*comm,0.208333333333333,0.375, 0.375, 0.875),
                          *strategy->getInvBQBt(bA,*bps), &fos );
   TEST_ASSERT(result3,
          std::endl << "   tLSCHIntegration::test_hScaling "
                    << ": Comparing inv(BQBt) operator");

   const bool result4 = tester.compare( *build2x2(*comm, 0.077777777777778, 0.177777777777778, 0.177777777777778, 0.477777777777778),
                          *strategy->getInvBHBt(bA,*bps), &fos );
   TEST_ASSERT(result4,
          std::endl << "   tLSCHIntegration::test_hScaling "
                    << ": Comparing inv(BHBt) operator");

   if(not allPassed || verbosity>=10) 
      os << ss.str(); 

   return allPassed;
}

} // end namespace Tests
} // end namespace PB
