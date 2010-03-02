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
#include "Teko_RepeatPreconditionerFactory.hpp"

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

RCP<Teuchos::ParameterList> buildLibPL(int numTerms,std::string scalingType)
{
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());
   RCP<Teuchos::ParameterList> nList = rcp(new Teuchos::ParameterList());

   nList->set("Type","Neumann Series");
   nList->set("Number of Terms",numTerms);
   nList->set("Scaling Type",scalingType);

   pl->set("Neumann",*nList);

   return pl;
}

TEUCHOS_UNIT_TEST(tRepeatPreconditionerFactory, constructor_test)
{
   // build global (or serial communicator)
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif

   Teko::LinearOp  A = build2x2(Comm,1,2,3,4);
   Teko::LinearOp iP = build2x2(Comm,1.0,0.0,0.0,1.0/4.0);
   Teko::LinearOp ImAiP = build2x2(Comm,0.0,-0.5,-3.0,0.0);

   Thyra::LinearOpTester<double> tester;
   tester.show_all_tests(true);

   {
      RCP<Teko::InverseFactory> precOpFact = rcp(new Teko::StaticOpInverseFactory(iP));
      RCP<Teko::RepeatPreconditionerFactory> precFact = rcp(new Teko::RepeatPreconditionerFactory(0,precOpFact));
      RCP<Teko::InverseFactory> invFact = rcp(new Teko::PreconditionerInverseFactory(precFact));

      Teko::LinearOp prec = Teko::buildInverse(*invFact,A);

      const bool result = tester.compare( *prec, *iP, &out);
      if (!result) {
         out << "Apply 0: FAILURE" << std::endl;
         success = false;
      }
      else
         out << "Apply 0: SUCCESS" << std::endl;
   }

   {
      using Teko::multiply;

      RCP<Teko::InverseFactory> precOpFact = rcp(new Teko::StaticOpInverseFactory(iP));
      RCP<Teko::RepeatPreconditionerFactory> precFact = rcp(new Teko::RepeatPreconditionerFactory(3,precOpFact));
      RCP<Teko::InverseFactory> invFact = rcp(new Teko::PreconditionerInverseFactory(precFact));

      Teko::LinearOp prec = Teko::buildInverse(*invFact,A);
      Teko::LinearOp exact = multiply(multiply(iP,ImAiP,ImAiP),ImAiP);

      const bool result = tester.compare( *prec, *iP, &out);
      if (!result) {
         out << "Apply 3: FAILURE" << std::endl;
         success = false;
      }
      else
         out << "Apply 3: SUCCESS" << std::endl;
   }
}
