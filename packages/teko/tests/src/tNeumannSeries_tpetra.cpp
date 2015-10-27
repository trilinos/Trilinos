/*
// @HEADER
// 
// ***********************************************************************
// 
//      Teko: A package for block and physics based preconditioning
//                  Copyright 2010 Sandia Corporation 
//  
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//  
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//  
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//  
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//  
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission. 
//  
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING 
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//  
// Questions? Contact Eric C. Cyr (eccyr@sandia.gov)
// 
// ***********************************************************************
// 
// @HEADER

*/

#include "tNeumannSeries_tpetra.hpp"

#include <string>

// Teko-Package includes
#include "Teko_Utilities.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_NeumannSeriesPreconditionerFactory.hpp"

// Thyra includes
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"

// TriUtils includes
#include "Trilinos_Util_CrsMatrixGallery.h"

// Test-rig
#include "Test_Utils.hpp"

namespace Teko {
namespace Test {

using Teuchos::rcp;
using Teuchos::RCP;
using Teuchos::rcpFromRef;
using Thyra::tpetraLinearOp;

static RCP<Teuchos::ParameterList> buildLibPL(int numTerms,std::string scalingType)
{
   RCP<Teuchos::ParameterList> pl = rcp(new Teuchos::ParameterList());
   RCP<Teuchos::ParameterList> nList = rcp(new Teuchos::ParameterList());

   nList->set("Type","Neumann Series");
   nList->set("Number of Terms",numTerms);
   nList->set("Scaling Type",scalingType);

   pl->set("Neumann",*nList);

   return pl;
}

RCP<Thyra::LinearOpBase<ST> > buildExampleOp(int type,RCP<const Teuchos::Comm<int> > comm)
{
   const RCP<const Tpetra::Map<LO,GO,NT> > map = rcp(new const Tpetra::Map<LO,GO,NT>(3,0,comm));
   const RCP<Tpetra::CrsMatrix<ST,LO,GO,NT> > mat  = Tpetra::createCrsMatrix<ST,LO,GO,NT>(map, 3);

   ST values[3];
   GO indices[3] = { 0, 1, 2 };

   switch(type) {
   case 2:
      values[0] = 5.0; values[1] = 0.0; values[2] = 0.0;
      mat->insertGlobalValues(0,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));
   
      values[0] = 2.0; values[1] = 6.0; values[2] = 0.0;
      mat->insertGlobalValues(1,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));
   
      values[0] = 0.0; values[1] = 3.0; values[2] = 7.0;
      mat->insertGlobalValues(2,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));
      break;
   case 1:
      values[0] = 1.0; values[1] = 0.0; values[2] = 0.0;
      mat->insertGlobalValues(0,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));
   
      values[0] = 2.0; values[1] = 1.0; values[2] = 0.0;
      mat->insertGlobalValues(1,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));
   
      values[0] = 0.0; values[1] = 3.0; values[2] = 1.0;
      mat->insertGlobalValues(2,Teuchos::ArrayView<GO>(indices,3),Teuchos::ArrayView<ST>(values,3));
   default:
      break;
   }

   mat->fillComplete();

   return Thyra::tpetraLinearOp<ST,LO,GO,NT>(Thyra::tpetraVectorSpace<ST,LO,GO,NT>(mat->getRangeMap()),Thyra::tpetraVectorSpace<ST,LO,GO,NT>(mat->getDomainMap()),mat);
}

void tNeumannSeries_tpetra::initializeTest()
{
   tolerance_ = 1.0e-15;
}

int tNeumannSeries_tpetra::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tNeumannSeries_tpetra";

   status = test_simpleOp(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"test_simpleOp\" ... PASSED","   \"test_simpleOp\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = test_scaledOp(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"test_scaledOp\" ... PASSED","   \"test_scaledOp\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tNeumannSeries_tpetra...PASSED","tNeumannSeries_tpetra...FAILED");
   }
   else {// Normal Operating Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tNeumannSeries_tpetra...FAILED");
   }

   return failcount;
}

bool tNeumannSeries_tpetra::test_simpleOp(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;


   // perform actual test
   Thyra::LinearOpTester<ST> tester;
   tester.show_all_tests(true);
   std::stringstream ss;
   Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");

   {
      // build library parameter list
      Teuchos::RCP<Teuchos::ParameterList> pl = buildLibPL(1,"None");
      if(verbosity>=10) {
         os << "   tNeumannSeries_tpetra::test_simpleOp :"
            << " printing library parameter list" << std::endl;
         pl->print(os);
      }

      RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(*pl);
      RCP<Teko::InverseFactory> neumann = invLib->getInverseFactory("Neumann");
      RCP<Teko::InverseFactory> direct = invLib->getInverseFactory("Ifpack2");

      Teko::LinearOp op = buildExampleOp(1,GetComm_tpetra());
   
      Teko::LinearOp neuInv = Teko::buildInverse(*neumann,op);
      Teko::LinearOp dirInv = Teko::buildInverse(*direct,op);

      const bool result = tester.compare( *neuInv, *dirInv, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(not result,
             std::endl << "   tNeumannSeries_tpetra::test_simpleOp "
             << ": Comparing underresolved factory generated operator to correct operator");
      if(result || verbosity>=10) 
         os << ss.str(); 
   }

   {
      // build library parameter list
      Teuchos::RCP<Teuchos::ParameterList> pl = buildLibPL(2,"None");
      if(verbosity>=10) {
         os << "   tNeumannSeries_tpetra::test_simpleOp :"
            << " printing library parameter list" << std::endl;
         pl->print(os);
      }

      RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(*pl);
      RCP<Teko::InverseFactory> neumann = invLib->getInverseFactory("Neumann");
      RCP<Teko::InverseFactory> direct = invLib->getInverseFactory("Ifpack2");

      Teko::LinearOp op = buildExampleOp(1,GetComm_tpetra());
   
      Teko::LinearOp neuInv = Teko::buildInverse(*neumann,op);
      Teko::LinearOp dirInv = Teko::buildInverse(*direct,op);

      const bool result = tester.compare( *neuInv, *dirInv, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tNeumannSeries_tpetra::test_simpleOp "
             << ": Comparing factory generated operator to correct operator");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }
 
   return allPassed;
}

bool tNeumannSeries_tpetra::test_scaledOp(int verbosity,std::ostream & os)
{
   bool status = false;
   bool allPassed = true;

   // perform actual test
   Thyra::LinearOpTester<ST> tester;
   tester.show_all_tests(true);
   std::stringstream ss;
   Teuchos::FancyOStream fos(rcpFromRef(ss),"      |||");

   {
      // build library parameter list
      Teuchos::RCP<Teuchos::ParameterList> pl = buildLibPL(2,"Diagonal");
      if(verbosity>=10) {
         os << "   tNeumannSeries_tpetra::test_scaledOp :"
            << " printing library parameter list" << std::endl;
         pl->print(os);
      }

      RCP<Teko::InverseLibrary> invLib = Teko::InverseLibrary::buildFromParameterList(*pl);
      RCP<Teko::InverseFactory> neumann = invLib->getInverseFactory("Neumann");
      RCP<Teko::InverseFactory> direct = invLib->getInverseFactory("Ifpack2");

      Teko::LinearOp op = buildExampleOp(2,GetComm_tpetra());
   
      Teko::LinearOp neuInv = Teko::buildInverse(*neumann,op);
      Teko::LinearOp dirInv = Teko::buildInverse(*direct,op);

      const bool result = tester.compare( *neuInv, *dirInv, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
             std::endl << "   tNeumannSeries_tpetra::test_scaledOp "
             << ": Comparing factory generated operator to correct operator");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }
 
   return allPassed;
}

} // end namespace Tests
} // end namespace Teko
