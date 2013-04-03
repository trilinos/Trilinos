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

#include "tLSCStabilized.hpp"
#include "Teko_LSCPreconditionerFactory.hpp"
#include "Teko_InvLSCStrategy.hpp"
#include "Teko_InverseLibrary.hpp"
#include "Teko_InverseFactory.hpp"

// Teuchos includes
#include "Teuchos_RCP.hpp"

// Epetra includes
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "Thyra_DefaultZeroLinearOp.hpp"
#include "Thyra_DefaultLinearOpSource.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_VectorStdOps.hpp"

#include <vector>

// This whole test rig is based on inverting the matrix
// 
//      [  1  2  1 -1 ]
//  A = [  2  1 -3  1 ]
//      [  1 -3  0  0 ]
//      [ -1  1  0  0 ]
//
// see the matlab file

namespace Teko {
namespace Test {

using namespace Teuchos;
using namespace Thyra;
using namespace Teko::NS;

void tLSCStabilized::initializeTest()
{
   tolerance_ = 1.0e-13;

   comm = rcp(new Epetra_SerialComm());
}

int tLSCStabilized::runTest(int verbosity,std::ostream & stdstrm,std::ostream & failstrm,int & totalrun)
{
   bool allTests = true;
   bool status;
   int failcount = 0;

   failstrm << "tLSCStabilized";

   status = test_diagonal(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"diagonal\" ... PASSED","   \"diagonal\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;

/*
   status = test_strategy(verbosity,failstrm);
   Teko_TEST_MSG(stdstrm,1,"   \"strategy\" ... PASSED","   \"strategy\" ... FAILED");
   allTests &= status;
   failcount += status ? 0 : 1;
   totalrun++;
*/

   status = allTests;
   if(verbosity >= 10) {
      Teko_TEST_MSG(failstrm,0,"tLSCStabilized...PASSED","tLSCStablilized...FAILED");
   }
   else {// Normal Operatoring Procedures (NOP)
      Teko_TEST_MSG(failstrm,0,"...PASSED","tLSCStabilized...FAILED");
   }

   return failcount;
}

bool tLSCStabilized::test_diagonal(int verbosity,std::ostream & os)
{
   // make sure the preconditioner is working by testing against the identity matrix
   typedef RCP<const Thyra::VectorBase<double> > Vector;
   typedef RCP<const Thyra::VectorSpaceBase<double> > VectorSpace;
   typedef RCP<const Thyra::LinearOpBase<double> > LinearOp;

   bool status = false;
   bool allPassed = true;
   double vec[2];
   double diff = 0.0;

   // build 4x4 matrix with block 2x2 diagonal subblocks
   //
   //            [ 1 0 7 0 ]
   // [ F G ] =  [ 0 2 0 8 ]
   // [ D C ]    [ 5 0 3 0 ]
   //            [ 0 6 0 4 ]
   //

   vec[0] = 1.0; vec[1] = 2.0;
   LinearOp F = Teko::Test::DiagMatrix(2,vec);

   vec[0] = 7.0; vec[1] = 8.0;
   LinearOp G = Teko::Test::DiagMatrix(2,vec);

   vec[0] = 5.0; vec[1] = 6.0;
   LinearOp D = Teko::Test::DiagMatrix(2,vec);

   vec[0] = 3.0; vec[1] = 4.0;
   LinearOp C = Teko::Test::DiagMatrix(2,vec);

   vec[0] = 1.0; vec[1] = 0.5;
   LinearOp iF = Teko::Test::DiagMatrix(2,vec);

   vec[0] = 0.030303030303030; vec[1] = 0.02205882352941;
   LinearOp iBBt = Teko::Test::DiagMatrix(2,vec);

   vec[0] = 0.026041666666667; vec[1] = 0.041666666666667;
   LinearOp aiD = Teko::Test::DiagMatrix(2,vec);

   LinearOp A = Thyra::block2x2(F,G,D,C);
 
   const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory 
         = rcp(new LSCPreconditionerFactory(iF,iBBt,aiD,Teuchos::null));
   RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory,A);

   // build linear operator
   RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp(); 

   const RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,*comm));
   // construct a couple of vectors
   Epetra_Vector ea(*map),eb(*map);
   Epetra_Vector ef(*map),eg(*map);
   const RCP<const Thyra::MultiVectorBase<double> > x = BlockVector(ea,eb,A->domain());
   const RCP<const Thyra::MultiVectorBase<double> > z = BlockVector(ef,eg,A->domain());
   const RCP<Thyra::MultiVectorBase<double> > y = Thyra::createMembers(A->range(),1); 

   // now checks of the preconditioner (should be exact!)
   /////////////////////////////////////////////////////////////////////////

   // test vector [0 1 1 3]
   ea[0] = 0.0; ea[1] = 1.0; eb[0] = 1.0; eb[1] = 3.0;
   ef[0] =  0.407268709825528; ef[1] =  1.560553633217993;
   eg[0] = -0.058181244260790; eg[1] = -0.265138408304498;
   Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
   TEST_ASSERT((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_,
               "   tLSCStabilized::test_diagonal " << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
            << diff << " <= " << tolerance_ << ")\n"
            << "      " << Print("x",x) 
            << "      " << Print("y",y) 
            << "      " << Print("z",z));

   // test vector [-2 4 7 9]
   ea[0] =-2.0; ea[1] = 4.0; eb[0] = 7.0; eb[1] = 9.0;
   ef[0] =  0.850880968778696; ef[1] =  5.181660899653979;
   eg[0] = -0.407268709825528; eg[1] = -0.795415224913495;
   Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
   TEST_ASSERT((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_,
               "   tLSCStabilized::test_diagonal " << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
            << diff << " <= " << tolerance_ << ")\n"
            << "      " << Print("x",x) 
            << "      " << Print("y",y) 
            << "      " << Print("z",z));

   // test vector [1 0 0 -5]
   ea[0] = 1.0; ea[1] = 0.0; eb[0] = 0.0; eb[1] =-5.0;
   ef[0] =  1.000000000000000; ef[1] = -1.767589388696655;
   eg[0] =  0.000000000000000; eg[1] =  0.441897347174164;
   Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
   TEST_ASSERT((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_,
               "   tLSCStabilized::test_diagonal " << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
            << diff << " <= " << tolerance_ << ")\n"
            << "      " << Print("x",x) 
            << "      " << Print("y",y) 
            << "      " << Print("z",z));

   // test vector [4 -4 6 12]
   ea[0] = 4.0; ea[1] =-4.0; eb[0] = 6.0; eb[1] =12.0;
   ef[0] =  6.443612258953168; ef[1] =  2.242214532871971;
   eg[0] = -0.349087465564738; eg[1] = -1.060553633217993;
   Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
   TEST_ASSERT((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_,
               "   tLSCStabilized::test_diagonal " << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
            << diff << " <= " << tolerance_ << ")\n"
            << "      " << Print("x",x) 
            << "      " << Print("y",y) 
            << "      " << Print("z",z));

   return allPassed;
}

bool tLSCStabilized::test_diagonalNotSym(int verbosity,std::ostream & os)
{
   // make sure the preconditioner is working by testing against the identity matrix
   typedef RCP<const Thyra::VectorBase<double> > Vector;
   typedef RCP<const Thyra::VectorSpaceBase<double> > VectorSpace;
   typedef RCP<const Thyra::LinearOpBase<double> > LinearOp;

   bool status = false;
   bool allPassed = true;
   double vec[2];
   double diff = 0.0;

   // build 4x4 matrix with block 2x2 diagonal subblocks
   //
   //            [ 1 0 7 0 ]
   // [ F G ] =  [ 0 2 0 8 ]
   // [ D C ]    [ 5 0 3 0 ]
   //            [ 0 6 0 4 ]
   //

   vec[0] = 1.0; vec[1] = 2.0;
   LinearOp F = Teko::Test::DiagMatrix(2,vec);

   vec[0] = 7.0; vec[1] = 8.0;
   LinearOp G = Teko::Test::DiagMatrix(2,vec);

   vec[0] = 5.0; vec[1] = 6.0;
   LinearOp D = Teko::Test::DiagMatrix(2,vec);

   vec[0] = 3.0; vec[1] = 4.0;
   LinearOp C = Teko::Test::DiagMatrix(2,vec);

   vec[0] = 1.0; vec[1] = 0.5;
   LinearOp iF = Teko::Test::DiagMatrix(2,vec);

   vec[0] = 0.030303030303030; vec[1] = 0.02205882352941;
   LinearOp iBBt = Teko::Test::DiagMatrix(2,vec);

   vec[0] = 0.026041666666667; vec[1] = 0.041666666666667;
   LinearOp aiD = Teko::Test::DiagMatrix(2,vec);

   LinearOp A = Thyra::block2x2(F,G,D,C);

   RCP<InverseLibrary> invLib = InverseLibrary::buildFromStratimikos();
   RCP<InverseFactory> invFact = invLib->getInverseFactory("Amesos");
 
   RCP<InvLSCStrategy> lscStrat = rcp(new InvLSCStrategy(invFact));
   // lscStrat->setSymmetric(false);
   const RCP<const Thyra::PreconditionerFactoryBase<double> > precFactory 
         = rcp(new LSCPreconditionerFactory(lscStrat));
   RCP<Thyra::PreconditionerBase<double> > prec = Thyra::prec<double>(*precFactory,A);

   // build linear operator
   RCP<const Thyra::LinearOpBase<double> > precOp = prec->getUnspecifiedPrecOp(); 

   const RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,*comm));
   // construct a couple of vectors
   Epetra_Vector ea(*map),eb(*map);
   Epetra_Vector ef(*map),eg(*map);
   const RCP<const Thyra::MultiVectorBase<double> > x = BlockVector(ea,eb,A->domain());
   const RCP<const Thyra::MultiVectorBase<double> > z = BlockVector(ef,eg,A->domain());
   const RCP<Thyra::MultiVectorBase<double> > y = Thyra::createMembers(A->range(),1); 

   // now checks of the preconditioner (should be exact!)
   /////////////////////////////////////////////////////////////////////////

   // test vector [0 1 1 3]
   ea[0] = 0.0; ea[1] = 1.0; eb[0] = 1.0; eb[1] = 3.0;
   ef[0] =  0.407268709825528; ef[1] =  1.560553633217993;
   eg[0] = -0.058181244260790; eg[1] = -0.265138408304498;
   Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
   TEST_ASSERT((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_,
               "   tLSCStabilized::test_diagonal " << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
            << diff << " <= " << tolerance_ << ")\n"
            << "      " << Print("x",x) 
            << "      " << Print("y",y) 
            << "      " << Print("z",z));

   // test vector [-2 4 7 9]
   ea[0] =-2.0; ea[1] = 4.0; eb[0] = 7.0; eb[1] = 9.0;
   ef[0] =  0.850880968778696; ef[1] =  5.181660899653979;
   eg[0] = -0.407268709825528; eg[1] = -0.795415224913495;
   Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
   TEST_ASSERT((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_,
               "   tLSCStabilized::test_diagonal " << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
            << diff << " <= " << tolerance_ << ")\n"
            << "      " << Print("x",x) 
            << "      " << Print("y",y) 
            << "      " << Print("z",z));

   // test vector [1 0 0 -5]
   ea[0] = 1.0; ea[1] = 0.0; eb[0] = 0.0; eb[1] =-5.0;
   ef[0] =  1.000000000000000; ef[1] = -1.767589388696655;
   eg[0] =  0.000000000000000; eg[1] =  0.441897347174164;
   Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
   TEST_ASSERT((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_,
               "   tLSCStabilized::test_diagonal " << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
            << diff << " <= " << tolerance_ << ")\n"
            << "      " << Print("x",x) 
            << "      " << Print("y",y) 
            << "      " << Print("z",z));

   // test vector [4 -4 6 12]
   ea[0] = 4.0; ea[1] =-4.0; eb[0] = 6.0; eb[1] =12.0;
   ef[0] =  6.443612258953168; ef[1] =  2.242214532871971;
   eg[0] = -0.349087465564738; eg[1] = -1.060553633217993;
   Thyra::apply(*precOp,Thyra::NOTRANS,*x,y.ptr());
   TEST_ASSERT((diff = Teko::Test::Difference(y,z)/Thyra::norm_2(*z->col(0)))<tolerance_,
               "   tLSCStabilized::test_diagonal " << toString(status) << ":(y=inv(A)*x) != z (|y-z|_2/|z|_2 = "
            << diff << " <= " << tolerance_ << ")\n"
            << "      " << Print("x",x) 
            << "      " << Print("y",y) 
            << "      " << Print("z",z));

   return allPassed;
}

bool tLSCStabilized::test_strategy(int verbosity,std::ostream & os)
{
   std::vector<int> indicies(2);
   std::vector<double> row0(2);
   int sz = 5;
   double vec[5];

   bool status = false;
   bool allPassed = true;

   vec[0] = 1.0; vec[1] = 2.0; vec[2] = 3.0; vec[3] = 4.0; vec[4] = 5.0;
   LinearOp F = Teko::Test::DiagMatrix(sz,vec);

   vec[0] = 7.0; vec[1] = 8.0; vec[2] = 9.0; vec[3] = 10.0; vec[4] = 11.0;
   LinearOp G = Teko::Test::DiagMatrix(sz,vec);

   vec[0] = 5.0; vec[1] = 6.0; vec[2] = 7.0; vec[3] = 8.0; vec[4] = 9.0;
   LinearOp D = Teko::Test::DiagMatrix(sz,vec);

   vec[0] = 3.0; vec[1] = 4.0; vec[2] = 5.0; vec[3] = 6.0; vec[4] = 7.0;
   LinearOp C = Teko::Test::DiagMatrix(sz,vec);

   vec[0] = 1.0; vec[1] = 1.0/2.0; vec[2] = 1.0/3.0; vec[3] = 1.0/4.0; vec[4] = 1.0/5.0;
   LinearOp iF = Teko::Test::DiagMatrix(sz,vec);

   vec[0] = 0.091304347826087;
   vec[1] = 0.090517241379310;
   vec[2] = 0.087646076794658;
   vec[3] = 0.084000000000000;
   vec[4] = 0.080152671755725;
   LinearOp iBQBtmC = Teko::Test::DiagMatrix(sz,vec);

   vec[0] = 0.020202020202020;
   vec[1] = 0.032323232323232;
   vec[2] = 0.040404040404040;
   vec[3] = 0.046176046176046;
   vec[4] = 0.050505050505051;
   LinearOp aiD = Teko::Test::DiagMatrix(sz,vec);

   LinearOp A = Thyra::block2x2(F,G,D,C);

   comm = rcp(new Epetra_SerialComm());
   const RCP<Epetra_Map> map = rcp(new Epetra_Map(sz,0,*comm));

   Epetra_Vector ea(*map),eb(*map);
   const RCP<const Thyra::MultiVectorBase<double> > x = BlockVector(ea,eb,A->domain());
   const RCP<Thyra::MultiVectorBase<double> > y = Thyra::createMembers(A->range(),1); 

   Teuchos::ParameterList paramList;
   paramList.set("Linear Solver Type","Amesos");
   RCP<Teko::InverseFactory> invFact = Teko::invFactoryFromParamList(paramList,"Amesos");

   // build Mass matrix
   vec[0] = 3.0; vec[1] = 4.0; vec[2] = 5.0; vec[3] = 6.0; vec[4] = 7.0;
   LinearOp mass = Teko::Test::DiagMatrix(sz,vec);

   vec[0] = 1.0/3.0; vec[1] = 1.0/4.0; vec[2] = 1.0/5.0; vec[3] = 1.0/6.0; vec[4] = 1.0/7.0;
   LinearOp invMass = Teko::Test::DiagMatrix(sz,vec);

   Thyra::LinearOpTester<double> tester;
   tester.set_all_error_tol(5e-4);
   tester.show_all_tests(true);
   std::stringstream ss;
   Teuchos::FancyOStream fos(Teuchos::rcpFromRef(ss),"      |||");

   Teko::BlockedLinearOp blkA = Teko::toBlockedLinearOp(A);

   // build preconditioner
   vec[0] = 1.0; vec[1] = 0.5; vec[2] = 1.0/3.0; vec[3] = 0.25; vec[4] = 0.2;
   LinearOp p00 = Teko::Test::DiagMatrix(sz,vec);

   vec[0] = 0.368351759561589;
   vec[1] = 0.325933832979017;
   vec[2] = 0.295436133965709;
   vec[3] = 0.272240115440115;
   vec[4] = 0.253891252128534;
   LinearOp p01 = Teko::Test::DiagMatrix(sz,vec);

   vec[0] = 0; vec[1] = 0; vec[2] = 0; vec[3] = 0; vec[4] = 0;
   LinearOp p10 = Teko::Test::DiagMatrix(sz,vec);

   vec[0] = -0.052621679937370;
   vec[1] = -0.081483458244754;
   vec[2] = -0.098478711321903;
   vec[3] = -0.108896046176046;
   vec[4] = -0.115405114603879;
   LinearOp p11 = Teko::Test::DiagMatrix(sz,vec);
   LinearOp P = Thyra::block2x2(p00,p01,p10,p11);
  
   // Kluge to get around problem with Anasazi
   // Teko::computeSpectralRad(Thyra::multiply(invMass,F),5e-2,false,3)/3.0;
   // Teko::computeSpectralRad(Thyra::multiply(invMass,F),5e-2,false,3)/3.0;
             
   // build inverse strategy
   { 
      bool result;
      Teko::NS::LSCPrecondState state;
      Teko::NS::InvLSCStrategy iStrat(invFact,mass,false);
      iStrat.setEigSolveParam(3);
      Teko::NS::LSCPreconditionerFactory factory(Teuchos::rcpFromRef(iStrat));
      LinearOp prec = factory.buildPreconditionerOperator(blkA,state);

      // test inverse mass
      ss.str("");
      result = tester.compare( *invMass, *iStrat.getInvMass(blkA,state), Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
            std::endl << "   tLSCStabilized::test_strategy " << toString(status)
                      << " : Comparing mass operators");
      if(not result || verbosity>=10) 
         os << ss.str(); 

      // test inverse F
      ss.str("");
      result = tester.compare( *iF, *iStrat.getInvF(blkA,state), Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
            std::endl << "   tLSCStabilized::test_strategy " << toString(status)
                      << " : Comparing F operators");
      if(not result || verbosity>=10) 
         os << ss.str(); 

      // test inverse B*Q*Bt-gamma*C
      ss.str("");
      result = tester.compare( *iBQBtmC, *iStrat.getInvBQBt(blkA,state), Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
            std::endl << "   tLSCStabilized::test_strategy " << toString(status)
                      << " : Comparing B*Q*Bt-C operators");
      if(not result || verbosity>=10) 
         os << ss.str(); 

      // test alpha*inv(D)
      ss.str("");
      // result = tester.compare( *aiD, *iStrat.getInvAlphaD(blkA,state), &fos );
      result = tester.compare( *aiD, *iStrat.getOuterStabilization(blkA,state), Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
            std::endl << "   tLSCStabilized::test_strategy " << toString(status)
                      << " : Comparing alpha*inv(D) operators");
      if(not result || verbosity>=10) 
         os << ss.str(); 

      // test full op
      ss.str("");
      result = tester.compare( *P, *prec, Teuchos::ptrFromRef(fos) );
      TEST_ASSERT(result,
            std::endl << "   tLSCStabilized::test_strategy " << toString(status)
                      << " : Comparing full op");
      if(not result || verbosity>=10) 
         os << ss.str(); 
   }

   return allPassed;
}

} // end namespace Test
} // end namespace Teko
