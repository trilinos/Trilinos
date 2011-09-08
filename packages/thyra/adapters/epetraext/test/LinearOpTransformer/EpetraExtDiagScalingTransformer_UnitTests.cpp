/*
// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
//                 Copyright (2004) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#include "Thyra_EpetraExtDiagScalingTransformer.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"
#ifdef _MSC_VER
// This is required for operator not to be defined
#include <iso646.h>
#endif

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif

#include "Teuchos_UnitTestHarness.hpp"


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Thyra::EpetraExtDiagScalingTransformer;
using Thyra::epetraExtDiagScalingTransformer;
using Thyra::VectorBase;
using Thyra::LinearOpBase;
using Thyra::createMember;
using Thyra::LinearOpTester;
using Thyra::adjoint;
using Thyra::multiply;
using Thyra::diagonal;


std::string matrixFile = "";


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "matrix-file", &matrixFile,
    "Defines the Epetra_CrsMatrix to read in."  );
}

// helper function to excercise all different versions of B*D*G
const Teuchos::RCP<const Thyra::LinearOpBase<double> >
buildADOperator(int scenario,const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A,
                              const double vecScale, std::ostream & out)
{
   // Create the implicit operator
   double scalea=-7.0;
   double scaled=10.0;

   RCP<const LinearOpBase<double> > M;
   RCP<VectorBase<double> > d;
   if(scenario<=2) 
      d = createMember(A->domain(), "d");
   else 
      d = createMember(A->range(), "d");
   V_S( d.ptr(), vecScale ); // ToDo: Set ton != 1.0 and generalize
 
   // create an operator based on requested scenario 
   // (these are the same as in EpetraExt)
   switch(scenario) {
   case 1: 
      M = multiply( scale(scalea,A), scale(scaled,diagonal(d)), "M" );
      out << " Testing A*D" << std::endl;
      break;
   case 2: 
      M = multiply( scale(scaled,diagonal(d)), scale(scalea,A), "M" );
      out << " Testing D*A" << std::endl;
      break;
   case 3: 
      M = multiply( A, scale(scaled,diagonal(d)), "M" );
      out << " Testing adj(A)*D" << std::endl;
      break;
   case 4: 
      M = multiply( scale(scaled,diagonal(d)), A, "M" );
      out << " Testing D*adj(A)" << std::endl;
      break;
   default:
      TEUCHOS_ASSERT(false);
      break;
   }


   out << "\nM = " << *M;

   return M;
}

TEUCHOS_UNIT_TEST( EpetraExtDiagScalingTransformer, isCompatible)
{
   out << "\nReading linear system in Epetra format from the file \'"<<matrixFile<<"\'\n";

#ifdef HAVE_MPI
   Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
   Epetra_SerialComm comm;
#endif

   RCP<Epetra_CrsMatrix> epetra_A;
   RCP<Epetra_CrsMatrix> epetra_B;
   EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A, NULL, NULL, NULL );
   EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_B, NULL, NULL, NULL );
   const RCP<const Thyra::LinearOpBase<double> > A = Thyra::epetraLinearOp(epetra_A);
   const RCP<const Thyra::LinearOpBase<double> > B = Thyra::epetraLinearOp(epetra_B);

   RCP<VectorBase<double> > d = createMember(B->domain(), "d");
   V_S( d.ptr(), 3.0 ); // ToDo: Set ton != 1.0 and generalize

   const RCP<EpetraExtDiagScalingTransformer> BD_transformer =
         epetraExtDiagScalingTransformer();
 
   TEST_ASSERT(BD_transformer != null);

   RCP<const LinearOpBase<double> > M;

   M = multiply(A,B);
   TEST_ASSERT(not BD_transformer->isCompatible(*M));

   M = multiply(A,scale(3.9,diagonal(d)));
   TEST_ASSERT(BD_transformer->isCompatible(*M));

   M = multiply(scale(3.0,diagonal(d)),scale(9.0,B));
   TEST_ASSERT(BD_transformer->isCompatible(*M));

   M = multiply(B,scale(3.0,diagonal(d)),scale(9.0,B));
   TEST_ASSERT(not BD_transformer->isCompatible(*M));
}

TEUCHOS_UNIT_TEST( EpetraExtDiagScalingTransformer, reapply_scaling)
{
  //
  // A) Read in problem matrices
  //
  
  out << "\nReading linear system in Epetra format from the file \'"<<matrixFile<<"\'\n";
    
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  RCP<Epetra_CrsMatrix> epetra_A;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A, NULL, NULL, NULL );

  //
  // B) Create the Thyra wrapped version
  //
 
  const RCP<const Thyra::LinearOpBase<double> > A = Thyra::epetraLinearOp(epetra_A);

  RCP<const Thyra::LinearOpBase<double> > M;

  // test scale on right
  for(int scenario=1;scenario<=2;scenario++) {
     out << "**********************************************************" << std::endl;
     out << "**************** Starting Scenario = " << scenario << std::endl;
     out << "**********************************************************" << std::endl;

     //
     // D) Check compatibility
     //
     const RCP<EpetraExtDiagScalingTransformer> BD_transformer =
           epetraExtDiagScalingTransformer();

     TEST_ASSERT(BD_transformer != null);

     //
     // E) Do the transformation (twice)
     //

     const RCP<LinearOpBase<double> > M_explicit = BD_transformer->createOutputOp();
     TEST_ASSERT(M_explicit != null);

     // do trans
     M = buildADOperator(scenario,A,4.5,out);
     BD_transformer->transform( *M, M_explicit.ptr() );
     const RCP<const Epetra_Operator> eOp_explicit0 = Thyra::get_Epetra_Operator(*M_explicit);

     M = buildADOperator(scenario,A,7.5,out);
     BD_transformer->transform( *M, M_explicit.ptr() );
     const RCP<const Epetra_Operator> eOp_explicit1 = Thyra::get_Epetra_Operator(*M_explicit);

     // Epetra operator should not change!
     TEST_EQUALITY(eOp_explicit0,eOp_explicit1);

     out << "\nM_explicit = " << *M_explicit;

     //
     // F) Check the explicit operator
     //

     LinearOpTester<double> M_explicit_tester;
     M_explicit_tester.show_all_tests(true);;

     const bool result = M_explicit_tester.compare( *M, *M_explicit, &out );
     if (!result) success = false;
  }
}

TEUCHOS_UNIT_TEST( EpetraExtDiagScalingTransformer, basic_BDG_GDB)
{
  
  //
  // A) Read in problem matrices
  //
  
  out << "\nReading linear system in Epetra format from the file \'"<<matrixFile<<"\'\n";
    
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  RCP<Epetra_CrsMatrix> epetra_A;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A, NULL, NULL, NULL );
  
  //
  // B) Create the Thyra wrapped version
  //
 
  const RCP<const Thyra::LinearOpBase<double> > A = Thyra::epetraLinearOp(epetra_A);

  //
  // C) Create implicit B*D*B operator
  //

  // build scenario=1 -> B*D*B, scenario=2-> B*D*B', scenario=3 -> B'*D*B
  for(int scenario=1;scenario<=4;scenario++) {
     out << "**********************************************************" << std::endl;
     out << "**************** Starting Scenario = " << scenario << std::endl;
     out << "**********************************************************" << std::endl;

     RCP<const Thyra::LinearOpBase<double> > M = buildADOperator(scenario,A,4.5,out);

     //
     // D) Check compatibility
     //
     const RCP<EpetraExtDiagScalingTransformer> BD_transformer =
           epetraExtDiagScalingTransformer();

     TEST_ASSERT(BD_transformer != null);

     BD_transformer->isCompatible(*M);

     //
     // E) Do the transformation
     //

     const RCP<LinearOpBase<double> > M_explicit = BD_transformer->createOutputOp();
     BD_transformer->transform( *M, M_explicit.ptr() );

     out << "\nM_explicit = " << *M_explicit;

     //
     // F) Check the explicit operator
     //

     LinearOpTester<double> M_explicit_tester;
     M_explicit_tester.show_all_tests(true);;

     const bool result = M_explicit_tester.compare( *M, *M_explicit, &out );
     if (!result) success = false;
  }

}

} // namespace
