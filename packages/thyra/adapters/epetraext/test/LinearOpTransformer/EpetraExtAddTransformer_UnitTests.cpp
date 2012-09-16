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
// Questions? Contact Roscoe A. Bartlett (bartlettra@ornl.gov) 
// 
// ***********************************************************************
// @HEADER
*/


#include "Thyra_EpetraExtAddTransformer.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultAddedLinearOp.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_DefaultIdentityLinearOp.hpp"
#include "EpetraExt_readEpetraLinearSystem.h"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_VerboseObject.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"
#include "Teuchos_CommandLineProcessor.hpp"
#include "Teuchos_StandardCatchMacros.hpp"

#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#else
#  include "Epetra_SerialComm.h"
#endif

#include "Teuchos_UnitTestHarness.hpp"


namespace {


using Teuchos::null;
using Teuchos::RCP;
using Thyra::EpetraExtAddTransformer;
using Thyra::epetraExtAddTransformer;
using Thyra::VectorBase;
using Thyra::LinearOpBase;
using Thyra::createMember;
using Thyra::LinearOpTester;
using Thyra::adjoint;
using Thyra::multiply;
using Thyra::diagonal;


std::string matrixFile = "";
std::string matrixFile2 = "";


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "matrix-file", &matrixFile,
    "Defines the Epetra_CrsMatrix to read in."  );
  Teuchos::UnitTestRepository::getCLP().setOption(
    "matrix-file-2", &matrixFile2,
    "Defines the Epetra_CrsMatrix to read in."  );
}

const Teuchos::RCP<const Thyra::LinearOpBase<double> >
buildAddOperator(int scenario,const Teuchos::RCP<const Thyra::LinearOpBase<double> > & A,
                               const Teuchos::RCP<const Thyra::LinearOpBase<double> > & B)
{
   // build operators for the various addition/adjoint scenarios
   RCP<const Thyra::LinearOpBase<double> > M;
   
   switch(scenario) {
   case 0:
      M = Thyra::add(A,B,"A+B");
      break;
   case 1:
      M = Thyra::add(A,B,"A+adj(B)");
      break;
   case 2:
      M = Thyra::add(A,B,"adj(A)+B");
      break;
   case 3:
      M = Thyra::add(A,B,"adb(A)+adb(B)");
      break;
   default:
      TEUCHOS_ASSERT(false);
      break;
   }

   return M;
}

TEUCHOS_UNIT_TEST( EpetraExtAddTransformer, diag_mat_Add )
{
  
  //
  // A) Read in problem matrices
  //
  
  out << "\nReading linear system in Epetra format from the file \'"<<matrixFile<<"\' ...\n";
    
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  RCP<Epetra_CrsMatrix> crsMat;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &crsMat, NULL, NULL, NULL );
  
  //
  // B) Create the Thyra wrapped version
  //
  double scaleMat=3.7;
  //double scaleDiag=-2.9;

 
  const RCP<const Thyra::LinearOpBase<double> > A = Thyra::scale<double>(scaleMat,Thyra::epetraLinearOp(crsMat));
  RCP<VectorBase<double> > d = createMember(A->domain(), "d");
  V_S( d.ptr(), 3.0 ); // ToDo: Set ton != 1.0 and generalize
  const RCP<const Thyra::LinearOpBase<double> > B = diagonal(d);

  out << "\nA = " << *A;
  out << "\nB = " << *B;

  for(int scenario=0;scenario<4;scenario++) {
     //
     // C) Create implicit A+B operator
     //
   
     const RCP<const Thyra::LinearOpBase<double> > M = buildAddOperator(scenario,A,B);
   
     //
     // D) Do the transformation
     //
   
     const RCP<EpetraExtAddTransformer> ApB_transformer = epetraExtAddTransformer();
   
     TEST_ASSERT(ApB_transformer != null);
   
     const RCP<LinearOpBase<double> > M_explicit = ApB_transformer->createOutputOp();
     ApB_transformer->transform( *M, M_explicit.ptr() );
   
     out << "\nM_explicit = " << *M_explicit;
   
     //
     // E) Check the explicit operator
     //
   
     LinearOpTester<double> M_explicit_tester;
     M_explicit_tester.show_all_tests(true);;
   
     const bool result = M_explicit_tester.compare( *M, *M_explicit, Teuchos::inOutArg(out) );
     if (!result) success = false;
  }

  for(int scenario=0;scenario<4;scenario++) {
     //
     // C) Create implicit A+B operator
     //
   
     const RCP<const Thyra::LinearOpBase<double> > M = buildAddOperator(scenario,B,A);
   
     //
     // D) Do the transformation
     //
   
     const RCP<EpetraExtAddTransformer> ApB_transformer = epetraExtAddTransformer();
   
     TEST_ASSERT(ApB_transformer != null);
   
     const RCP<LinearOpBase<double> > M_explicit = ApB_transformer->createOutputOp();
     ApB_transformer->transform( *M, M_explicit.ptr() );
   
     out << "\nM_explicit = " << *M_explicit;
   
     //
     // E) Check the explicit operator
     //
   
     LinearOpTester<double> M_explicit_tester;
     M_explicit_tester.show_all_tests(true);;
   
     const bool result = M_explicit_tester.compare( *M, *M_explicit, Teuchos::inOutArg(out) );
     if (!result) success = false;
  }
}

TEUCHOS_UNIT_TEST( EpetraExtAddTransformer, id_mat_Add )
{
  
  //
  // A) Read in problem matrices
  //
  
  out << "\nReading linear system in Epetra format from the file \'"<<matrixFile<<"\' ...\n";
    
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  RCP<Epetra_CrsMatrix> crsMat;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &crsMat, NULL, NULL, NULL );
  
  //
  // B) Create the Thyra wrapped version
  //
  double scaleMat=3.7;
  //double scaleDiag=-2.9;

 
  const RCP<const Thyra::LinearOpBase<double> > A = Thyra::scale<double>(scaleMat,Thyra::epetraLinearOp(crsMat));
  const RCP<const Thyra::LinearOpBase<double> > B = identity(A->domain(),"id");

  out << "\nA = " << *A;
  out << "\nB = " << *B;

  for(int scenario=0;scenario<4;scenario++) {
     //
     // C) Create implicit A+B operator
     //
   
     const RCP<const Thyra::LinearOpBase<double> > M = buildAddOperator(scenario,A,B);
   
     //
     // D) Do the transformation
     //
   
     const RCP<EpetraExtAddTransformer> ApB_transformer = epetraExtAddTransformer();
   
     TEST_ASSERT(ApB_transformer != null);
   
     const RCP<LinearOpBase<double> > M_explicit = ApB_transformer->createOutputOp();
     ApB_transformer->transform( *M, M_explicit.ptr() );
   
     out << "\nM_explicit = " << *M_explicit;
   
     //
     // E) Check the explicit operator
     //
   
     LinearOpTester<double> M_explicit_tester;
     M_explicit_tester.show_all_tests(true);;
   
     const bool result = M_explicit_tester.compare( *M, *M_explicit, Teuchos::inOutArg(out) );
     if (!result) success = false;
  }

  for(int scenario=0;scenario<4;scenario++) {
     //
     // C) Create implicit A+B operator
     //
   
     const RCP<const Thyra::LinearOpBase<double> > M = buildAddOperator(scenario,B,A);
   
     //
     // D) Do the transformation
     //
   
     const RCP<EpetraExtAddTransformer> ApB_transformer = epetraExtAddTransformer();
   
     TEST_ASSERT(ApB_transformer != null);
   
     const RCP<LinearOpBase<double> > M_explicit = ApB_transformer->createOutputOp();
     ApB_transformer->transform( *M, M_explicit.ptr() );
   
     out << "\nM_explicit = " << *M_explicit;
   
     //
     // E) Check the explicit operator
     //
   
     LinearOpTester<double> M_explicit_tester;
     M_explicit_tester.show_all_tests(true);;
   
     const bool result = M_explicit_tester.compare( *M, *M_explicit, Teuchos::inOutArg(out) );
     if (!result) success = false;
  }
}

TEUCHOS_UNIT_TEST( EpetraExtAddTransformer, basic_Add )
{
  
  //
  // A) Read in problem matrices
  //
  
  out << "\nReading linear system in Epetra format from the file \'"<<matrixFile<<"\' ...\n";
  out << "\nReading linear system in Epetra format from the file \'"<<matrixFile2<<"\' ...\n";
    
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  RCP<Epetra_CrsMatrix> epetra_A;
  RCP<Epetra_CrsMatrix> epetra_B;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A, NULL, NULL, NULL );
  EpetraExt::readEpetraLinearSystem( matrixFile2, comm, &epetra_B, NULL, NULL, NULL );
  
  //
  // B) Create the Thyra wrapped version
  //
  double scaleA=3.7;
  double scaleB=-2.9;
 
  const RCP<const Thyra::LinearOpBase<double> > A = Thyra::scale<double>(scaleA,Thyra::epetraLinearOp(epetra_B));
  const RCP<const Thyra::LinearOpBase<double> > B = Thyra::scale<double>(scaleB,Thyra::epetraLinearOp(epetra_B));

  out << "\nA = " << *A;
  out << "\nB = " << *B;

  for(int scenario=0;scenario<4;scenario++) {
     //
     // C) Create implicit A+B operator
     //
   
     const RCP<const Thyra::LinearOpBase<double> > M = buildAddOperator(scenario,A,B);
   
     //
     // D) Do the transformation
     //
   
     const RCP<EpetraExtAddTransformer> ApB_transformer = epetraExtAddTransformer();
   
     TEST_ASSERT(ApB_transformer != null);
   
     const RCP<LinearOpBase<double> > M_explicit = ApB_transformer->createOutputOp();
     ApB_transformer->transform( *M, M_explicit.ptr() );
   
     out << "\nM_explicit = " << *M_explicit;
   
     //
     // E) Check the explicit operator
     //
   
     LinearOpTester<double> M_explicit_tester;
     M_explicit_tester.show_all_tests(true);;
   
     const bool result = M_explicit_tester.compare( *M, *M_explicit, Teuchos::inOutArg(out) );
     if (!result) success = false;
  }
}

TEUCHOS_UNIT_TEST( EpetraExtAddTransformer, mod_Add )
{
  
  //
  // A) Read in problem matrices
  //
  
  out << "\nReading linear system in Epetra format from the file \'"<<matrixFile<<"\' ...\n";
  out << "\nReading linear system in Epetra format from the file \'"<<matrixFile2<<"\' ...\n";
    
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  RCP<Epetra_CrsMatrix> epetra_A;
  RCP<Epetra_CrsMatrix> epetra_B;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_A, NULL, NULL, NULL );
  EpetraExt::readEpetraLinearSystem( matrixFile2, comm, &epetra_B, NULL, NULL, NULL );
  
  //
  // B) Create the Thyra wrapped version
  //
  double scaleA=3.7;
  double scaleB=-2.9;
 
  const RCP<const Thyra::LinearOpBase<double> > A = Thyra::scale<double>(scaleA,Thyra::epetraLinearOp(epetra_B));
  const RCP<const Thyra::LinearOpBase<double> > B = Thyra::scale<double>(scaleB,Thyra::epetraLinearOp(epetra_B));

  out << "\nA = " << *A;
  out << "\nB = " << *B;

  for(int scenario=0;scenario<4;scenario++) {
     //
     // C) Create implicit A+B operator
     //
   
     const RCP<const Thyra::LinearOpBase<double> > M = buildAddOperator(scenario,A,B);
   
     //
     // D) Do the transformation
     //
   
     const RCP<EpetraExtAddTransformer> ApB_transformer = epetraExtAddTransformer();
   
     TEST_ASSERT(ApB_transformer != null);
   
     const RCP<LinearOpBase<double> > M_explicit = ApB_transformer->createOutputOp();
     const RCP<LinearOpBase<double> > M_explicit_orig = M_explicit;
     ApB_transformer->transform( *M, M_explicit.ptr() );

     // do some violence to one of the operators
     double * view; int numEntries;
     epetra_B->Scale(3.2);
     TEUCHOS_ASSERT(epetra_B->ExtractMyRowView(3,numEntries,view)==0);
     for(int i=0;i<numEntries;i++) view[i] += view[i]*double(i+1);

     // compute multiply again
     ApB_transformer->transform( *M, M_explicit.ptr() );

     out << "\nM_explicit = " << *M_explicit;

     TEUCHOS_ASSERT(Thyra::get_Epetra_Operator(*M_explicit)
                  ==Thyra::get_Epetra_Operator(*M_explicit_orig));
   
     //
     // E) Check the explicit operator
     //
   
     LinearOpTester<double> M_explicit_tester;
     M_explicit_tester.show_all_tests(true);;
   
     const bool result = M_explicit_tester.compare( *M, *M_explicit, Teuchos::inoutArg(out) );
     if (!result) success = false;
  }
}

} // end namespace
