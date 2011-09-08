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


#include "Thyra_EpetraExtDiagScaledMatProdTransformer.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_TestingTools.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_EpetraLinearOp.hpp"
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
using Thyra::EpetraExtDiagScaledMatProdTransformer;
using Thyra::epetraExtDiagScaledMatProdTransformer;
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
    "Defines the second Epetra_CrsMatrix to read in."  );
}

// helper function to excercise all different versions of B*D*G
const Teuchos::RCP<const Thyra::LinearOpBase<double> >
buildBDGOperator(int scenario,const Teuchos::RCP<const Thyra::LinearOpBase<double> > & B,
                              const Teuchos::RCP<const Thyra::LinearOpBase<double> > & G,
                              const double vecScale, std::ostream & out)
{
   // Create the implicit operator
   double scalea=10.0;
   double scaleb=-7.0;
   double scaled=52.0;

   RCP<const LinearOpBase<double> > M;
   RCP<VectorBase<double> > d;
   if(scenario<=2) 
      d = createMember(B->domain(), "d");
   else 
      d = createMember(B->range(), "d");
   V_S( d.ptr(), vecScale ); // ToDo: Set ton != 1.0 and generalize
 
   // create an operator based on requested scenario 
   // (these are the same as in EpetraExt)
   switch(scenario) {
   case 1: 
      M = multiply( scale(scalea,B), diagonal(d), scale(scaleb,G), "M" );
      out << " Testing B*D*G" << std::endl;
      break;
   case 2: 
      M = multiply( scale(scalea,B), diagonal(d), adjoint(G), "M" );
      out << " Testing B*D*adj(G)" << std::endl;
      break;
   case 3: 
      M = multiply( adjoint(B), diagonal(d), scale(scaleb,G), "M" );
      out << " Testing adj(B)*D*G" << std::endl;
      break;
   case 4: 
      M = multiply( adjoint(B), diagonal(d), adjoint(G), "M" );
      out << " Testing adj(B)*D*adj(G)" << std::endl;
      break;
   case 5: 
      M = multiply( B, scale(scaled,diagonal(d)), G, "M" );
      out << " Testing B*(52.0*D)*G" << std::endl;
      break;
   default:
      TEUCHOS_ASSERT(false);
      break;
   }


   out << "\nM = " << *M;

   return M;
}

// helper function to excercise all different versions of B*D*G
const Teuchos::RCP<const Thyra::LinearOpBase<double> >
buildBGOperator(int scenario,const Teuchos::RCP<const Thyra::LinearOpBase<double> > & B,
                             const Teuchos::RCP<const Thyra::LinearOpBase<double> > & G,
                             std::ostream & out)
{
   // Create the implicit operator
   double scalea=10.0;
   double scaleb=-7.0;
   RCP<const LinearOpBase<double> > M;

   // create an operator based on requested scenario 
   // (these are the same as in EpetraExt)
   switch(scenario) {
   case 1: 
      M = multiply( scale(scalea,B), scale(scaleb,G), "M" );
      out << " Testing B*G" << std::endl;
      break;
   case 2: 
      M = multiply( scale(scalea,B), adjoint(G), "M" );
      out << " Testing B*adj(G)" << std::endl;
      break;
   case 3: 
      M = multiply( adjoint(B), scale(scaleb,G), "M" );
      out << " Testing adj(B)*G" << std::endl;
      break;
   case 4: 
      M = multiply( adjoint(B), adjoint(G), "M" );
      out << " Testing adj(B)*adj(G)" << std::endl;
      break;
   default:
      TEUCHOS_ASSERT(false);
      break;
   }


   out << "\nM = " << *M;

   return M;
}

const Teuchos::RCP<const Thyra::LinearOpBase<double> >
buildBDBOperator(int scenario,const Teuchos::RCP<const Thyra::LinearOpBase<double> > & B,
                 const double vecScale, std::ostream & out)
{
   return buildBDGOperator(scenario,B,B,vecScale,out);
}



TEUCHOS_UNIT_TEST( EpetraExtDiagScaledMatProdTransformer, basic_BDB )
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
  RCP<Epetra_CrsMatrix> epetra_B;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_B, NULL, NULL, NULL );
  
  //
  // B) Create the Thyra wrapped version
  //
 
  const RCP<const Thyra::LinearOpBase<double> > B = Thyra::epetraLinearOp(epetra_B);

  //
  // C) Create implicit B*D*B operator
  //

  // build scenario=1 -> B*D*B, scenario=2 -> B*D*B',
  //       scenario=3 -> B'*D*B, scenario=4 -> B'*D*B'
  //int scenario = 3;
  for(int scenario=1;scenario<=5;scenario++) {
     const RCP<const Thyra::LinearOpBase<double> > M = buildBDBOperator(scenario,B,4.5,out);

     //
     // D) Do the transformation
     //

     const RCP<EpetraExtDiagScaledMatProdTransformer> BtDB_transformer =
       epetraExtDiagScaledMatProdTransformer();

     TEST_ASSERT(BtDB_transformer != null);

     const RCP<LinearOpBase<double> > M_explicit = BtDB_transformer->createOutputOp();

     BtDB_transformer->transform( *M, M_explicit.ptr() );

     out << "\nM_explicit = " << *M_explicit;

     //
     // E) Check the explicit operator
     //

     LinearOpTester<double> M_explicit_tester;
     M_explicit_tester.show_all_tests(true);;

     const bool result = M_explicit_tester.compare( *M, *M_explicit, &out );
     if (!result) success = false;
  }

}

TEUCHOS_UNIT_TEST( EpetraExtDiagScaledMatProdTransformer, basic_BDG_GDB)
{
  
  //
  // A) Read in problem matrices
  //
  
  out << "\nReading linear system in Epetra format from the file \'"<<matrixFile<<"\'";
  out << " and from the file \'"<<matrixFile2<<"\' ...\n";
    
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  RCP<Epetra_CrsMatrix> epetra_B;
  RCP<Epetra_CrsMatrix> epetra_G;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_B, NULL, NULL, NULL );
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_G, NULL, NULL, NULL );
  
  //
  // B) Create the Thyra wrapped version
  //
 
  const RCP<const Thyra::LinearOpBase<double> > B = Thyra::epetraLinearOp(epetra_B);
  const RCP<const Thyra::LinearOpBase<double> > G = Thyra::epetraLinearOp(epetra_G);

  //
  // C) Create implicit B*D*B operator
  //

  // build scenario=1 -> B*D*B, scenario=2-> B*D*B', scenario=3 -> B'*D*B
  for(int scenario=1;scenario<=3;scenario++) {
     RCP<const Thyra::LinearOpBase<double> > M;
     if(scenario==1 || scenario==3)
        M = buildBDGOperator(scenario,B,G,4.5,out);
     else
        M = buildBDGOperator(scenario,G,B,4.5,out);

     //
     // D) Do the transformation
     //

     const RCP<EpetraExtDiagScaledMatProdTransformer> BtDB_transformer =
       epetraExtDiagScaledMatProdTransformer();

     TEST_ASSERT(BtDB_transformer != null);

     const RCP<LinearOpBase<double> > M_explicit = BtDB_transformer->createOutputOp();

     BtDB_transformer->transform( *M, M_explicit.ptr() );

     out << "\nM_explicit = " << *M_explicit;

     //
     // E) Check the explicit operator
     //

     LinearOpTester<double> M_explicit_tester;
     M_explicit_tester.show_all_tests(true);;

     const bool result = M_explicit_tester.compare( *M, *M_explicit, &out );
     if (!result) success = false;
  }

}

TEUCHOS_UNIT_TEST( EpetraExtDiagScaledMatProdTransformer, basic_GDG )
{
  
  //
  // A) Read in problem matrices
  //
  
  out << "\nReading linear system in Epetra format from the file \'"<<matrixFile2<<"\' ...\n";
    
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  RCP<Epetra_CrsMatrix> epetra_G;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_G, NULL, NULL, NULL );
  
  //
  // B) Create the Thyra wrapped version
  //
 
  const RCP<const Thyra::LinearOpBase<double> > G = Thyra::epetraLinearOp(epetra_G);

  //
  // C) Create implicit B*D*B operator
  //

  // build scenario=1 -> B*D*B, scenario=3 -> B'*D*B
  int scenes[] = {1,4};
  for(int i=0;i<2;i++) {
     int scenario = scenes[i];
     const RCP<const Thyra::LinearOpBase<double> > M = buildBDGOperator(scenario,G,G,4.5,out);

     //
     // D) Do the transformation
     //

     const RCP<EpetraExtDiagScaledMatProdTransformer> BtDB_transformer =
       epetraExtDiagScaledMatProdTransformer();

     TEST_ASSERT(BtDB_transformer != null);

     const RCP<LinearOpBase<double> > M_explicit = BtDB_transformer->createOutputOp();

     BtDB_transformer->transform( *M, M_explicit.ptr() );

     out << "\nM_explicit = " << *M_explicit;

     //
     // E) Check the explicit operator
     //

     LinearOpTester<double> M_explicit_tester;
     M_explicit_tester.show_all_tests(true);;

     const bool result = M_explicit_tester.compare( *M, *M_explicit, &out );
     if (!result) success = false;
  }

}

TEUCHOS_UNIT_TEST( EpetraExtDiagScaledMatProdTransformer, basic_BG_GB_GG)
{
  
  //
  // A) Read in problem matrices
  //
  
  out << "\nReading linear system in Epetra format from the file \'"<<matrixFile<<"\'";
  out << " and from the file \'"<<matrixFile2<<"\' ...\n";
    
#ifdef HAVE_MPI
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  RCP<Epetra_CrsMatrix> epetra_B;
  RCP<Epetra_CrsMatrix> epetra_G;
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_B, NULL, NULL, NULL );
  EpetraExt::readEpetraLinearSystem( matrixFile, comm, &epetra_G, NULL, NULL, NULL );
  
  //
  // B) Create the Thyra wrapped version
  //
 
  const RCP<const Thyra::LinearOpBase<double> > B = Thyra::epetraLinearOp(epetra_B);
  const RCP<const Thyra::LinearOpBase<double> > G = Thyra::epetraLinearOp(epetra_G);

  //
  // C) Create implicit B*D*B operator
  //

  // build scenario=1 -> B*D*B, scenario=2-> B*D*B', scenario=3 -> B'*D*B
  for(int scenario=1;scenario<=4;scenario++) {
     RCP<const Thyra::LinearOpBase<double> > M;
     if(scenario==1 || scenario==3)
        M = buildBGOperator(scenario,B,G,out);
     else if(scenario==2)
        M = buildBGOperator(scenario,G,B,out);
     else
        M = buildBGOperator(scenario,G,G,out);

     //
     // D) Do the transformation
     //

     const RCP<EpetraExtDiagScaledMatProdTransformer> BtDB_transformer =
       epetraExtDiagScaledMatProdTransformer();

     TEST_ASSERT(BtDB_transformer != null);

     const RCP<LinearOpBase<double> > M_explicit = BtDB_transformer->createOutputOp();

     BtDB_transformer->transform( *M, M_explicit.ptr() );

     out << "\nM_explicit = " << *M_explicit;

     //
     // E) Check the explicit operator
     //

     LinearOpTester<double> M_explicit_tester;
     M_explicit_tester.show_all_tests(true);;

     const bool result = M_explicit_tester.compare( *M, *M_explicit, &out );
     if (!result) success = false;
  }

}

} // namespace
