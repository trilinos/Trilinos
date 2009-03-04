
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


TEUCHOS_STATIC_SETUP()
{
  Teuchos::UnitTestRepository::getCLP().setOption(
    "matrix-file", &matrixFile,
    "Defines the Epetra_CrsMatrix to read in."  );
}



TEUCHOS_UNIT_TEST( EpetraExtDiagScaledMatProdTransformer, basic )
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

  RCP<const Thyra::LinearOpBase<double> > B = Thyra::epetraLinearOp( epetra_B, "B" );

  out << "\nB = " << *B;

  //
  // C) Create the implicit operator
  //

  const RCP<VectorBase<double> > d = createMember(B->range(), "d");
  V_S( d.ptr(), 1.0 ); // ToDo: Set ton != 1.0 and generalize
  const RCP<const LinearOpBase<double> > M =
    //multiply( adjoint(B), diagonal(d), B, "M" );
    multiply( B, diagonal(d), B, "M" );

  out << "\nM = " << *M;

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


} // namespace


