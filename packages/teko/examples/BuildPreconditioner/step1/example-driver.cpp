#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_FancyOStream.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Thyra includes
#include "Thyra_LinearOpBase.hpp"
#include "Thyra_LinearOpWithSolveFactoryHelpers.hpp"
#include "Thyra_PreconditionerFactoryHelpers.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultScaledAdjointLinearOp.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultPreconditioner.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"

// include basic Epetra information
#ifdef HAVE_MPI
   #include "Epetra_MpiComm.h"
   #include "mpi.h"
#else
   #include "Epetra_SerialComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"

// EpetraExt 
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// PB-Package includes
#include "PB_Utilities.hpp"
#include "PB_InverseFactory.hpp"
#include "Epetra/PB_StridedEpetraOperator.hpp"
#include "Epetra/PB_EpetraBlockPreconditioner.hpp"

// Aztec includes
#include "AztecOO.h"
#include "AztecOO_Operator.h"

#include <iostream>

// use would include your header
#include "ExamplePreconditionerFactory.cpp"

using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc,char * argv[])
{
   // calls MPI_Init and MPI_Finalize
   Teuchos::GlobalMPISession mpiSession(&argc,&argv);

   // build global (or serial communicator
   #ifdef HAVE_MPI
      Epetra_MpiComm Comm(MPI_COMM_WORLD);
   #else
      Epetra_SerialComm Comm;
   #endif


   // EpetraExt::MatrixMarketFileToCrsMatrix("./jac_test.mm",map,map,map,ptrA);
   // Read in the matrix, store pointer as an RCP
   Epetra_CrsMatrix * ptrA = 0;
   EpetraExt::MatrixMarketFileToCrsMatrix("../../data/jac_test.mm",Comm,ptrA);
   RCP<Epetra_CrsMatrix> A = rcp(ptrA);

   // allocate vectors
   RCP<Epetra_Vector> b = rcp(new Epetra_Vector(A->OperatorRangeMap()));
   RCP<Epetra_Vector> x = rcp(new Epetra_Vector(A->OperatorDomainMap()));

   b->Random();
   x->PutScalar(0.0);

   // Block the linear system using a strided epetra operator
   std::vector<int> vec(2); vec[0] = 2; vec[1] = 1;
   PB::Epetra::StridedEpetraOperator sA(vec,A);

   // Set parameters for the inverse factory
   RCP<const Teuchos::ParameterList> paramList = PB::invFactoryValidParameters();
   
   // build the inverse factory needed by the example preconditioner
   RCP<const PB::InverseFactory> inverse = PB::invFactoryFromParamList(*paramList,"ML");

   // build the preconditioner factory
   double alpha = 0.9;
   RCP<PB::BlockPreconditionerFactory> precFact 
         = rcp(new ExamplePreconditionerFactory(inverse,alpha));

   // using the preconditioner factory construct an Epetra_Operator
   PB::Epetra::EpetraBlockPreconditioner prec(precFact);
   prec.buildPreconditioner(sA); // fill the Epetra_Operator based on the strided operator

   // Setup the linear solve: notice A is used directly
   Epetra_LinearProblem problem(&*A,&*x,&*b);

   // build the solver
   AztecOO::AztecOO solver(problem);
   solver.SetAztecOption(AZ_solver,AZ_gmres);
   solver.SetAztecOption(AZ_precond,AZ_none);
   solver.SetAztecOption(AZ_kspace,50);
   solver.SetAztecOption(AZ_output,10);
   solver.SetPrecOperator(&prec);

   // solve the linear system
   solver.Iterate(1000,1e-5);

   return 0;
}
