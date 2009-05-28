// Teuchos includes /*@ \label{lned:being-includes} @*/
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Epetra includes
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"

// EpetraExt includes
#include "EpetraExt_CrsMatrixIn.h"
#include "EpetraExt_VectorIn.h"

// PB-Package includes
#include "PB_Utilities.hpp"
#include "PB_InverseFactory.hpp"
#include "PB_InverseLibrary.hpp"
#include "Epetra/PB_StridedEpetraOperator.hpp"
#include "Epetra/PB_EpetraBlockPreconditioner.hpp"
#include "NS/PB_LSCPreconditionerFactory.hpp"
#include "NS/PB_SIMPLEPreconditionerFactory.hpp"

// Aztec includes
#include "AztecOO.h"
#include "AztecOO_Operator.h"

#include <iostream> /*@ \label{lned:end-includes} @*/

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

int main(int argc,char * argv[])
{
   // calls MPI_Init and MPI_Finalize
   Teuchos::GlobalMPISession mpiSession(&argc,&argv);

   // build global communicator
   Epetra_MpiComm Comm(MPI_COMM_WORLD);

   // Read in the matrix, store pointer as an RCP
   Epetra_CrsMatrix * ptrA = 0;
   EpetraExt::MatrixMarketFileToCrsMatrix("../data/nsjac_test.mm",Comm,ptrA);
   RCP<Epetra_CrsMatrix> A = rcp(ptrA);

   // allocate vectors
   RCP<Epetra_Vector> b = rcp(new Epetra_Vector(A->OperatorRangeMap()));
   RCP<Epetra_Vector> x = rcp(new Epetra_Vector(A->OperatorDomainMap()));

   b->Random();
   x->PutScalar(0.0);

   // Break apart the strided linear system
   /////////////////////////////////////////////////////////

   // Block the linear system using a strided epetra operator
   std::vector<int> vec(2); vec[0] = 2; vec[1] = 1; /*@ \label{lned:define-strided} @*/
   PB::Epetra::StridedEpetraOperator sA(vec,A);

   // Build the preconditioner /*@ \label{lned:construct-prec} @*/
   /////////////////////////////////////////////////////////

   // build an InverseLibrary
   RCP<PB::InverseLibrary> invLib = PB::InverseLibrary::buildFromStratimikos();
   
   // build the inverse factory needed by the example preconditioner
   RCP<const PB::InverseFactory> inverse  /*@ \label{lned:define-inv-fact} @*/
         = invLib->getInverseFactory("Amesos");

   // build the preconditioner factory
   // RCP<PB::NS::LSCStrategy> strategy = rcp(new PB::NS::InvLSCStrategy(inverse,true)); /*@ \label{lned:const-prec-strategy} @*/
   // RCP<PB::BlockPreconditionerFactory> precFact /*@ \label{lned:const-prec-fact} @*/
   //       = rcp(new PB::NS::LSCPreconditionerFactory(strategy));
   RCP<PB::BlockPreconditionerFactory> precFact /*@ \label{lned:const-prec-fact} @*/
          = rcp(new PB::NS::SIMPLEPreconditionerFactory(inverse,0.9));

   // using the preconditioner factory construct an Epetra_Operator
   PB::Epetra::EpetraBlockPreconditioner prec(precFact); /*@ \label{lned:const-epetra-prec} @*/
   prec.buildPreconditioner(sA); // fill epetra preconditioner using the strided operator

   // Build and solve the linear system
   /////////////////////////////////////////////////////////

   // Setup the linear solve: notice A is used directly 
   Epetra_LinearProblem problem(&*A,&*x,&*b); /*@ \label{lned:aztec-solve} @*/

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
