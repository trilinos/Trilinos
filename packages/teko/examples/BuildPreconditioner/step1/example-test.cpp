// Teuchos includes /*@ \label{lnet:being-includes} @*/
#include "Teuchos_ConfigDefs.hpp"
#include "Teuchos_MPISession.hpp"
#include "Teuchos_GlobalMPISession.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_XMLParameterListHelpers.hpp"

// Thyra includes
#include "Thyra_EpetraLinearOp.hpp"

// Epetra includes
#include "mpi.h"
#include "Epetra_MpiComm.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_Vector.h"
#include "Epetra_LinearProblem.h"

// PB-Package includes
#include "PB_Utilities.hpp"
#include "PB_InverseFactory.hpp"
#include "PB_InverseLibrary.hpp"
#include "Epetra/PB_EpetraOperatorWrapper.hpp"
#include "Epetra/PB_EpetraBlockPreconditioner.hpp"

#include "ExamplePreconditionerFactory.cpp"

#include <iostream> /*@ \label{lnet:end-includes} @*/

// for simplicity
using Teuchos::RCP;
using Teuchos::rcp;

// utility function to construct Epetra operators
RCP<Epetra_CrsMatrix> build2x2(double a,double b,double c,double d,Epetra_Comm & Comm)
{
   RCP<Epetra_Map> map = rcp(new Epetra_Map(2,0,Comm));
   RCP<Epetra_CrsMatrix> matrix = rcp(new Epetra_CrsMatrix(Copy,*map,2));

   int indicies[2] = { 0, 1 };
   double values[2];

   // build first row
   values[0] = a;
   values[1] = b;
   matrix->InsertGlobalValues(0,2,values,indicies);

   // build second row
   values[0] = c;
   values[1] = d;
   matrix->InsertGlobalValues(1,2,values,indicies);

   // pack it up
   matrix->FillComplete();

   return matrix;
}

int main(int argc,char * argv[])
{
   // calls MPI_Init and MPI_Finalize
   Teuchos::GlobalMPISession mpiSession(&argc,&argv);

   // build global communicator
   Epetra_MpiComm Comm(MPI_COMM_WORLD);

   // build the sub blocks
   PB::LinearOp A_00 = Thyra::epetraLinearOp(build2x2(1,2,2,1,Comm));
   PB::LinearOp A_01 = Thyra::epetraLinearOp(build2x2(0,-1,3,4,Comm));
   PB::LinearOp A_10 = Thyra::epetraLinearOp(build2x2(1,6,-3,2,Comm));
   PB::LinearOp A_11 = Thyra::epetraLinearOp(build2x2(2,1,1,2,Comm));

   // build the Epetra operator
   PB::Epetra::EpetraOperatorWrapper A(PB::block2x2(A_00,A_01,A_10,A_11));

   // build the Epetra vector
   Epetra_Vector b(A.OperatorRangeMap());
   Epetra_Vector x(A.OperatorDomainMap());
   x.PutScalar(0.0);

   // build the RHS vector
   b[0] = 1;
   b[1] = 2;
   b[2] = 3;
   b[3] = 4;

   // Build the preconditioner /*@ \label{lnet:construct-prec} @*/
   /////////////////////////////////////////////////////////

   // build an InverseLibrary
   RCP<PB::InverseLibrary> invLib = PB::InverseLibrary::buildFromStratimikos(); /*@ \label{lnet:define-inv-params} @*/
   
   // build the inverse factory needed by the example preconditioner
   RCP<const PB::InverseFactory> inverse  /*@ \label{lnet:define-inv-fact} @*/
         = invLib->getInverseFactory("Amesos");

   // build the preconditioner factory
   RCP<PB::BlockPreconditionerFactory> precFact /*@ \label{lnet:const-prec-fact} @*/
          = rcp(new ExamplePreconditionerFactory(inverse,0.9));

   // using the preconditioner factory construct an Epetra_Operator
   PB::Epetra::EpetraBlockPreconditioner prec(precFact); /*@ \label{lnet:const-epetra-prec} @*/
   prec.buildPreconditioner(A); // fill epetra preconditioner using the strided operator

   // apply the precondtioner
   prec.ApplyInverse(b,x);

   x.Print(std::cout);

   return 0;
}
