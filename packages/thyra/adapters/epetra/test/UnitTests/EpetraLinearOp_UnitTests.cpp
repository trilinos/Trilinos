#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Epetra_SerialComm.h"
#ifdef HAVE_MPI
#  include "Epetra_MpiComm.h"
#endif
#include "Epetra_Map.h"
#include "Epetra_CrsMatrix.h"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_as.hpp"

#include "Teuchos_UnitTestHarness.hpp"

#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_MultiVectorBase.hpp"


namespace {


//
// Helper code and declarations
//

using Teuchos::as;
using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::Array;


const int g_localDim = 4; // ToDo: Make variable!


RCP<const Epetra_Comm> getEpetraComm()
{
#ifdef HAVE_MPI
  return rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
#else
  return rcp(new Epetra_SerialComm());
#endif
}


const int numLocalRows = g_localDim;

RCP<Epetra_CrsMatrix> getEpetraMatrix(int numRows,int numCols,double shift=0.0) 
{
  const RCP<const Epetra_Comm> comm = getEpetraComm();

  const Epetra_Map rowMap(numRows, 0, *comm);
  const Epetra_Map domainMap(numCols, numCols, 0, *comm);
 
  const RCP<Epetra_CrsMatrix> epetraCrsM =
    rcp(new Epetra_CrsMatrix(Copy, rowMap, numCols));

  Array<double> rowEntries(numCols);
  Array<int> columnIndices(numCols);
  for (int j = 0; j < numCols; ++j) {
    columnIndices[j] = j;
  }

  for (int i = 0; i < numLocalRows; ++i) {
    
    for (int j = 0; j < numCols; ++j) {
      rowEntries[j] = as<double>(i+1) + as<double>(j+1) / 10 + shift;
    }

    epetraCrsM->InsertMyValues( i, numCols, &rowEntries[0], &columnIndices[0] );

  }

  epetraCrsM->FillComplete(domainMap, rowMap);

  return epetraCrsM;
}

//
// Unit Tests
//


TEUCHOS_UNIT_TEST( EpetraLinearOp, rectangular )
{
  const RCP<const Epetra_Comm> comm = getEpetraComm();

  const int numLocalRows = g_localDim;
  const int numRows = numLocalRows * comm->NumProc();
  const int numCols = numLocalRows / 2;

  RCP<Epetra_CrsMatrix> epetraCrsM = getEpetraMatrix(numRows,numCols);
/*
  const Epetra_Map rowMap(-1, numLocalRows, 0, *comm);

  const int numCols = rowMap.NumGlobalElements() / 2;

  const Epetra_Map domainMap(numCols, numCols, 0, *comm);
 
  const RCP<Epetra_CrsMatrix> epetraCrsM =
    rcp(new Epetra_CrsMatrix(Copy, rowMap, numCols));

  Array<double> rowEntries(numCols);
  Array<int> columnIndices(numCols);
  for (int j = 0; j < numCols; ++j) {
    columnIndices[j] = j;
  }

  for (int i = 0; i < numLocalRows; ++i) {
    
    for (int j = 0; j < numCols; ++j) {
      rowEntries[j] = as<double>(i+1) + as<double>(j+1) / 10;
    }

    epetraCrsM->InsertMyValues( i, numCols, &rowEntries[0], &columnIndices[0] );

  }

  epetraCrsM->FillComplete(domainMap, rowMap);
*/

  const RCP<const Thyra::LinearOpBase<double> > epetraOp =
    Thyra::epetraLinearOp(epetraCrsM);

  Thyra::LinearOpTester<double> linearOpTester;

  const bool result = linearOpTester.check(*epetraOp, &out);

  if (!result)
    success = false;
   
}

/*
TEUCHOS_UNIT_TEST( EpetraLinearOp, blocked_op )
{
   // build sub operators
   RCP<const Thyra::LinearOpBase<double> > A00 = Thyra::epetraLinearOp(getEpetraMatrix(4,4,0));
   RCP<const Thyra::LinearOpBase<double> > A01 = Thyra::epetraLinearOp(getEpetraMatrix(4,3,1));
   RCP<const Thyra::LinearOpBase<double> > A02 = Thyra::epetraLinearOp(getEpetraMatrix(4,2,2));
   RCP<const Thyra::LinearOpBase<double> > A10 = Thyra::epetraLinearOp(getEpetraMatrix(3,4,3));
   RCP<const Thyra::LinearOpBase<double> > A11 = Thyra::epetraLinearOp(getEpetraMatrix(3,3,4));
   RCP<const Thyra::LinearOpBase<double> > A12 = Thyra::epetraLinearOp(getEpetraMatrix(3,2,5));
   RCP<const Thyra::LinearOpBase<double> > A20 = Thyra::epetraLinearOp(getEpetraMatrix(2,4,6));
   RCP<const Thyra::LinearOpBase<double> > A21 = Thyra::epetraLinearOp(getEpetraMatrix(2,3,8));
   RCP<const Thyra::LinearOpBase<double> > A22 = Thyra::epetraLinearOp(getEpetraMatrix(2,2,8));

   out << "Sub operators built" << std::endl;

   // build composite operator
   RCP<const Thyra::LinearOpBase<double> > A 
         = Thyra::block2x2<double>(Thyra::block2x2<double>(A00,A01,A10,A11), Thyra::block2x1<double>(A02,A12),
                                   Thyra::block1x2<double>(A20,A21),         A22);

   out << "Composite operator built" << std::endl;

   // build vectors for use in apply
   RCP<Thyra::MultiVectorBase<double> > x = Thyra::createMembers<double>(A->domain(),3);
   RCP<Thyra::MultiVectorBase<double> > y = Thyra::createMembers<double>(A->range(),3);
 
   Thyra::randomize(-1.0,1.0,x.ptr().get());

   out << "A = \n" << Teuchos::describe(*A,Teuchos::VERB_HIGH) << std::endl;
   out << "x = \n" << Teuchos::describe(*x,Teuchos::VERB_HIGH) << std::endl;
   out << "y = \n" << Teuchos::describe(*y,Teuchos::VERB_HIGH) << std::endl;

   try {
      // perform a matrix vector multiply
      A->apply(Thyra::NONCONJ_ELE,*x,&*y);
   }
   catch(std::exception & exp) {
      out << "\"apply\" caused an exception" << std::endl;
      out << exp.what() << std::endl; 
      success = false;
   }

   out << "Test complete" << std::endl;
}
*/


} // namespace
