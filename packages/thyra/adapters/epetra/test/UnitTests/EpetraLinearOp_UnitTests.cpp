#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_LinearOpTester.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_MultiVectorBase.hpp"
#include "EpetraThyraAdaptersTestHelpers.hpp"

#include "Teuchos_UnitTestHarness.hpp"


namespace {


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

  const RCP<const Thyra::LinearOpBase<double> > epetraOp =
    Thyra::epetraLinearOp(epetraCrsM);

  Thyra::LinearOpTester<double> linearOpTester;

  const bool result = linearOpTester.check(*epetraOp, &out);

  if (!result)
    success = false;
   
}


TEUCHOS_UNIT_TEST( EpetraLinearOp, blocked_op )
{

  using Teuchos::describe;
  using Thyra::block2x2;
  using Thyra::block2x1;
  using Thyra::block1x2;
  using Thyra::LinearOpBase;
  using Thyra::MultiVectorBase;
  using Thyra::createMembers;
  
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

  {
     // build composite operator
     RCP<const LinearOpBase<double> > A = block2x2<double>(
       block2x2<double>(A00, A01, A10, A11),   block2x1<double>(A02,A12),
       block1x2<double>(A20, A21),             A22
       );
   
     out << "First composite operator built" << std::endl;
     
     // build vectors for use in apply
     RCP<MultiVectorBase<double> > x = createMembers<double>(A->domain(), 3);
     RCP<MultiVectorBase<double> > y = createMembers<double>(A->range(), 3);
     
     Thyra::randomize(-1.0,1.0,x.ptr().get());
   
     out << "A = \n" << describe(*A, Teuchos::VERB_HIGH) << std::endl;
     out << "x = \n" << describe(*x, Teuchos::VERB_HIGH) << std::endl;
     out << "y = \n" << describe(*y, Teuchos::VERB_HIGH) << std::endl;
     
     // perform a matrix vector multiply
     A->apply(Thyra::NONCONJ_ELE,*x,&*y);

     out << "First composite operator completed" << std::endl;
  }
/*
  {
     RCP<const LinearOpBase<double> > A = block2x2<double>(
       A11, block1x2<double>(A10,A12),block2x1<double>(A01,A21),
       block2x2<double>(A00,A02,A20,A22));
     
     out << "Second composite operator built" << std::endl;
     
     // build vectors for use in apply
     RCP<MultiVectorBase<double> > x = createMembers<double>(A->domain(), 3);
     RCP<MultiVectorBase<double> > y = createMembers<double>(A->range(), 3);
     
     Thyra::randomize(-1.0,1.0,x.ptr().get());
   
     out << "A = \n" << describe(*A, Teuchos::VERB_HIGH) << std::endl;
     out << "x = \n" << describe(*x, Teuchos::VERB_HIGH) << std::endl;
     out << "y = \n" << describe(*y, Teuchos::VERB_HIGH) << std::endl;
     
     // perform a matrix vector multiply
     A->apply(Thyra::NONCONJ_ELE,*x,&*y);

     out << "Second composite operator completed" << std::endl;
  }
*/

  out << "Test complete" << std::endl;

}


} // namespace
