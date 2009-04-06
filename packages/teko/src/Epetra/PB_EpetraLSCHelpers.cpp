#include "PB_EpetraLSCHelpers.hpp"

// Thyra Includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"

// Epetra includes
#include "Epetra_Vector.h"

// EpetraExt includes
#include "EpetraExt_ProductOperator.h"
#include "EpetraExt_MatrixMatrix.h"

// PB includes
#include "PB_EpetraOperatorWrapper.hpp"
#include "PB_Helpers.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::null;

namespace PB {
namespace Epetra {

// This routine will construct the LSC matrices needed for both stabilized and
// stable discretizations. The difference will depend on if the A operator has
// a stabilization matrix or not. If it does not the out_aD parameter will not be
// set. For a reference on the LSC preconditioners see
//
//    Elman, Howle, Shadid, Silvester, and Tuminaro, "Least Squares Preconditioners
//    for Stabilized Discretizations of the Navier-Stokes Euqations," SISC-2007.
void buildLSCOperators(const Epetra_Operator & in_A,const Epetra_RowMatrix & in_Qu,
                       RCP<const Epetra_CrsMatrix> & out_BQBtmgC,RCP<const Epetra_Vector> & out_aiD)
{
   using namespace Thyra;

   // cast it to an EpetraOperatorWrapper: if it fails complain!
   const RCP<const PB::Epetra::EpetraOperatorWrapper> wrapA 
      = rcp_dynamic_cast<const PB::Epetra::EpetraOperatorWrapper>(rcp(&in_A,false),true);

   std::vector<Teuchos::RCP<const Epetra_CrsMatrix> > blocks;
   std::pair<int,int> shape = thyraMatrixToCrsVector(wrapA->getThyraOp(),true,blocks);

   // this better be a 2x2 matrix!
   TEST_FOR_EXCEPTION(shape!=std::make_pair(2,2),std::runtime_error,"\"buildLSCOperators\" requires a 2x2 block matrix");

   // cast everything to a Epetra_CrsMatrix
   buildLSCOperators(blocks[0],blocks[1],blocks[2],blocks[3],rcpFromRef(in_Qu),out_BQBtmgC,out_aiD); 
}

void buildLSCOperators(const Epetra_Operator & in_A,Teuchos::RCP<const Epetra_CrsMatrix> & out_BQBtmgC,
                       Teuchos::RCP<const Epetra_Vector> & out_aiD)
{
   // cast it to an EpetraOperatorWrapper: if it fails complain!
   const RCP<const PB::Epetra::EpetraOperatorWrapper> wrapA 
      = rcp_dynamic_cast<const PB::Epetra::EpetraOperatorWrapper>(rcp(&in_A,false),true);

   std::vector<Teuchos::RCP<const Epetra_CrsMatrix> > blocks;
   std::pair<int,int> shape = thyraMatrixToCrsVector(wrapA->getThyraOp(),true,blocks);

   // this better be a 2x2 matrix!
   TEST_FOR_EXCEPTION(shape!=std::make_pair(2,2),std::runtime_error,"\"buildLSCOperators\" requires a 2x2 block matrix");

   // cast everything to a Epetra_CrsMatrix
   buildLSCOperators(blocks[0],blocks[1],blocks[2],blocks[3],Teuchos::null,out_BQBtmgC,out_aiD); 
}

// function to be used internally, handles both stable and stabilized discretizations
void buildLSCOperators(const Teuchos::RCP<const Epetra_CrsMatrix> & F,const Teuchos::RCP<const Epetra_CrsMatrix> & Bt,
                       const Teuchos::RCP<const Epetra_CrsMatrix> & B,const Teuchos::RCP<const Epetra_CrsMatrix> & C,
                       const Teuchos::RCP<const Epetra_RowMatrix> & in_Qu,
                       RCP<const Epetra_CrsMatrix> & out_BQBtmgC,RCP<const Epetra_Vector> & out_aiD)
{
   using namespace Thyra;

   // setup some human readable boolean flags
   bool haveMass = not is_null(in_Qu);
   bool stableDiscretization = is_null(C);
   #ifdef STUPID_EPETRA_ERRS_OFF
   int oldTrace = 0;
   #endif

   // compute what is needed for a stable (and stabilized!) discretization
   ////////////////////////////////////////////////////////////////////////////////

   RCP<Epetra_Vector> idQu; // inv(diag(Qu)) - as necessary

   if(haveMass) {
      // use the mass
      idQu = rcp(new Epetra_Vector(in_Qu->OperatorRangeMap()));
      TEST_FOR_EXCEPT(in_Qu->ExtractDiagonalCopy(*idQu));
      TEST_FOR_EXCEPT(idQu->Reciprocal(*idQu));
   }

   // now compute B*idQu*Bt
   // if this is a stable discretization we are done!
   const RCP<Epetra_CrsMatrix> BBt = rcp(new Epetra_CrsMatrix(Copy,B->RowMap(),0));
   RCP<const Epetra_CrsMatrix> BQ;

   if(haveMass) {
      // use the mass

      RCP<Epetra_CrsMatrix> tempBQ = rcp(new Epetra_CrsMatrix(*B)); 
         // this extra variable is needed becuase the argument B is a "const Epetra_CrsMatrix"

      TEST_FOR_EXCEPT(tempBQ->RightScale(*idQu));  // B = B * inv(diag(Q_u))

      BQ = tempBQ;
   }
   else {
      // don't use a mass
      BQ = B;
   }

   // compute BQBt = B*inv(diag(Q))*Bt
   
   #ifdef STUPID_EPETRA_ERRS_OFF
   oldTrace = BBt->GetTracebackMode(); 
   BBt->SetTracebackMode(0); 
   #endif
   //TEST_FOR_EXCEPT
   (EpetraExt::MatrixMatrix::Multiply(*BQ,false,*Bt,false,*BBt));
   #ifdef STUPID_EPETRA_ERRS_OFF
   BBt->SetTracebackMode(oldTrace); 
   #endif
 
   if(stableDiscretization) {
      // stable discretization...we are done!
      out_BQBtmgC = BBt;
      out_aiD = Teuchos::null;

      return;
   }

   // not a stable discretization: compute what is needed for stabilized discretization
   ////////////////////////////////////////////////////////////////////////////////

   // compute inv(diag(F))
   RCP<Epetra_Vector> idF = rcp(new Epetra_Vector(F->OperatorRangeMap()));
   TEST_FOR_EXCEPT(F->ExtractDiagonalCopy(*idF));
   TEST_FOR_EXCEPT(idF->Reciprocal(*idF));

   // compute diag(C)
   Epetra_Vector dC(C->OperatorRangeMap());
   TEST_FOR_EXCEPT(C->ExtractDiagonalCopy(dC));

   // construct D = diag(B*idF*Bt-C) = diag(B*idF*Bt)-diag(C)
   RCP<Epetra_CrsMatrix> BidFBt = rcp(new Epetra_CrsMatrix(Copy,B->OperatorRangeMap(),Bt->OperatorDomainMap(),0));
   Epetra_CrsMatrix BidF(*B);
   TEST_FOR_EXCEPT(BidF.RightScale(*idF));  // B = B * inv(diag(F))
   #ifdef STUPID_EPETRA_ERRS_OFF
   oldTrace = BidFBt->GetTracebackMode(); 
   BidFBt->SetTracebackMode(0); 
   #endif
   //TEST_FOR_EXCEPT(EpetraExt::MatrixMatrix::Multiply(BidF,false,*Bt,false,*BidFBt));
   (EpetraExt::MatrixMatrix::Multiply(BidF,false,*Bt,false,*BidFBt));
   #ifdef STUPID_EPETRA_ERRS_OFF
   BidFBt->SetTracebackMode(oldTrace); 
   #endif

   // compute diag(B*idF*Bt) and compute diag(B*idF*Bt)-diag(C)
   RCP<Epetra_Vector> pDiag = rcp(new Epetra_Vector(C->OperatorRangeMap())); // this vector can be reused
   TEST_FOR_EXCEPT(BidFBt->ExtractDiagonalCopy(*pDiag)); // diag(B*idF*Bt)
   pDiag->Update(-1.0,dC,1.0);                           // D = diag(B*idF*Bt)-diag(C)
   pDiag->Reciprocal(*pDiag);                            // inv(D)

   // compute alpha and gamma parameters
   ////////////////////////////////////////////////////////////////////////////////

   // multiply diagonal matrix times F matrix
   RCP<const LinearOpBase<double> > iQupF;

   if(haveMass) {
      // we got mass!
      const RCP<VectorBase<double> > thyraIDQu = create_Vector(idQu,create_VectorSpace(rcp(&in_Qu->OperatorRangeMap(),false)));
      iQupF = multiply<double>(diagonal(thyraIDQu),epetraLinearOp(F));
   }
   else {
      // no mass
      iQupF = epetraLinearOp(F);
   }

   // do 6 power iterations to compute spectral radius: EHSST2007 Eq. 4.28
   double gamma = std::abs(PB::computeSpectralRad(iQupF,5e-2,false,5))/3.0; 

   // do 6 power iterations to compute spectral radius: EHSST2007 Eq. 4.29 (w/ correct sign!)
   const RCP<VectorBase<double> > thyraIDD
         = create_Vector(pDiag,create_VectorSpace(Teuchos::rcpFromRef(C->OperatorRangeMap())));
   const RCP<const LinearOpBase<double> > BidFBtidD = multiply<double>(epetraLinearOp(BidFBt),diagonal(thyraIDD));
   double alpha = 1.0/std::abs(PB::computeSpectralRad(BidFBtidD,5e-2,false,5));
   
   // scale: pDiag = alpha*pDiag
   pDiag->Scale(alpha);

   // build B*idQ*Bt - gamma*C
   ////////////////////////////////////////////////////////////////////////////////
   RCP<Epetra_CrsMatrix> BQBtmgC = rcp(new Epetra_CrsMatrix(Copy,B->OperatorRangeMap(),0));
   #ifdef STUPID_EPETRA_ERRS_OFF
   oldTrace = BQBtmgC->GetTracebackMode(); 
   BQBtmgC->SetTracebackMode(0); 
   #endif
   TEST_FOR_EXCEPT(EpetraExt::MatrixMatrix::Add(*C,false,-gamma,*BQBtmgC,0.0)); // BQBtmgC = -gamma*C
   TEST_FOR_EXCEPT(EpetraExt::MatrixMatrix::Add(*BBt,false,1.0,*BQBtmgC,1.0)); // BQBtmgC += B*Q*Bt
   BQBtmgC->FillComplete();
   #ifdef STUPID_EPETRA_ERRS_OFF
   BQBtmgC->SetTracebackMode(oldTrace); 
   #endif

   // set output
   out_BQBtmgC = BQBtmgC;
   out_aiD = pDiag;
}

// convenience function for breaking a Thyra matrix (possibly blocked)
// composed of Epetra_CrsMatrix objects.  This is primarily intended to
// be used internally by the library.
std::pair<int,int> thyraMatrixToCrsVector(const RCP<const Thyra::LinearOpBase<double> > & A,bool failOnNonZero,
                                          std::vector<Teuchos::RCP<const Epetra_CrsMatrix> > & blocks)
{
   using namespace Thyra;

   std::vector<RCP<const LinearOpBase<double> > > thyraBlocks;
   std::pair<int,int> shape;

   // figure out if A is a block matrix, fill temporary
   // vector with components
   /////////////////////////////////////////////////////

   RCP<const BlockedLinearOpBase<double> > blockOp = rcp_dynamic_cast<const BlockedLinearOpBase<double> >(A,false);
   if(blockOp!=Teuchos::null) {
      // this is a block matrix
      shape = std::make_pair(blockOp->productRange()->numBlocks(),
                             blockOp->productDomain()->numBlocks());

      // build a vector of Thyra blocks
      for(int i=0;i<shape.first;i++)
         for(int j=0;j<shape.second;j++)
            thyraBlocks.push_back(blockOp->getBlock(i,j));
   } 
   else {
      // this is not a block matrix
      shape = std::make_pair(1,1);

      thyraBlocks.push_back(A);
   }

   // things needed for getEpetraOpView
   EOpTransp eOpTransp;
   EApplyEpetraOpAs eOpApplyAs;
   EAdjointEpetraOp eOpAdjointSupport;

   // loop over Thyra blocks and convert them to Epetra_CrsMatrix or null
   std::vector<RCP<const LinearOpBase<double> > >::const_iterator itr; 
   for(itr=thyraBlocks.begin();itr!=thyraBlocks.end();++itr) {
      // we want to fill this
      RCP<const Epetra_CrsMatrix> eCrsOp;

      // try to get Thyra wrapped Epetra operator
      const RCP<const EpetraLinearOpBase> eThyraOp = rcp_dynamic_cast<const EpetraLinearOpBase>(*itr);

      // work out which type it is
      if(eThyraOp!=Teuchos::null) { 
         // this is an Epetra_Operator => it should be Epetra_CrsMatrix
         RCP<const Epetra_Operator> eGeneralOp;

         // extract Epetra_Operator (we are ignoring all the modifiers!)
         eThyraOp->getEpetraOpView(&eGeneralOp,&eOpTransp,&eOpApplyAs,&eOpAdjointSupport);

         // if failOnNonZero==false this may return null!
         eCrsOp = rcp_dynamic_cast<const Epetra_CrsMatrix>(eGeneralOp,failOnNonZero);
      }
      else { 
         // this is not an Epetra_Operator => it should be Thyra::zero
         rcp_dynamic_cast<const ZeroLinearOpBase<double> >(*itr);
         
         // if we get this far...this is a zero operator
         eCrsOp = Teuchos::null;
      }

      blocks.push_back(eCrsOp);
   }

   return shape;
}


} // end namespace Epetra
} // end namespace PB
