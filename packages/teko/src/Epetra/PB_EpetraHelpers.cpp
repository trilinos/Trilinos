#include "PB_EpetraHelpers.hpp"

// Thyra Includes
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_DefaultMultipliedLinearOp.hpp"
#include "Thyra_DefaultDiagonalLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_SpmdVectorSpaceBase.hpp"

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

// build an epetra operator from a block 2x2 matrix
Epetra_Operator * block2x2(const Epetra_Operator * A,const Epetra_Operator * B,
                         const Epetra_Operator * C,const Epetra_Operator * D,const std::string & str)
{
   using Teuchos::rcpFromRef;
   using Thyra::zero;

   // not all things can be zero
   assert(not (A==0 && B==0));
   assert(not (A==0 && C==0));
   assert(not (B==0 && D==0));
   assert(not (C==0 && D==0));


   RCP<const Thyra::LinearOpBase<double> > ptrA = A!=0 ? Thyra::epetraLinearOp(rcp(A,false),str+"_00"): null;
   RCP<const Thyra::LinearOpBase<double> > ptrB = B!=0 ? Thyra::epetraLinearOp(rcp(B,false),str+"_01"): null;
   RCP<const Thyra::LinearOpBase<double> > ptrC = C!=0 ? Thyra::epetraLinearOp(rcp(C,false),str+"_10"): null;
   RCP<const Thyra::LinearOpBase<double> > ptrD = D!=0 ? Thyra::epetraLinearOp(rcp(D,false),str+"_11"): null;

   // get some ranges so the null matrices can be set to zero
   const RCP<const Thyra::VectorSpaceBase<double> > row0Range  = A!=0 ? ptrA->range()  : ptrB->range();
   const RCP<const Thyra::VectorSpaceBase<double> > row1Range  = C!=0 ? ptrC->range()  : ptrD->range();
   const RCP<const Thyra::VectorSpaceBase<double> > col0Domain = A!=0 ? ptrA->domain() : ptrC->domain();
   const RCP<const Thyra::VectorSpaceBase<double> > col1Domain = B!=0 ? ptrB->domain() : ptrD->domain();

   // set previously null pointers to 0
   if(A==0) ptrA = Thyra::zero(row0Range,col0Domain);
   if(B==0) ptrB = Thyra::zero(row0Range,col1Domain);
   if(C==0) ptrC = Thyra::zero(row1Range,col0Domain);
   if(D==0) ptrD = Thyra::zero(row1Range,col1Domain);

   const RCP<const Thyra::LinearOpBase<double> > mat = Thyra::block2x2(ptrA,ptrB,ptrC,ptrD,str);

   return new PB::Epetra::EpetraOperatorWrapper(mat);
}

// Convert an Epetra "inverse" operator to a "forward" operator
Epetra_Operator * mechanicalInverse(const Epetra_Operator * inverse)
{
   // convert to a forward operator
   Teuchos::RCP<const Epetra_Operator> opAr[] = {rcp(inverse,false)};
   Teuchos::ETransp epetraOpsTransp[] = { Teuchos::NO_TRANS };
   EpetraExt::ProductOperator::EApplyMode epetraOpsApplyMode[] = { EpetraExt::ProductOperator::APPLY_MODE_APPLY_INVERSE };

   return new EpetraExt::ProductOperator(1,opAr,epetraOpsTransp,epetraOpsApplyMode);
}

// here we assume that Epetra_Vector is blocked and each block is distributed across all processors
void epetraToThyra(const Epetra_MultiVector & e,const Teuchos::Ptr<Thyra::VectorBase<double> > & t)
{
   using namespace Thyra;
   using namespace Teuchos;

   TEST_FOR_EXCEPTION(e.NumVectors()!=1,std::runtime_error,
                      "epetraToThyra not available for Epetra_MultiVector with dimension greater than 1");

   const double * epetraData = e[0];

   // dynamic cast to product space
   ProductVectorBase<double> & t_prod = dyn_cast<ProductVectorBase<double> >(*t);
   
   // get some information about the blocks
   int numBlocks = rcp_dynamic_cast<const ProductVectorSpaceBase<double> >(t_prod.space())->numBlocks();
 
   // loop over blocks and extract what is needed from the Epetra_Vector
   int blockOffset = 0;
   for(int i=0;i<numBlocks;i++) { 
      int stride = 0;
      double * localData = 0;

      // get sub-block vector and space
      const RCP<SpmdVectorBase<double> > blk = rcp_dynamic_cast<SpmdVectorBase<double> >(t_prod.getNonconstVectorBlock(i),true);
      const RCP<const SpmdVectorSpaceBase<double> > blk_spc = blk->spmdSpace();

      // get information about blocks
      int localSubDim = blk_spc->localSubDim();
      int localOffset = blk_spc->localOffset();

      // get local data vector
      blk->getLocalData(&localData,&stride);

      // perform copy
      for(int j=0;j<localSubDim;j+=stride) {
         localData[j] = epetraData[j+blockOffset];
      }

      // commit changes
      blk->commitLocalData(localData);       

      blockOffset += localSubDim; 
   }
}

// here we assume that Epetra_Vector is blocked and each block is distributed across all processors
void thyraToEpetra(const Teuchos::RCP<const Thyra::VectorBase<double> > & t, Epetra_MultiVector & e)
{
   using namespace Thyra;
   using namespace Teuchos;

   TEST_FOR_EXCEPTION(e.NumVectors()!=1,std::runtime_error,
                      "epetraToThyra not available for Epetra_MultiVector with dimension greater than 1");

   double * epetraData = e[0];

   const RCP<Teuchos::FancyOStream> out = Teuchos::VerboseObjectBase::getDefaultOStream();

   // dynamic cast to product space
   const RCP<const ProductVectorBase<double> > t_prod = rcp_dynamic_cast<const ProductVectorBase<double> >(t,true);
   
   // get some information about the blocks
   int numBlocks = rcp_dynamic_cast<const ProductVectorSpaceBase<double> >(t_prod->space())->numBlocks();
 
   // loop over blocks and extract what is needed from the Epetra_Vector
   int blockOffset = 0;
   for(int i=0;i<numBlocks;i++) { 
      int stride = 0;
      const double * localData = 0;

      // get sub-block vector and space
      const RCP<const SpmdVectorBase<double> > blk = rcp_dynamic_cast<const SpmdVectorBase<double> >(t_prod->getVectorBlock(i),true);
      const RCP<const SpmdVectorSpaceBase<double> > blk_spc = blk->spmdSpace();

      // get information about blocks
      int localSubDim = blk_spc->localSubDim();
      int localOffset = blk_spc->localOffset();

      // get local data vector
      blk->getLocalData(&localData,&stride);       

      // perform copy
      for(int j=0;j<localSubDim;j+=stride)
         epetraData[j+blockOffset] = localData[j];

      blockOffset += localSubDim; 
   }
}

} // end namespace Epetra
} // end namespace PB
