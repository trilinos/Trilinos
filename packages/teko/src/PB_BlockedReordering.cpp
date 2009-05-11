#include <iostream>

#include "PB_BlockedReordering.hpp"
#include "PB_Utilities.hpp"

#include "Thyra_DefaultProductMultiVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::rcp_dynamic_cast;
using Teuchos::Array;

namespace PB {

void BlockReorderManager::SetBlock(int blockIndex,int reorder)
{ 
   TEUCHOS_ASSERT(blockIndex<children_.size());

   RCP<BlockReorderManager> child = rcp(new BlockReorderLeaf(reorder));

   children_[blockIndex] = child;
}

const Teuchos::RCP<BlockReorderManager> BlockReorderManager::GetBlock(int blockIndex)
{
   TEUCHOS_ASSERT(blockIndex<children_.size());

   if(children_[blockIndex]==Teuchos::null)
      children_[blockIndex] = rcp(new BlockReorderManager());

   return children_[blockIndex];
}

const Teuchos::RCP<const BlockReorderManager> BlockReorderManager::GetBlock(int blockIndex) const
{
   TEUCHOS_ASSERT(blockIndex<children_.size());

   return children_[blockIndex];
}

std::string BlockReorderManager::toString() const
{
   // build the string by recursively calling each child
   std::stringstream ss;
   ss << "[";
   for(int i=0;i<children_.size();i++) {
      if(children_[i]==Teuchos::null) 
         ss << " <NULL> ";
      else
         ss << " " << children_[i]->toString() << " "; 
   }
   ss << "]";

   return ss.str();
}

int BlockReorderManager::LargestIndex() const
{
   int max = 0;
   for(int i=0;i<children_.size();i++) {
      // see if current child is larger 
      if(children_[i]!=Teuchos::null) {
         int subMax = children_[i]->LargestIndex();
         max = max > subMax ? max : subMax;
      }
   }

   return max;
} 

Teuchos::RCP<const Thyra::LinearOpBase<double> >
buildReorderedLinearOp(const BlockReorderManager & bmm,
                       const Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> > & blkOp)
{
   return buildReorderedLinearOp(bmm,bmm,blkOp);
}

Teuchos::RCP<const Thyra::LinearOpBase<double> >
buildReorderedLinearOp(const BlockReorderManager & rowMgr,const BlockReorderManager & colMgr,
                       const Teuchos::RCP<const Thyra::BlockedLinearOpBase<double> > & blkOp)
{
   typedef RCP<const BlockReorderManager> BRMptr;

   int rowSz = rowMgr.GetNumBlocks();
   int colSz = colMgr.GetNumBlocks();

   if(rowSz==0 && colSz==0) {
      // both are leaf nodes
      const BlockReorderLeaf & rowLeaf = dynamic_cast<const BlockReorderLeaf &>(rowMgr);
      const BlockReorderLeaf & colLeaf = dynamic_cast<const BlockReorderLeaf &>(colMgr);

      // simply return entry in matrix
      return blkOp->getBlock(rowLeaf.GetIndex(),colLeaf.GetIndex());
   }
   else if(rowSz==0) {
      // only row is a leaf node
      const BlockReorderLeaf & rowLeaf = dynamic_cast<const BlockReorderLeaf &>(rowMgr);

      PB::BlockedLinearOp reBlkOp = PB::createBlockedOp();

      // operator will be rowSz by colSz
      reBlkOp->beginBlockFill(1,colSz);   

      // fill the column entries
      for(int col=0;col<colSz;col++) {
         BRMptr colPtr = colMgr.GetBlock(col);

         reBlkOp->setBlock(0,col,buildReorderedLinearOp(rowMgr,*colPtr,blkOp));
      } 

      // done building
      reBlkOp->endBlockFill();   

      return reBlkOp;
   }
   else if(colSz==0) {
      // only row is a leaf node
      const BlockReorderLeaf & colLeaf = dynamic_cast<const BlockReorderLeaf &>(colMgr);

      PB::BlockedLinearOp reBlkOp = PB::createBlockedOp();

      // operator will be rowSz by colSz
      reBlkOp->beginBlockFill(rowSz,1);   

      // fill the row entries
      for(int row=0;row<rowSz;row++) {
         BRMptr rowPtr = rowMgr.GetBlock(row);

         reBlkOp->setBlock(row,0,buildReorderedLinearOp(*rowPtr,colMgr,blkOp));
      } 

      // done building
      reBlkOp->endBlockFill();   

      return reBlkOp;
   }
   else {
      PB::BlockedLinearOp reBlkOp = PB::createBlockedOp();
  
      // this is the general case
      TEUCHOS_ASSERT(rowSz>0);
      TEUCHOS_ASSERT(colSz>0);

      // operator will be rowSz by colSz
      reBlkOp->beginBlockFill(rowSz,colSz);   
   
      for(int row=0;row<rowSz;row++) {
         BRMptr rowPtr = rowMgr.GetBlock(row);
   
         for(int col=0;col<colSz;col++) {
            BRMptr colPtr = colMgr.GetBlock(col);
   
            reBlkOp->setBlock(row,col,buildReorderedLinearOp(*rowPtr,*colPtr,blkOp));
         } 
      }

      // done building
      reBlkOp->endBlockFill();   

      return reBlkOp;
   }
}

/** \brief Convert a flat multi vector into a reordered multivector.
  *
  * Convert a flat multi vector into a reordered multivector.
  */
Teuchos::RCP<Thyra::MultiVectorBase<double> >
buildReorderedMultiVector(const BlockReorderManager & mgr,
                          const Teuchos::RCP<Thyra::ProductMultiVectorBase<double> > & blkVec)
{
   using Teuchos::rcp_const_cast;

   // give vector const so that the right function is called
   const Teuchos::RCP<const Thyra::ProductMultiVectorBase<double> > blkVecConst
      = rcp_const_cast<const Thyra::ProductMultiVectorBase<double> >(blkVec);
   
   // its not really const, so take it away
   const Teuchos::RCP<Thyra::MultiVectorBase<double> > result
         = rcp_const_cast<Thyra::MultiVectorBase<double> >(buildReorderedMultiVector(mgr,blkVecConst));

   return result; 
}

/** \brief Convert a flat multi vector into a reordered multivector.
  *
  * Convert a flat multi vector into a reordered multivector.
  */
Teuchos::RCP<const Thyra::MultiVectorBase<double> >
buildReorderedMultiVector(const BlockReorderManager & mgr,
                          const Teuchos::RCP<const Thyra::ProductMultiVectorBase<double> > & blkVec)
{
   typedef RCP<const BlockReorderManager> BRMptr;

   int sz = mgr.GetNumBlocks();

   if(sz==0) {
      // its a  leaf nodes
      const BlockReorderLeaf & leaf = dynamic_cast<const BlockReorderLeaf &>(mgr);

      // simply return entry in matrix
      return blkVec->getMultiVectorBlock(leaf.GetIndex());
   } 
   else {
      Array<RCP<const Thyra::MultiVectorBase<double> > > multiVecs;
      Array<RCP<const Thyra::VectorSpaceBase<double> > > vecSpaces;

      // loop over each row
      for(int i=0;i<sz;i++) {
         BRMptr blkMgr = mgr.GetBlock(i);

         const RCP<const Thyra::MultiVectorBase<double> > lmv = buildReorderedMultiVector(*blkMgr,blkVec);
         const RCP<const Thyra::VectorSpaceBase<double> > lvs = lmv->range();

         multiVecs.push_back(lmv);
         vecSpaces.push_back(lvs);
      }

      // build a vector space
      const RCP<const Thyra::DefaultProductVectorSpace<double> > vs 
            = Thyra::productVectorSpace<double>(vecSpaces);

      // build the vector
      return Thyra::defaultProductMultiVector<double>(vs,multiVecs);
   }
}

/** Helper function to assist with the non-constant
  * version of buildFlatMultiVector.
  */
void buildNonconstFlatMultiVector(const BlockReorderManager & mgr,
                                  const RCP<Thyra::MultiVectorBase<double> > & blkVec,
                                  Array<RCP<Thyra::MultiVectorBase<double> > > & multivecs,
                                  Array<RCP<const Thyra::VectorSpaceBase<double> > > & vecspaces)
{
   typedef RCP<const BlockReorderManager> BRMptr;

   int sz = mgr.GetNumBlocks();

   if(sz==0) {
      // its a  leaf nodes
      const BlockReorderLeaf & leaf = dynamic_cast<const BlockReorderLeaf &>(mgr);
      int index = leaf.GetIndex();

      // simply return entry in matrix
      multivecs[index] = blkVec;
      vecspaces[index] = blkVec->range();
   } 
   else {
      const RCP<Thyra::ProductMultiVectorBase<double> > prodMV
            = rcp_dynamic_cast<Thyra::ProductMultiVectorBase<double> >(blkVec);
 
      // get flattened elements from each child
      for(int i=0;i<sz;i++) {
         const RCP<Thyra::MultiVectorBase<double> > mv = prodMV->getNonconstMultiVectorBlock(i);
         buildNonconstFlatMultiVector(*mgr.GetBlock(i),mv,multivecs,vecspaces);
      }
   }
   
}

/** Helper function to assist with the function
  * of the same name.
  */
void buildFlatMultiVector(const BlockReorderManager & mgr,
                          const RCP<const Thyra::MultiVectorBase<double> > & blkVec,
                          Array<RCP<const Thyra::MultiVectorBase<double> > > & multivecs,
                          Array<RCP<const Thyra::VectorSpaceBase<double> > > & vecspaces)
{
   typedef RCP<const BlockReorderManager> BRMptr;

   int sz = mgr.GetNumBlocks();

   if(sz==0) {
      // its a  leaf nodes
      const BlockReorderLeaf & leaf = dynamic_cast<const BlockReorderLeaf &>(mgr);
      int index = leaf.GetIndex();

      // simply return entry in matrix
      multivecs[index] = blkVec;
      vecspaces[index] = blkVec->range();
   } 
   else {
      const RCP<const Thyra::ProductMultiVectorBase<double> > prodMV
            = rcp_dynamic_cast<const Thyra::ProductMultiVectorBase<double> >(blkVec);
 
      // get flattened elements from each child
      for(int i=0;i<sz;i++) {
         RCP<const Thyra::MultiVectorBase<double> > mv = prodMV->getMultiVectorBlock(i);
         buildFlatMultiVector(*mgr.GetBlock(i),mv,multivecs,vecspaces);
      }
   }
   
}

/** \brief Convert a reordered multivector into a flat multivector.
  *
  * Convert a reordered multivector into a flat multivector.
  */
Teuchos::RCP<Thyra::MultiVectorBase<double> >
buildFlatMultiVector(const BlockReorderManager & mgr,
                     const Teuchos::RCP<Thyra::ProductMultiVectorBase<double> > & blkVec)
{
   int numBlocks = mgr.LargestIndex()+1;
 
   Array<RCP<Thyra::MultiVectorBase<double> > > multivecs(numBlocks);
   Array<RCP<const Thyra::VectorSpaceBase<double> > > vecspaces(numBlocks);

   // flatten everything into a vector first
   buildNonconstFlatMultiVector(mgr,blkVec,multivecs,vecspaces);

   // build a vector space
   const RCP<Thyra::DefaultProductVectorSpace<double> > vs 
         = Thyra::productVectorSpace<double>(vecspaces);

   // build the vector
   return Thyra::defaultProductMultiVector<double>(vs,multivecs);
}

/** \brief Convert a reordered multivector into a flat multivector.
  *
  * Convert a reordered multivector into a flat multivector.
  */
Teuchos::RCP<const Thyra::MultiVectorBase<double> >
buildFlatMultiVector(const BlockReorderManager & mgr,
                     const Teuchos::RCP<const Thyra::ProductMultiVectorBase<double> > & blkVec)
{
   int numBlocks = mgr.LargestIndex()+1;
 
   Array<RCP<const Thyra::MultiVectorBase<double> > > multivecs(numBlocks);
   Array<RCP<const Thyra::VectorSpaceBase<double> > > vecspaces(numBlocks);

   // flatten everything into a vector first
   buildFlatMultiVector(mgr,blkVec,multivecs,vecspaces);

   // build a vector space
   const RCP<const Thyra::DefaultProductVectorSpace<double> > vs 
         = Thyra::productVectorSpace<double>(vecspaces);

   // build the vector
   return Thyra::defaultProductMultiVector<double>(vs,multivecs);
}

} // end namespace PB
