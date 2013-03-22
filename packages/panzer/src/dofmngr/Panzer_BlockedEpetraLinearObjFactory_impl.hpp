// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#ifndef PANZER_BLOCKED_EPETRA_LINEAR_OBJ_FACTORY_IMPL_HPP
#define PANZER_BLOCKED_EPETRA_LINEAR_OBJ_FACTORY_IMPL_HPP

#include "Panzer_UniqueGlobalIndexer.hpp"

#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"

#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

using Teuchos::RCP;

namespace panzer {

// ************************************************************
// class BlockedEpetraLinearObjFactory
// ************************************************************

template <typename Traits,typename LocalOrdinalT>
BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
BlockedEpetraLinearObjFactory(const Teuchos::RCP<const Epetra_Comm> & comm,
                              const Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,std::pair<int,int> > > & blkProvider,
                              const std::vector<Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > > & gidProviders)
   : blockProvider_(blkProvider), blockedDOFManager_(Teuchos::null), comm_(comm)
{ 
   gidProviders_ = gidProviders;

   makeRoomForBlocks(gidProviders_.size());

   // build and register the gather/scatter evaluators with 
   // the base class.
   this->buildGatherScatterEvaluators(*this);
}

template <typename Traits,typename LocalOrdinalT>
BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
BlockedEpetraLinearObjFactory(const Teuchos::RCP<const Epetra_Comm> & comm,
                              const Teuchos::RCP<const BlockedDOFManager<LocalOrdinalT,int> > & gidProvider)
   : blockProvider_(gidProvider), blockedDOFManager_(gidProvider), comm_(comm), rawMpiComm_(Teuchos::null)
{ 
   for(std::size_t i=0;i<gidProvider->getFieldDOFManagers().size();i++)
      gidProviders_.push_back(gidProvider->getFieldDOFManagers()[i]);

   makeRoomForBlocks(gidProviders_.size());

   // build and register the gather/scatter evaluators with 
   // the base class.
   this->buildGatherScatterEvaluators(*this);
}

template <typename Traits,typename LocalOrdinalT>
BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
BlockedEpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                              const Teuchos::RCP<const BlockedDOFManager<LocalOrdinalT,int> > & gidProvider)
   : blockProvider_(gidProvider), blockedDOFManager_(gidProvider), comm_(Teuchos::null), rawMpiComm_(comm->getRawMpiComm())
{ 
   comm_ = Teuchos::rcp(new Epetra_MpiComm(*rawMpiComm_));

   for(std::size_t i=0;i<gidProvider->getFieldDOFManagers().size();i++)
      gidProviders_.push_back(gidProvider->getFieldDOFManagers()[i]);

   makeRoomForBlocks(gidProviders_.size());

   // build and register the gather/scatter evaluators with 
   // the base class.
   this->buildGatherScatterEvaluators(*this);
}

template <typename Traits,typename LocalOrdinalT>
BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::~BlockedEpetraLinearObjFactory()
{ }

// LinearObjectFactory functions 
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<LinearObjContainer> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::buildLinearObjContainer() const
{
   std::vector<Teuchos::RCP<const Epetra_Map> > blockMaps;
   std::size_t blockDim = gidProviders_.size();
   for(std::size_t i=0;i<blockDim;i++)
      blockMaps.push_back(getMap(i));

   Teuchos::RCP<BlockedEpetraLinearObjContainer > container = Teuchos::rcp(new BlockedEpetraLinearObjContainer);
   container->setMapsForBlocks(blockMaps);

   return container;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<LinearObjContainer> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::buildGhostedLinearObjContainer() const
{
   std::vector<Teuchos::RCP<const Epetra_Map> > blockMaps;
   std::size_t blockDim = gidProviders_.size();
   for(std::size_t i=0;i<blockDim;i++)
      blockMaps.push_back(getGhostedMap(i));

   Teuchos::RCP<BlockedEpetraLinearObjContainer > container = Teuchos::rcp(new BlockedEpetraLinearObjContainer);
   container->setMapsForBlocks(blockMaps);

   return container;
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::globalToGhostContainer(const LinearObjContainer & in,
                                                                          LinearObjContainer & out,int mem) const
{
   using Teuchos::is_null;

   typedef LinearObjContainer LOC;
   typedef BlockedEpetraLinearObjContainer BLOC;
   const BLOC & b_in = Teuchos::dyn_cast<const BLOC>(in); 
   BLOC & b_out = Teuchos::dyn_cast<BLOC>(out); 
  
   // Operations occur if the GLOBAL container has the correct targets!
   // Users set the GLOBAL continer arguments
   if ( !is_null(b_in.get_x()) && !is_null(b_out.get_x()) && ((mem & LOC::X)==LOC::X))
     globalToGhostThyraVector(b_in.get_x(),b_out.get_x());
  
   if ( !is_null(b_in.get_dxdt()) && !is_null(b_out.get_dxdt()) && ((mem & LOC::DxDt)==LOC::DxDt))
     globalToGhostThyraVector(b_in.get_dxdt(),b_out.get_dxdt());

   if ( !is_null(b_in.get_f()) && !is_null(b_out.get_f()) && ((mem & LOC::F)==LOC::F))
      globalToGhostThyraVector(b_in.get_f(),b_out.get_f());
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::ghostToGlobalContainer(const LinearObjContainer & in,
                                                                          LinearObjContainer & out,int mem) const
{
   using Teuchos::is_null;

   typedef LinearObjContainer LOC;
   typedef BlockedEpetraLinearObjContainer BLOC;
   const BLOC & b_in = Teuchos::dyn_cast<const BLOC>(in); 
   BLOC & b_out = Teuchos::dyn_cast<BLOC>(out); 

   // Operations occur if the GLOBAL container has the correct targets!
   // Users set the GLOBAL continer arguments
   if ( !is_null(b_in.get_x()) && !is_null(b_out.get_x()) && ((mem & LOC::X)==LOC::X))
     ghostToGlobalThyraVector(b_in.get_x(),b_out.get_x());

   if ( !is_null(b_in.get_f()) && !is_null(b_out.get_f()) && ((mem & LOC::F)==LOC::F))
     ghostToGlobalThyraVector(b_in.get_f(),b_out.get_f());

   if ( !is_null(b_in.get_A()) && !is_null(b_out.get_A()) && ((mem & LOC::Mat)==LOC::Mat))
     ghostToGlobalThyraMatrix(*b_in.get_A(),*b_out.get_A());
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
adjustForDirichletConditions(const LinearObjContainer & localBCRows,
                             const LinearObjContainer & globalBCRows,
                             LinearObjContainer & ghostedObjs) const
{
   typedef BlockedEpetraLinearObjContainer BLOC;

   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::LinearOpBase;
   using Thyra::PhysicallyBlockedLinearOpBase;
   using Thyra::VectorBase;
   using Thyra::ProductVectorBase;
   using Thyra::get_Epetra_Vector;
   using Thyra::get_Epetra_Operator;

   std::size_t blockDim = gidProviders_.size();

   // first cast to block LOCs
   const BLOC & b_localBCRows = Teuchos::dyn_cast<const BLOC>(localBCRows); 
   const BLOC & b_globalBCRows = Teuchos::dyn_cast<const BLOC>(globalBCRows); 
   BLOC & b_ghosted = Teuchos::dyn_cast<BLOC>(ghostedObjs); 

   // TEUCHOS_ASSERT(b_ghosted.get_A()!=Teuchos::null);
   // TEUCHOS_ASSERT(b_ghosted.get_f()!=Teuchos::null);
   TEUCHOS_ASSERT(b_localBCRows.get_x()!=Teuchos::null);
   TEUCHOS_ASSERT(b_globalBCRows.get_x()!=Teuchos::null);

   // cast each component as needed to their product form
   RCP<PhysicallyBlockedLinearOpBase<double> > A = rcp_dynamic_cast<PhysicallyBlockedLinearOpBase<double> >(b_ghosted.get_A());
   RCP<ProductVectorBase<double> > f = rcp_dynamic_cast<ProductVectorBase<double> >(b_ghosted.get_f());
   RCP<ProductVectorBase<double> > local_bcs  = rcp_dynamic_cast<ProductVectorBase<double> >(b_localBCRows.get_x(),true);
   RCP<ProductVectorBase<double> > global_bcs = rcp_dynamic_cast<ProductVectorBase<double> >(b_globalBCRows.get_x(),true);

   // sanity check!
   if(A!=Teuchos::null) TEUCHOS_ASSERT(A->productRange()->numBlocks()==(int) blockDim);
   if(A!=Teuchos::null) TEUCHOS_ASSERT(A->productDomain()->numBlocks()==(int) blockDim);
   if(f!=Teuchos::null) TEUCHOS_ASSERT(f->productSpace()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(local_bcs->productSpace()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(global_bcs->productSpace()->numBlocks()==(int) blockDim);

   for(std::size_t i=0;i<blockDim;i++) {
      // grab epetra vector
      RCP<const Epetra_Vector> e_local_bcs = get_Epetra_Vector(*getGhostedMap(i),local_bcs->getVectorBlock(i));
      RCP<const Epetra_Vector> e_global_bcs = get_Epetra_Vector(*getGhostedMap(i),global_bcs->getVectorBlock(i));

      // pull out epetra values
      RCP<VectorBase<double> > th_f = (f==Teuchos::null) ? Teuchos::null : f->getNonconstVectorBlock(i);
      RCP<Epetra_Vector> e_f;
      if(th_f==Teuchos::null)
        e_f = Teuchos::null;
      else
        e_f = get_Epetra_Vector(*getGhostedMap(i),th_f);

      for(std::size_t j=0;j<blockDim;j++) {

         // pull out epetra values
         RCP<LinearOpBase<double> > th_A = (A== Teuchos::null)? Teuchos::null : A->getNonconstBlock(i,j);
 
         // don't do anyting if opertor is null
         RCP<Epetra_CrsMatrix> e_A;
         if(th_A==Teuchos::null)
            e_A = Teuchos::null;
         else 
            e_A = rcp_dynamic_cast<Epetra_CrsMatrix>(get_Epetra_Operator(*th_A),true);

         // adjust Block operator
         adjustForDirichletConditions(*e_local_bcs,*e_global_bcs,e_f.ptr(),e_A.ptr());

         e_f = Teuchos::null; // this is so we only adjust it once on the first pass
      }
   }
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
adjustForDirichletConditions(const Epetra_Vector & local_bcs,
                             const Epetra_Vector & global_bcs,
                             const Teuchos::Ptr<Epetra_Vector> & f,
                             const Teuchos::Ptr<Epetra_CrsMatrix> & A) const
{
   if(f==Teuchos::null && A==Teuchos::null)
      return;

   TEUCHOS_ASSERT(local_bcs.MyLength()==global_bcs.MyLength());
   for(int i=0;i<local_bcs.MyLength();i++) {
      if(global_bcs[i]==0.0)
         continue;

      int numEntries = 0;
      double * values = 0;
      int * indices = 0;

      if(local_bcs[i]==0.0) { 
         // this boundary condition was NOT set by this processor

         // if they exist put 0.0 in each entry
         if(!Teuchos::is_null(f))
            (*f)[i] = 0.0;
         if(!Teuchos::is_null(A)) {
            A->ExtractMyRowView(i,numEntries,values,indices);
            for(int c=0;c<numEntries;c++) 
               values[c] = 0.0;
         }
      }
      else {
         // this boundary condition was set by this processor

         double scaleFactor = global_bcs[i];

         // if they exist scale linear objects by scale factor
         if(!Teuchos::is_null(f))
            (*f)[i] /= scaleFactor;
         if(!Teuchos::is_null(A)) {
            A->ExtractMyRowView(i,numEntries,values,indices);
            for(int c=0;c<numEntries;c++) 
               values[c] /= scaleFactor;
         }
      }
   }
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::MpiComm<int> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getComm() const
{
   return Teuchos::MpiComm<int>(Teuchos::opaqueWrapper(dynamic_cast<const Epetra_MpiComm &>(*getEpetraComm()).Comm()));
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
initializeContainer(int mem,LinearObjContainer & loc) const
{
   typedef BlockedEpetraLinearObjContainer BLOC;

   BLOC & bloc = Teuchos::dyn_cast<BLOC>(loc);
   initializeContainer(mem,bloc);
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
initializeGhostedContainer(int mem,LinearObjContainer & loc) const
{
   typedef BlockedEpetraLinearObjContainer BLOC;

   BLOC & bloc = Teuchos::dyn_cast<BLOC>(loc);
   initializeGhostedContainer(mem,bloc);
}

// Generic methods 
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
initializeContainer(int mem,BlockedEpetraLinearObjContainer & loc) const
{
   typedef BlockedEpetraLinearObjContainer BLOC;

   loc.clear();

   if((mem & BLOC::X) == BLOC::X)
      loc.set_x(getThyraDomainVector());

   if((mem & BLOC::DxDt) == BLOC::DxDt)
      loc.set_dxdt(getThyraDomainVector());
    
   if((mem & BLOC::F) == BLOC::F)
      loc.set_f(getThyraRangeVector());

   if((mem & BLOC::Mat) == BLOC::Mat)
      loc.set_A(getThyraMatrix());
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
initializeGhostedContainer(int mem,BlockedEpetraLinearObjContainer & loc) const
{
   typedef BlockedEpetraLinearObjContainer BLOC;

   loc.clear();

   if((mem & BLOC::X) == BLOC::X)
      loc.set_x(getGhostedThyraDomainVector());

   if((mem & BLOC::DxDt) == BLOC::DxDt)
      loc.set_dxdt(getGhostedThyraDomainVector());
    
   if((mem & BLOC::F) == BLOC::F)
      loc.set_f(getGhostedThyraRangeVector());

   if((mem & BLOC::Mat) == BLOC::Mat)
      loc.set_A(getGhostedThyraMatrix());
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
addExcludedPair(int rowBlock,int colBlock)
{
   excludedPairs_.insert(std::make_pair(rowBlock,colBlock));
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
addExcludedPairs(const std::vector<std::pair<int,int> > & exPairs)
{
   for(std::size_t i=0;i<exPairs.size();i++)
      excludedPairs_.insert(exPairs[i]);
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGlobalIndexer(int i) const
{
   return gidProviders_[i];
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
makeRoomForBlocks(std::size_t blockCnt)
{
   maps_.resize(blockCnt); 
   ghostedMaps_.resize(blockCnt); 
   importers_.resize(blockCnt); 
   exporters_.resize(blockCnt); 
}

// Thyra methods 
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getThyraDomainSpace() const
{
   if(domainSpace_==Teuchos::null) {
      // loop over all vectors and build the vector space
      std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double> > > vsArray;
      for(std::size_t i=0;i<gidProviders_.size();i++)  
         vsArray.push_back(Thyra::create_VectorSpace(getMap(i)));

      domainSpace_ = Thyra::productVectorSpace<double>(vsArray);
   }
   
   return domainSpace_;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getThyraRangeSpace() const
{
   if(rangeSpace_==Teuchos::null) {
      // loop over all vectors and build the vector space
      std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double> > > vsArray;
      for(std::size_t i=0;i<gidProviders_.size();i++)  
         vsArray.push_back(Thyra::create_VectorSpace(getMap(i)));

      rangeSpace_ = Thyra::productVectorSpace<double>(vsArray);
   }
   
   return rangeSpace_;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Thyra::VectorBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getThyraDomainVector() const
{
   Teuchos::RCP<Thyra::VectorBase<double> > vec =
      Thyra::createMember<double>(*getThyraDomainSpace());
   Thyra::assign(vec.ptr(),0.0);

   Teuchos::RCP<Thyra::ProductVectorBase<double> > p_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(vec);
   for(std::size_t i=0;i<gidProviders_.size();i++) {
      TEUCHOS_ASSERT(Teuchos::rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(p_vec->getNonconstVectorBlock(i))->spmdSpace()->localSubDim()==getMap(i)->NumMyElements());
   }

   return vec;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Thyra::VectorBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getThyraRangeVector() const
{
   Teuchos::RCP<Thyra::VectorBase<double> > vec =
      Thyra::createMember<double>(*getThyraRangeSpace());
   Thyra::assign(vec.ptr(),0.0);

   return vec;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Thyra::LinearOpBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getThyraMatrix() const
{
   Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > blockedOp = Thyra::defaultBlockedLinearOp<double>();

   // get the block dimension
   std::size_t blockDim = gidProviders_.size();

   // this operator will be square
   blockedOp->beginBlockFill(blockDim,blockDim);

   // loop over each block
   for(std::size_t i=0;i<blockDim;i++) { 
      for(std::size_t j=0;j<blockDim;j++) {
         if(excludedPairs_.find(std::make_pair(i,j))==excludedPairs_.end()) {
            // build (i,j) block matrix and add it to blocked operator
            Teuchos::RCP<Thyra::LinearOpBase<double> > block = Thyra::nonconstEpetraLinearOp(getEpetraMatrix(i,j));
            blockedOp->setNonconstBlock(i,j,block);
         }
      }
   }

   // all done
   blockedOp->endBlockFill();

   return blockedOp;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Thyra::VectorSpaceBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGhostedThyraDomainSpace() const
{
   if(ghostedDomainSpace_==Teuchos::null) {
      // loop over all vectors and build the vector space
      std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double> > > vsArray;
      for(std::size_t i=0;i<gidProviders_.size();i++)  
         vsArray.push_back(Thyra::create_VectorSpace(getGhostedMap(i)));

      ghostedDomainSpace_ = Thyra::productVectorSpace<double>(vsArray);
   }
   
   return ghostedDomainSpace_;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Thyra::VectorSpaceBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGhostedThyraRangeSpace() const
{
   if(ghostedRangeSpace_==Teuchos::null) {
      // loop over all vectors and build the vector space
      std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double> > > vsArray;
      for(std::size_t i=0;i<gidProviders_.size();i++)  
         vsArray.push_back(Thyra::create_VectorSpace(getGhostedMap(i)));

      ghostedRangeSpace_ = Thyra::productVectorSpace<double>(vsArray);
   }
   
   return ghostedRangeSpace_;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Thyra::VectorBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGhostedThyraDomainVector() const
{
   Teuchos::RCP<Thyra::VectorBase<double> > vec =
      Thyra::createMember<double>(*getGhostedThyraDomainSpace());
   Thyra::assign(vec.ptr(),0.0);

   return vec;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Thyra::VectorBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGhostedThyraRangeVector() const
{
   Teuchos::RCP<Thyra::VectorBase<double> > vec =
      Thyra::createMember<double>(*getGhostedThyraRangeSpace());
   Thyra::assign(vec.ptr(),0.0);
 
   return vec;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Thyra::BlockedLinearOpBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGhostedThyraMatrix() const
{
   Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > blockedOp = Thyra::defaultBlockedLinearOp<double>();

   // get the block dimension
   std::size_t blockDim = gidProviders_.size();

   // this operator will be square
   blockedOp->beginBlockFill(blockDim,blockDim);

   // loop over each block
   for(std::size_t i=0;i<blockDim;i++) { 
      for(std::size_t j=0;j<blockDim;j++) {
         if(excludedPairs_.find(std::make_pair(i,j))==excludedPairs_.end()) {
            // build (i,j) block matrix and add it to blocked operator
            Teuchos::RCP<Thyra::LinearOpBase<double> > block = Thyra::nonconstEpetraLinearOp(getGhostedEpetraMatrix(i,j));
            blockedOp->setNonconstBlock(i,j,block);
         }
      }
   }

   // all done
   blockedOp->endBlockFill();

   return blockedOp;
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
ghostToGlobalThyraVector(const Teuchos::RCP<const Thyra::VectorBase<double> > & in,
                         const Teuchos::RCP<Thyra::VectorBase<double> > & out) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;
   using Thyra::get_Epetra_Vector;

   std::size_t blockDim = gidProviders_.size();

   // get product vectors
   RCP<const ProductVectorBase<double> > prod_in = rcp_dynamic_cast<const ProductVectorBase<double> >(in,true);
   RCP<ProductVectorBase<double> > prod_out      = rcp_dynamic_cast<ProductVectorBase<double> >(out,true);

   TEUCHOS_ASSERT(prod_in->productSpace()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(prod_out->productSpace()->numBlocks()==(int) blockDim);

   for(std::size_t i=0;i<blockDim;i++) {
      // first get each Epetra vector
      RCP<const Epetra_Vector> ep_in = get_Epetra_Vector(*getGhostedMap(i),prod_in->getVectorBlock(i));
      RCP<Epetra_Vector> ep_out      = get_Epetra_Vector(*getMap(i),prod_out->getNonconstVectorBlock(i));

      // use Epetra to do global communication
      ghostToGlobalEpetraVector(i,*ep_in,*ep_out);
   }
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
ghostToGlobalThyraMatrix(const Thyra::LinearOpBase<double> & in,Thyra::LinearOpBase<double> & out) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Teuchos::dyn_cast;
   using Thyra::LinearOpBase;
   using Thyra::PhysicallyBlockedLinearOpBase;
   using Thyra::get_Epetra_Operator;

   std::size_t blockDim = gidProviders_.size();

   // get product vectors
   const PhysicallyBlockedLinearOpBase<double> & prod_in = dyn_cast<const PhysicallyBlockedLinearOpBase<double> >(in);
   PhysicallyBlockedLinearOpBase<double> & prod_out      = dyn_cast<PhysicallyBlockedLinearOpBase<double> >(out);

   TEUCHOS_ASSERT(prod_in.productRange()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(prod_in.productDomain()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(prod_out.productRange()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(prod_out.productDomain()->numBlocks()==(int) blockDim);

   for(std::size_t i=0;i<blockDim;i++) {
      for(std::size_t j=0;j<blockDim;j++) {
         if(excludedPairs_.find(std::make_pair(i,j))==excludedPairs_.end()) {
            // extract the blocks
            RCP<const LinearOpBase<double> > th_in = prod_in.getBlock(i,j);
            RCP<LinearOpBase<double> > th_out = prod_out.getNonconstBlock(i,j);
   
            // sanity check
            TEUCHOS_ASSERT(th_in!=Teuchos::null);
            TEUCHOS_ASSERT(th_out!=Teuchos::null);
   
            // get the epetra version of the blocks
            RCP<const Epetra_CrsMatrix> ep_in = rcp_dynamic_cast<const Epetra_CrsMatrix>(get_Epetra_Operator(*th_in),true);
            RCP<Epetra_CrsMatrix> ep_out      = rcp_dynamic_cast<Epetra_CrsMatrix>(get_Epetra_Operator(*th_out),true);
   
            // use Epetra to do global communication
            ghostToGlobalEpetraMatrix(i,*ep_in,*ep_out);
         }
      }
   }
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
globalToGhostThyraVector(const Teuchos::RCP<const Thyra::VectorBase<double> > & in,
                         const Teuchos::RCP<Thyra::VectorBase<double> > & out) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;
   using Thyra::get_Epetra_Vector;

   std::size_t blockDim = gidProviders_.size();

   // get product vectors
   RCP<const ProductVectorBase<double> > prod_in = rcp_dynamic_cast<const ProductVectorBase<double> >(in,true);
   RCP<ProductVectorBase<double> > prod_out      = rcp_dynamic_cast<ProductVectorBase<double> >(out,true);

   TEUCHOS_ASSERT(prod_in->productSpace()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(prod_out->productSpace()->numBlocks()==(int) blockDim);

   for(std::size_t i=0;i<blockDim;i++) {
      // first get each Epetra vector
      RCP<const Epetra_Vector> ep_in = get_Epetra_Vector(*getMap(i),prod_in->getVectorBlock(i));
      RCP<Epetra_Vector> ep_out      = get_Epetra_Vector(*getGhostedMap(i),prod_out->getNonconstVectorBlock(i));

      // use Epetra to do global communication
      globalToGhostEpetraVector(i,*ep_in,*ep_out);
   }
}

// Epetra methods 
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
ghostToGlobalEpetraVector(int i,const Epetra_Vector & in,Epetra_Vector & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<Epetra_Export> exporter = getGhostedExport(i);
   out.PutScalar(0.0);
   int err = out.Export(in,*exporter,Add);
   TEUCHOS_ASSERT_EQUALITY(err,0);
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
ghostToGlobalEpetraMatrix(int blockRow,const Epetra_CrsMatrix & in,Epetra_CrsMatrix & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<Epetra_Export> exporter = getGhostedExport(blockRow);
   out.PutScalar(0.0);
   int err = out.Export(in,*exporter,Add);
   TEUCHOS_ASSERT_EQUALITY(err,0);
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
globalToGhostEpetraVector(int i,const Epetra_Vector & in,Epetra_Vector & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<Epetra_Import> importer = getGhostedImport(i);
   out.PutScalar(0.0);
   int err = out.Import(in,*importer,Insert);
   TEUCHOS_ASSERT_EQUALITY(err,0);
}

// get the map from the matrix
template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getMap(int i) const
{
   if(maps_[i]==Teuchos::null) 
      maps_[i] = buildEpetraMap(i);

   return maps_[i];
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGhostedMap(int i) const
{
   if(ghostedMaps_[i]==Teuchos::null) 
      ghostedMaps_[i] = buildEpetraGhostedMap(i);

   return ghostedMaps_[i];
}

// get the graph of the crs matrix
template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGraph(int i,int j) const
{
   typedef boost::unordered_map<std::pair<int,int>,Teuchos::RCP<Epetra_CrsGraph> > GraphMap;
   
   GraphMap::const_iterator itr = graphs_.find(std::make_pair(i,j));
   Teuchos::RCP<Epetra_CrsGraph> graph;
   if(itr==graphs_.end()) {
      graph = buildEpetraGraph(i,j);
      graphs_[std::make_pair(i,j)] = graph;
   }
   else
      graph = itr->second;

   TEUCHOS_ASSERT(graph!=Teuchos::null);
   return graph;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGhostedGraph(int i,int j) const
{
   typedef boost::unordered_map<std::pair<int,int>,Teuchos::RCP<Epetra_CrsGraph> > GraphMap;
   
   GraphMap::const_iterator itr = ghostedGraphs_.find(std::make_pair(i,j));
   Teuchos::RCP<Epetra_CrsGraph> ghostedGraph;
   if(itr==ghostedGraphs_.end()) {
      ghostedGraph = buildEpetraGhostedGraph(i,j);
      ghostedGraphs_[std::make_pair(i,j)] = ghostedGraph;
   }
   else
      ghostedGraph = itr->second;

   TEUCHOS_ASSERT(ghostedGraph!=Teuchos::null);
   return ghostedGraph;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Import> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGhostedImport(int i) const
{
   if(importers_[i]==Teuchos::null)
      importers_[i] = Teuchos::rcp(new Epetra_Import(*getGhostedMap(i),*getMap(i)));

   return importers_[i];
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Export> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGhostedExport(int i) const
{
   if(exporters_[i]==Teuchos::null)
      exporters_[i] = Teuchos::rcp(new Epetra_Export(*getGhostedMap(i),*getMap(i)));

   return exporters_[i];
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
buildEpetraMap(int i) const
{
   Teuchos::RCP<Epetra_Map> map; // result
   std::vector<int> indices;

   // get the global indices
   getGlobalIndexer(i)->getOwnedIndices(indices);

   map = Teuchos::rcp(new Epetra_Map(-1,indices.size(),&indices[0],0,*comm_));

   return map;
}

// build the ghosted map
template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
buildEpetraGhostedMap(int i) const
{
   std::vector<int> indices;

   // get the global indices
   getGlobalIndexer(i)->getOwnedAndSharedIndices(indices);

   return Teuchos::rcp(new Epetra_Map(-1,indices.size(),&indices[0],0,*comm_));
}

// get the graph of the crs matrix
template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
buildEpetraGraph(int i,int j) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // build the map and allocate the space for the graph and
   // grab the ghosted graph
   RCP<Epetra_Map> map_i = getMap(i);
   RCP<Epetra_Map> map_j = getMap(j);
   RCP<Epetra_CrsGraph> graph  = rcp(new Epetra_CrsGraph(Copy,*map_i,0));
   RCP<Epetra_CrsGraph> oGraph = getGhostedGraph(i,j);

   // perform the communication to finish building graph
   RCP<Epetra_Export> exporter = getGhostedExport(i);
   int err = graph->Export( *oGraph, *exporter, Insert );
   TEUCHOS_ASSERT_EQUALITY(err,0);
   graph->FillComplete(*map_j,*map_i);
   graph->OptimizeStorage();
  
   return graph;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
buildEpetraGhostedGraph(int i,int j) const
{
   // build the map and allocate the space for the graph
   Teuchos::RCP<Epetra_Map> rowMap = getGhostedMap(i);
   Teuchos::RCP<Epetra_Map> colMap = getGhostedMap(j);
   Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*rowMap,*colMap,0));

   std::vector<std::string> elementBlockIds;
   
   Teuchos::RCP<const UniqueGlobalIndexer<LocalOrdinalT,int> > rowProvider, colProvider;
 
   rowProvider = getGlobalIndexer(i);
   colProvider = getGlobalIndexer(j);

   blockProvider_->getElementBlockIds(elementBlockIds); // each sub provider "should" have the
                                                        // same element blocks

   // graph information about the mesh
   std::vector<std::string>::const_iterator blockItr;
   for(blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;

      // grab elements for this block
      const std::vector<LocalOrdinalT> & elements = blockProvider_->getElementBlock(blockId); // each sub provider "should" have the
                                                                                              // same elements in each element block

      // get information about number of indicies
      std::vector<int> row_gids;
      std::vector<int> col_gids;

      // loop over the elemnts
      for(std::size_t elmt=0;elmt<elements.size();elmt++) {

         rowProvider->getElementGIDs(elements[elmt],row_gids);
         colProvider->getElementGIDs(elements[elmt],col_gids);
         for(std::size_t row=0;row<row_gids.size();row++)
            graph->InsertGlobalIndices(row_gids[row],col_gids.size(),&col_gids[0]);
      }
   }

   // finish filling the graph: Make sure the colmap and row maps coincide to 
   //                           minimize calls to LID lookups
   graph->FillComplete(*colMap,*rowMap);
   graph->OptimizeStorage();

   return graph;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Epetra_CrsMatrix> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getEpetraMatrix(int i,int j) const
{
   Teuchos::RCP<Epetra_CrsGraph> eGraph = getGraph(i,j);
   Teuchos::RCP<Epetra_CrsMatrix> mat = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *eGraph));
   TEUCHOS_ASSERT(mat->Filled());
   return mat;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<Epetra_CrsMatrix> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGhostedEpetraMatrix(int i,int j) const
{
   Teuchos::RCP<Epetra_CrsGraph> eGraph = getGhostedGraph(i,j); 
   Teuchos::RCP<Epetra_CrsMatrix> mat = Teuchos::rcp(new Epetra_CrsMatrix(Copy, *eGraph));
   TEUCHOS_ASSERT(mat->Filled());
   return mat;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<const Epetra_Comm> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getEpetraComm() const
{
   return comm_;
}

template <typename Traits,typename LocalOrdinalT>
int BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getBlockRowCount() const
{
   return gidProviders_.size();
}

template <typename Traits,typename LocalOrdinalT>
int BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getBlockColCount() const
{
   return gidProviders_.size();
}

}

#endif
