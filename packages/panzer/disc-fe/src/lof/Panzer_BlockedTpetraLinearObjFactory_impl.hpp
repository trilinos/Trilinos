// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_BlockedTpetraLinearObjFactory_impl_hpp__
#define   __Panzer_BlockedTpetraLinearObjFactory_impl_hpp__

// Panzer
#include "Panzer_BlockedVector_ReadOnly_GlobalEvaluationData.hpp"
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_EpetraVector_Write_GlobalEvaluationData.hpp"                    // JMG:  Remove this eventually.
#endif
#include "Panzer_TpetraVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_GlobalIndexer.hpp"

// Thyra
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Thyra_TpetraThyraWrappers.hpp"

// Tpetra
#include "Tpetra_CrsMatrix.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_Vector.hpp"

namespace panzer {

using Teuchos::RCP;

// ************************************************************
// class BlockedTpetraLinearObjFactory
// ************************************************************

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
BlockedTpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                              const Teuchos::RCP<const BlockedDOFManager> & gidProvider)
   : blockProvider_(gidProvider), blockedDOFManager_(gidProvider), comm_(comm)
{
  for(std::size_t i=0;i<gidProvider->getFieldDOFManagers().size();i++)
    gidProviders_.push_back(gidProvider->getFieldDOFManagers()[i]);

  makeRoomForBlocks(gidProviders_.size());

  // build and register the gather/scatter evaluators with
  // the base class.
  this->buildGatherScatterEvaluators(*this);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
BlockedTpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                              const std::vector<Teuchos::RCP<const panzer::GlobalIndexer>> & gidProviders)
  : gidProviders_(gidProviders), comm_(comm)
{
  makeRoomForBlocks(gidProviders_.size());
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
~BlockedTpetraLinearObjFactory()
{ }

// LinearObjectFactory functions
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<LinearObjContainer>
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildLinearObjContainer() const
{
   std::vector<Teuchos::RCP<const MapType> > blockMaps;
   std::size_t blockDim = gidProviders_.size();
   for(std::size_t i=0;i<blockDim;i++)
      blockMaps.push_back(getMap(i));

   Teuchos::RCP<BTLOC> container = Teuchos::rcp(new BTLOC);
   container->setMapsForBlocks(blockMaps);

   return container;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<LinearObjContainer>
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildGhostedLinearObjContainer() const
{
   std::vector<Teuchos::RCP<const MapType> > blockMaps;
   std::size_t blockDim = gidProviders_.size();
   for(std::size_t i=0;i<blockDim;i++)
      blockMaps.push_back(getGhostedMap(i));

   Teuchos::RCP<BTLOC> container = Teuchos::rcp(new BTLOC);
   container->setMapsForBlocks(blockMaps);

   return container;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
globalToGhostContainer(const LinearObjContainer & in,LinearObjContainer & out,int mem) const
{
   using Teuchos::is_null;

   typedef LinearObjContainer LOC;

   const BTLOC & b_in = Teuchos::dyn_cast<const BTLOC>(in);
   BTLOC & b_out = Teuchos::dyn_cast<BTLOC>(out);

   // Operations occur if the GLOBAL container has the correct targets!
   // Users set the GLOBAL continer arguments
   if ( !is_null(b_in.get_x()) && !is_null(b_out.get_x()) && ((mem & LOC::X)==LOC::X))
     globalToGhostThyraVector(b_in.get_x(),b_out.get_x());

   if ( !is_null(b_in.get_dxdt()) && !is_null(b_out.get_dxdt()) && ((mem & LOC::DxDt)==LOC::DxDt))
     globalToGhostThyraVector(b_in.get_dxdt(),b_out.get_dxdt());

   if ( !is_null(b_in.get_f()) && !is_null(b_out.get_f()) && ((mem & LOC::F)==LOC::F))
      globalToGhostThyraVector(b_in.get_f(),b_out.get_f());
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
ghostToGlobalContainer(const LinearObjContainer & in,LinearObjContainer & out,int mem) const
{
   using Teuchos::is_null;

   typedef LinearObjContainer LOC;

   const BTLOC & b_in = Teuchos::dyn_cast<const BTLOC>(in);
   BTLOC & b_out = Teuchos::dyn_cast<BTLOC>(out);

   // Operations occur if the GLOBAL container has the correct targets!
   // Users set the GLOBAL continer arguments
   if ( !is_null(b_in.get_x()) && !is_null(b_out.get_x()) && ((mem & LOC::X)==LOC::X))
     ghostToGlobalThyraVector(b_in.get_x(),b_out.get_x());

   if ( !is_null(b_in.get_f()) && !is_null(b_out.get_f()) && ((mem & LOC::F)==LOC::F))
     ghostToGlobalThyraVector(b_in.get_f(),b_out.get_f());

   if ( !is_null(b_in.get_A()) && !is_null(b_out.get_A()) && ((mem & LOC::Mat)==LOC::Mat))
     ghostToGlobalThyraMatrix(*b_in.get_A(),*b_out.get_A());
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
adjustForDirichletConditions(const LinearObjContainer & localBCRows,
                             const LinearObjContainer & globalBCRows,
                             LinearObjContainer & ghostedObjs,
                             bool zeroVectorRows, bool adjustX) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::LinearOpBase;
   using Thyra::PhysicallyBlockedLinearOpBase;
   using Thyra::VectorBase;
   using Thyra::ProductVectorBase;

   std::size_t blockDim = gidProviders_.size();

   // first cast to block LOCs
   const BTLOC & b_localBCRows = Teuchos::dyn_cast<const BTLOC>(localBCRows);
   const BTLOC & b_globalBCRows = Teuchos::dyn_cast<const BTLOC>(globalBCRows);
   BTLOC & b_ghosted = Teuchos::dyn_cast<BTLOC>(ghostedObjs);

   TEUCHOS_ASSERT(b_localBCRows.get_f()!=Teuchos::null);
   TEUCHOS_ASSERT(b_globalBCRows.get_f()!=Teuchos::null);

   // cast each component as needed to their product form
   RCP<PhysicallyBlockedLinearOpBase<ScalarT> > A = rcp_dynamic_cast<PhysicallyBlockedLinearOpBase<ScalarT> >(b_ghosted.get_A());
   RCP<ProductVectorBase<ScalarT> > f = rcp_dynamic_cast<ProductVectorBase<ScalarT> >(b_ghosted.get_f());
   RCP<ProductVectorBase<ScalarT> > local_bcs  = rcp_dynamic_cast<ProductVectorBase<ScalarT> >(b_localBCRows.get_f(),true);
   RCP<ProductVectorBase<ScalarT> > global_bcs = rcp_dynamic_cast<ProductVectorBase<ScalarT> >(b_globalBCRows.get_f(),true);

   if(adjustX) f = rcp_dynamic_cast<ProductVectorBase<ScalarT> >(b_ghosted.get_x());

   // sanity check!
   if(A!=Teuchos::null) TEUCHOS_ASSERT(A->productRange()->numBlocks()==(int) blockDim);
   if(A!=Teuchos::null) TEUCHOS_ASSERT(A->productDomain()->numBlocks()==(int) blockDim);
   if(f!=Teuchos::null) TEUCHOS_ASSERT(f->productSpace()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(local_bcs->productSpace()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(global_bcs->productSpace()->numBlocks()==(int) blockDim);

   for(std::size_t i=0;i<blockDim;i++) {
      // grab epetra vector
      RCP<const VectorType> t_local_bcs = rcp_dynamic_cast<const ThyraVector>(local_bcs->getVectorBlock(i),true)->getConstTpetraVector();
      RCP<const VectorType> t_global_bcs = rcp_dynamic_cast<const ThyraVector>(global_bcs->getVectorBlock(i),true)->getConstTpetraVector();

      // pull out epetra values
      RCP<VectorBase<ScalarT> > th_f = (f==Teuchos::null) ? Teuchos::null : f->getNonconstVectorBlock(i);
      RCP<VectorType> t_f;
      if(th_f==Teuchos::null)
        t_f = Teuchos::null;
      else
        t_f = rcp_dynamic_cast<ThyraVector>(th_f,true)->getTpetraVector();

      for(std::size_t j=0;j<blockDim;j++) {
        RCP<const MapType> map_i = getGhostedMap(i);
        RCP<const MapType> map_j = getGhostedMap(j);

         // pull out epetra values
         RCP<LinearOpBase<ScalarT> > th_A = (A== Teuchos::null)? Teuchos::null : A->getNonconstBlock(i,j);

         // don't do anyting if opertor is null
         RCP<CrsMatrixType> t_A;
         if(th_A==Teuchos::null)
            t_A = Teuchos::null;
         else {
            RCP<OperatorType> t_A_op = rcp_dynamic_cast<ThyraLinearOp>(th_A,true)->getTpetraOperator();
            t_A = rcp_dynamic_cast<CrsMatrixType>(t_A_op,true);
         }

         // adjust Block operator
         if(t_A!=Teuchos::null) {
           t_A->resumeFill();
         }

         adjustForDirichletConditions(*t_local_bcs,*t_global_bcs,t_f.ptr(),t_A.ptr(),zeroVectorRows);

         if(t_A!=Teuchos::null) {
           //t_A->fillComplete(map_j,map_i);
         }

         t_f = Teuchos::null; // this is so we only adjust it once on the first pass
      }
   }
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
adjustForDirichletConditions(const VectorType & local_bcs,
                             const VectorType & global_bcs,
                             const Teuchos::Ptr<VectorType> & f,
                             const Teuchos::Ptr<CrsMatrixType> & A,
                             bool zeroVectorRows) const
{
   if(f==Teuchos::null && A==Teuchos::null)
      return;

   Teuchos::ArrayRCP<ScalarT> f_array = f!=Teuchos::null ? f->get1dViewNonConst() : Teuchos::null;

   Teuchos::ArrayRCP<const ScalarT> local_bcs_array = local_bcs.get1dView();
   Teuchos::ArrayRCP<const ScalarT> global_bcs_array = global_bcs.get1dView();

   TEUCHOS_ASSERT(local_bcs.getLocalLength()==global_bcs.getLocalLength());
   for(std::size_t i=0;i<local_bcs.getLocalLength();i++) {
      if(global_bcs_array[i]==0.0)
         continue;

      if(local_bcs_array[i]==0.0 || zeroVectorRows) {
         // this boundary condition was NOT set by this processor

         // if they exist put 0.0 in each entry
         if(!Teuchos::is_null(f))
            f_array[i] = 0.0;
         if(!Teuchos::is_null(A)) {
            std::size_t numEntries = 0;
            std::size_t sz = A->getNumEntriesInLocalRow(i);
	    typename CrsMatrixType::nonconst_local_inds_host_view_type indices("indices", sz);
	    typename CrsMatrixType::nonconst_values_host_view_type values("values", sz);

            A->getLocalRowCopy(i,indices,values,numEntries);

            for(std::size_t c=0;c<numEntries;c++)
	      values(c) = 0.0;

            A->replaceLocalValues(i,indices,values);
         }
      }
      else {
         // this boundary condition was set by this processor

         ScalarT scaleFactor = global_bcs_array[i];

         // if they exist scale linear objects by scale factor
         if(!Teuchos::is_null(f))
            f_array[i] /= scaleFactor;
         if(!Teuchos::is_null(A)) {
            std::size_t numEntries = 0;
            std::size_t sz = A->getNumEntriesInLocalRow(i);
	    typename CrsMatrixType::nonconst_local_inds_host_view_type indices("indices", sz);
	    typename CrsMatrixType::nonconst_values_host_view_type values("values", sz);

            A->getLocalRowCopy(i,indices,values,numEntries);

            for(std::size_t c=0;c<numEntries;c++)
	      values(c) /= scaleFactor;

            A->replaceLocalValues(i,indices,values);
         }
      }
   }
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
applyDirichletBCs(const LinearObjContainer & /* counter */,
                  LinearObjContainer & /* result */) const
{
  TEUCHOS_ASSERT(false); // not yet implemented
}

///////////////////////////////////////////////////////////////////////////////
//
//  buildReadOnlyDomainContainer()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename ScalarT, typename LocalOrdinalT,
  typename GlobalOrdinalT, typename NodeT>
Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData>
BlockedTpetraLinearObjFactory<Traits, ScalarT, LocalOrdinalT, GlobalOrdinalT,
  NodeT>::
buildReadOnlyDomainContainer() const
{
  using std::vector;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using BVROGED = panzer::BlockedVector_ReadOnly_GlobalEvaluationData;
  using TVROGED = panzer::TpetraVector_ReadOnly_GlobalEvaluationData<ScalarT,
    LocalOrdinalT, GlobalOrdinalT, NodeT>;
  vector<RCP<ReadOnlyVector_GlobalEvaluationData>> gedBlocks;
  for (int i(0); i < getBlockColCount(); ++i)
  {
    auto tvroged = rcp(new TVROGED);
    tvroged->initialize(getGhostedImport(i), getGhostedMap(i), getMap(i));
    gedBlocks.push_back(tvroged);
  }
  auto ged = rcp(new BVROGED);
  ged->initialize(getGhostedThyraDomainSpace(), getThyraDomainSpace(),
    gedBlocks);
  return ged;
} // end of buildReadOnlyDomainContainer()

#ifdef PANZER_HAVE_EPETRA_STACK
///////////////////////////////////////////////////////////////////////////////
//
//  buildWriteDomainContainer()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename ScalarT, typename LocalOrdinalT,
  typename GlobalOrdinalT, typename NodeT>
Teuchos::RCP<WriteVector_GlobalEvaluationData>
BlockedTpetraLinearObjFactory<Traits, ScalarT, LocalOrdinalT, GlobalOrdinalT,
  NodeT>::
buildWriteDomainContainer() const
{
  using std::logic_error;
  using Teuchos::rcp;
  using EVWGED = panzer::EpetraVector_Write_GlobalEvaluationData;
  auto ged = rcp(new EVWGED);
  TEUCHOS_TEST_FOR_EXCEPTION(true, logic_error, "NOT YET IMPLEMENTED")
  return ged;
} // end of buildWriteDomainContainer()
#endif // PANZER_HAVE_EPETRA_STACK

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::MpiComm<int> BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getComm() const
{
   return *comm_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
initializeContainer(int mem,LinearObjContainer & loc) const
{
   BTLOC & bloc = Teuchos::dyn_cast<BTLOC>(loc);
   initializeContainer(mem,bloc);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
initializeGhostedContainer(int mem,LinearObjContainer & loc) const
{
   BTLOC & bloc = Teuchos::dyn_cast<BTLOC>(loc);
   initializeGhostedContainer(mem,bloc);
}

// Generic methods
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
initializeContainer(int mem,BTLOC & loc) const
{
   typedef LinearObjContainer LOC;

   loc.clear();

   if((mem & LOC::X) == LOC::X)
      loc.set_x(getThyraDomainVector());

   if((mem & LOC::DxDt) == LOC::DxDt)
      loc.set_dxdt(getThyraDomainVector());

   if((mem & LOC::F) == LOC::F)
      loc.set_f(getThyraRangeVector());

   if((mem & LOC::Mat) == LOC::Mat)
      loc.set_A(getThyraMatrix());
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
initializeGhostedContainer(int mem,BTLOC & loc) const
{
   typedef LinearObjContainer LOC;

   loc.clear();

   if((mem & LOC::X) == LOC::X)
      loc.set_x(getGhostedThyraDomainVector());

   if((mem & LOC::DxDt) == LOC::DxDt)
      loc.set_dxdt(getGhostedThyraDomainVector());

   if((mem & LOC::F) == LOC::F) {
      loc.set_f(getGhostedThyraRangeVector());
      loc.setRequiresDirichletAdjustment(true);
   }

   if((mem & LOC::Mat) == LOC::Mat) {
      loc.set_A(getGhostedThyraMatrix());
      loc.setRequiresDirichletAdjustment(true);
   }
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
addExcludedPair(int rowBlock,int colBlock)
{
   excludedPairs_.insert(std::make_pair(rowBlock,colBlock));
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
addExcludedPairs(const std::vector<std::pair<int,int> > & exPairs)
{
   for(std::size_t i=0;i<exPairs.size();i++)
      excludedPairs_.insert(exPairs[i]);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const GlobalIndexer>
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGlobalIndexer(int i) const
{
   return gidProviders_[i];
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
makeRoomForBlocks(std::size_t blockCnt)
{
   maps_.resize(blockCnt);
   ghostedMaps_.resize(blockCnt);
   importers_.resize(blockCnt);
   exporters_.resize(blockCnt);
}

// Thyra methods
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getThyraDomainSpace() const
{
   if(domainSpace_==Teuchos::null) {
      // loop over all vectors and build the vector space
      std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > > vsArray;
      for(std::size_t i=0;i<gidProviders_.size();i++)
         vsArray.push_back(Thyra::createVectorSpace<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getMap(i)));

      domainSpace_ = Thyra::productVectorSpace<ScalarT>(vsArray);
   }

   return domainSpace_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getThyraRangeSpace() const
{
   if(rangeSpace_==Teuchos::null) {
      // loop over all vectors and build the vector space
      std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > > vsArray;
      for(std::size_t i=0;i<gidProviders_.size();i++)
         vsArray.push_back(Thyra::createVectorSpace<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getMap(i)));

      rangeSpace_ = Thyra::productVectorSpace<ScalarT>(vsArray);
   }

   return rangeSpace_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getThyraDomainSpace(int blk) const
{
   if(domainSpace_==Teuchos::null) {
     getThyraDomainSpace();
   }

   return domainSpace_->getBlock(blk);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getThyraRangeSpace(int blk) const
{
   if(rangeSpace_==Teuchos::null) {
     getThyraRangeSpace();
   }

   return rangeSpace_->getBlock(blk);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Thyra::VectorBase<ScalarT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getThyraDomainVector() const
{
   Teuchos::RCP<Thyra::VectorBase<ScalarT> > vec =
      Thyra::createMember<ScalarT>(*getThyraDomainSpace());
   Thyra::assign(vec.ptr(),0.0);

   Teuchos::RCP<Thyra::ProductVectorBase<ScalarT> > p_vec = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<ScalarT> >(vec);
   for(std::size_t i=0;i<gidProviders_.size();i++) {
      TEUCHOS_ASSERT(Teuchos::rcp_dynamic_cast<Thyra::SpmdVectorBase<ScalarT> >(p_vec->getNonconstVectorBlock(i))->spmdSpace()->localSubDim()==
                     Teuchos::as<int>(getMap(i)->getLocalNumElements()));
   }

   return vec;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Thyra::VectorBase<ScalarT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getThyraRangeVector() const
{
   Teuchos::RCP<Thyra::VectorBase<ScalarT> > vec =
      Thyra::createMember<ScalarT>(*getThyraRangeSpace());
   Thyra::assign(vec.ptr(),0.0);

   return vec;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Thyra::LinearOpBase<ScalarT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getThyraMatrix() const
{
   Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<ScalarT> > blockedOp = Thyra::defaultBlockedLinearOp<ScalarT>();

   // get the block dimension
   std::size_t blockDim = gidProviders_.size();

   // this operator will be square
   blockedOp->beginBlockFill(blockDim,blockDim);

   // loop over each block
   for(std::size_t i=0;i<blockDim;i++) {
      for(std::size_t j=0;j<blockDim;j++) {
         if(excludedPairs_.find(std::make_pair(i,j))==excludedPairs_.end()) {
            // build (i,j) block matrix and add it to blocked operator
            Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > block = Thyra::createLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getTpetraMatrix(i,j));
            blockedOp->setNonconstBlock(i,j,block);
         }
      }
   }

   // all done
   blockedOp->endBlockFill();

   return blockedOp;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedThyraDomainSpace() const
{
   if(ghostedDomainSpace_==Teuchos::null) {
      // loop over all vectors and build the vector space
      std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > > vsArray;
      for(std::size_t i=0;i<gidProviders_.size();i++)
         vsArray.push_back(Thyra::createVectorSpace<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getGhostedMap(i)));

      ghostedDomainSpace_ = Thyra::productVectorSpace<ScalarT>(vsArray);
   }

   return ghostedDomainSpace_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedThyraRangeSpace() const
{
   if(ghostedRangeSpace_==Teuchos::null) {
      // loop over all vectors and build the vector space
      std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<ScalarT> > > vsArray;
      for(std::size_t i=0;i<gidProviders_.size();i++)
         vsArray.push_back(Thyra::createVectorSpace<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getGhostedMap(i)));

      ghostedRangeSpace_ = Thyra::productVectorSpace<ScalarT>(vsArray);
   }

   return ghostedRangeSpace_;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Thyra::VectorBase<ScalarT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedThyraDomainVector() const
{
   Teuchos::RCP<Thyra::VectorBase<ScalarT> > vec =
      Thyra::createMember<ScalarT>(*getGhostedThyraDomainSpace());
   Thyra::assign(vec.ptr(),0.0);

   return vec;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Thyra::VectorBase<ScalarT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedThyraRangeVector() const
{
   Teuchos::RCP<Thyra::VectorBase<ScalarT> > vec =
      Thyra::createMember<ScalarT>(*getGhostedThyraRangeSpace());
   Thyra::assign(vec.ptr(),0.0);

   return vec;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Thyra::BlockedLinearOpBase<ScalarT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedThyraMatrix() const
{
   Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<ScalarT> > blockedOp = Thyra::defaultBlockedLinearOp<ScalarT>();

   // get the block dimension
   std::size_t blockDim = gidProviders_.size();

   // this operator will be square
   blockedOp->beginBlockFill(blockDim,blockDim);

   // loop over each block
   for(std::size_t i=0;i<blockDim;i++) {
      for(std::size_t j=0;j<blockDim;j++) {
         if(excludedPairs_.find(std::make_pair(i,j))==excludedPairs_.end()) {
            // build (i,j) block matrix and add it to blocked operator
            Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > block = Thyra::createLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getGhostedTpetraMatrix(i,j));
            blockedOp->setNonconstBlock(i,j,block);
         }
      }
   }

   // all done
   blockedOp->endBlockFill();

   return blockedOp;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
ghostToGlobalThyraVector(const Teuchos::RCP<const Thyra::VectorBase<ScalarT> > & in,
                         const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & out) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;

   std::size_t blockDim = gidProviders_.size();

   // get product vectors
   RCP<const ProductVectorBase<ScalarT> > prod_in = rcp_dynamic_cast<const ProductVectorBase<ScalarT> >(in,true);
   RCP<ProductVectorBase<ScalarT> > prod_out      = rcp_dynamic_cast<ProductVectorBase<ScalarT> >(out,true);

   TEUCHOS_ASSERT(prod_in->productSpace()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(prod_out->productSpace()->numBlocks()==(int) blockDim);

   for(std::size_t i=0;i<blockDim;i++) {
      // first get each Tpetra vector
      RCP<const VectorType> tp_in = rcp_dynamic_cast<const ThyraVector>(prod_in->getVectorBlock(i),true)->getConstTpetraVector();
      RCP<VectorType> tp_out      = rcp_dynamic_cast<ThyraVector>(prod_out->getNonconstVectorBlock(i),true)->getTpetraVector();

      // use Tpetra to do global communication
      ghostToGlobalTpetraVector(i,*tp_in,*tp_out);
   }
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
ghostToGlobalThyraMatrix(const Thyra::LinearOpBase<ScalarT> & in,Thyra::LinearOpBase<ScalarT> & out) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Teuchos::dyn_cast;
   using Thyra::LinearOpBase;
   using Thyra::PhysicallyBlockedLinearOpBase;

   std::size_t blockDim = gidProviders_.size();

   // get product vectors
   const PhysicallyBlockedLinearOpBase<ScalarT> & prod_in = dyn_cast<const PhysicallyBlockedLinearOpBase<ScalarT> >(in);
   PhysicallyBlockedLinearOpBase<ScalarT> & prod_out      = dyn_cast<PhysicallyBlockedLinearOpBase<ScalarT> >(out);

   TEUCHOS_ASSERT(prod_in.productRange()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(prod_in.productDomain()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(prod_out.productRange()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(prod_out.productDomain()->numBlocks()==(int) blockDim);

   for(std::size_t i=0;i<blockDim;i++) {
      for(std::size_t j=0;j<blockDim;j++) {
         if(excludedPairs_.find(std::make_pair(i,j))==excludedPairs_.end()) {
            // extract the blocks
            RCP<const LinearOpBase<ScalarT> > th_in = prod_in.getBlock(i,j);
            RCP<LinearOpBase<ScalarT> > th_out = prod_out.getNonconstBlock(i,j);

            // sanity check
            TEUCHOS_ASSERT(th_in!=Teuchos::null);
            TEUCHOS_ASSERT(th_out!=Teuchos::null);

            // get the epetra version of the blocks
            RCP<const OperatorType> tp_op_in = rcp_dynamic_cast<const ThyraLinearOp>(th_in,true)->getConstTpetraOperator();
            RCP<OperatorType> tp_op_out      = rcp_dynamic_cast<ThyraLinearOp>(th_out,true)->getTpetraOperator();

            RCP<const CrsMatrixType> tp_in = rcp_dynamic_cast<const CrsMatrixType>(tp_op_in,true);
            RCP<CrsMatrixType> tp_out      = rcp_dynamic_cast<CrsMatrixType>(tp_op_out,true);

            // use Tpetra to do global communication
            ghostToGlobalTpetraMatrix(i,*tp_in,*tp_out);
         }
      }
   }
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
globalToGhostThyraVector(const Teuchos::RCP<const Thyra::VectorBase<ScalarT> > & in,
                         const Teuchos::RCP<Thyra::VectorBase<ScalarT> > & out) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;

   std::size_t blockDim = gidProviders_.size();

   // get product vectors
   RCP<const ProductVectorBase<ScalarT> > prod_in = rcp_dynamic_cast<const ProductVectorBase<ScalarT> >(in,true);
   RCP<ProductVectorBase<ScalarT> > prod_out      = rcp_dynamic_cast<ProductVectorBase<ScalarT> >(out,true);

   TEUCHOS_ASSERT(prod_in->productSpace()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(prod_out->productSpace()->numBlocks()==(int) blockDim);

   for(std::size_t i=0;i<blockDim;i++) {
      // first get each Tpetra vector
      RCP<const VectorType> tp_in = rcp_dynamic_cast<const ThyraVector>(prod_in->getVectorBlock(i),true)->getConstTpetraVector();
      RCP<VectorType> tp_out      = rcp_dynamic_cast<ThyraVector>(prod_out->getNonconstVectorBlock(i),true)->getTpetraVector();

      // use Tpetra to do global communication
      globalToGhostTpetraVector(i,*tp_in,*tp_out);
   }
}

// Tpetra methods
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
ghostToGlobalTpetraVector(int i,const VectorType & in,VectorType & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<const ExportType> exporter = getGhostedExport(i);
   out.putScalar(0.0);
   out.doExport(in,*exporter,Tpetra::ADD);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
ghostToGlobalTpetraMatrix(int blockRow,const CrsMatrixType & in,CrsMatrixType & out) const
{
   using Teuchos::RCP;

   RCP<const MapType> map_i = out.getRangeMap();
   RCP<const MapType> map_j = out.getDomainMap();

   // do the global distribution
   RCP<const ExportType> exporter = getGhostedExport(blockRow);

   out.resumeFill();
   out.setAllToScalar(0.0);
   out.doExport(in,*exporter,Tpetra::ADD);
   out.fillComplete(map_j,map_i);
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
globalToGhostTpetraVector(int i,const VectorType & in,VectorType & out) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<const ImportType> importer = getGhostedImport(i);
   out.putScalar(0.0);
   out.doImport(in,*importer,Tpetra::INSERT);
}

// get the map from the matrix
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getMap(int i) const
{
   if(maps_[i]==Teuchos::null)
      maps_[i] = buildTpetraMap(i);

   return maps_[i];
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedMap(int i) const
{
   if(ghostedMaps_[i]==Teuchos::null)
      ghostedMaps_[i] = buildTpetraGhostedMap(i);

   return ghostedMaps_[i];
}

// get the graph of the crs matrix
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGraph(int i,int j) const
{
  typedef std::unordered_map<std::pair<int,int>,Teuchos::RCP<const CrsGraphType>,panzer::pair_hash> GraphMap;

   typename GraphMap::const_iterator itr = graphs_.find(std::make_pair(i,j));
   Teuchos::RCP<const CrsGraphType> graph;
   if(itr==graphs_.end()) {
      graph = buildTpetraGraph(i,j);
      graphs_[std::make_pair(i,j)] = graph;
   }
   else
      graph = itr->second;

   TEUCHOS_ASSERT(graph!=Teuchos::null);
   return graph;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedGraph(int i,int j) const
{
  typedef std::unordered_map<std::pair<int,int>,Teuchos::RCP<const CrsGraphType>,panzer::pair_hash> GraphMap;

   typename GraphMap::const_iterator itr = ghostedGraphs_.find(std::make_pair(i,j));
   Teuchos::RCP<const CrsGraphType> ghostedGraph;
   if(itr==ghostedGraphs_.end()) {
      ghostedGraph = buildTpetraGhostedGraph(i,j);
      ghostedGraphs_[std::make_pair(i,j)] = ghostedGraph;
   }
   else
      ghostedGraph = itr->second;

   TEUCHOS_ASSERT(ghostedGraph!=Teuchos::null);
   return ghostedGraph;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const  Tpetra::Import<LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedImport(int i) const
{
   if(importers_[i]==Teuchos::null)
      importers_[i] = Teuchos::rcp(new ImportType(getMap(i),getGhostedMap(i)));

   return importers_[i];
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const  Tpetra::Export<LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedExport(int i) const
{
   if(exporters_[i]==Teuchos::null)
      exporters_[i] = Teuchos::rcp(new ExportType(getGhostedMap(i),getMap(i)));

   return exporters_[i];
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildTpetraMap(int i) const
{
   std::vector<GlobalOrdinalT> indices;

   // get the global indices
   getGlobalIndexer(i)->getOwnedIndices(indices);

   return Teuchos::rcp(new MapType(Teuchos::OrdinalTraits<GlobalOrdinalT>::invalid(),indices,0,comm_));
}

// build the ghosted map
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Tpetra::Map<LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildTpetraGhostedMap(int i) const
{
   std::vector<GlobalOrdinalT> indices;

   // get the global indices
   getGlobalIndexer(i)->getOwnedAndGhostedIndices(indices);

   return Teuchos::rcp(new MapType(Teuchos::OrdinalTraits<GlobalOrdinalT>::invalid(),indices,0,comm_));
}

// get the graph of the crs matrix
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildTpetraGraph(int i,int j) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // build the map and allocate the space for the graph and
   // grab the ghosted graph
   RCP<const MapType> map_i = getMap(i);
   RCP<const MapType> map_j = getMap(j);

   RCP<CrsGraphType> graph  = rcp(new CrsGraphType(map_i,0));
   RCP<const CrsGraphType> oGraph = getGhostedGraph(i,j);

   // perform the communication to finish building graph
   RCP<const ExportType> exporter = getGhostedExport(i);
   graph->doExport( *oGraph, *exporter, Tpetra::INSERT );
   graph->fillComplete(map_j,map_i);

   return graph;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildTpetraGhostedGraph(int i,int j) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // build the map and allocate the space for the graph and
   // grab the ghosted graph
   RCP<const MapType> map_i = getGhostedMap(i);
   RCP<const MapType> map_j = getGhostedMap(j);

   std::vector<std::string> elementBlockIds;

   Teuchos::RCP<const GlobalIndexer> rowProvider, colProvider;

   rowProvider = getGlobalIndexer(i);
   colProvider = getGlobalIndexer(j);

   gidProviders_[0]->getElementBlockIds(elementBlockIds); // each sub provider "should" have the
                                                          // same element blocks

   // Count number of entries in each row of graph; needed for graph constructor
   std::vector<size_t> nEntriesPerRow(map_i->getLocalNumElements(), 0);
   std::vector<std::string>::const_iterator blockItr;
   for(blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;
      // grab elements for this block
      const std::vector<LocalOrdinalT> & elements = gidProviders_[0]->getElementBlock(blockId); // each sub provider "should" have the
                                                                                                // same elements in each element block

      // get information about number of indicies
      std::vector<GlobalOrdinalT> row_gids;
      std::vector<GlobalOrdinalT> col_gids;

      // loop over the elemnts
      for(std::size_t elmt=0;elmt<elements.size();elmt++) {

         rowProvider->getElementGIDs(elements[elmt],row_gids);
         colProvider->getElementGIDs(elements[elmt],col_gids);
         for(std::size_t row=0;row<row_gids.size();row++) {
            LocalOrdinalT lid = map_i->getLocalElement(row_gids[row]);
            nEntriesPerRow[lid] += col_gids.size();
         }
      }
   }
   Teuchos::ArrayView<const size_t> nEntriesPerRowView(nEntriesPerRow);
   RCP<CrsGraphType> graph  = rcp(new CrsGraphType(map_i,map_j, nEntriesPerRowView));



   // graph information about the mesh
   for(blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;

      // grab elements for this block
      const std::vector<LocalOrdinalT> & elements = gidProviders_[0]->getElementBlock(blockId); // each sub provider "should" have the
                                                                                                // same elements in each element block

      // get information about number of indicies
      std::vector<GlobalOrdinalT> row_gids;
      std::vector<GlobalOrdinalT> col_gids;

      // loop over the elemnts
      for(std::size_t elmt=0;elmt<elements.size();elmt++) {

         rowProvider->getElementGIDs(elements[elmt],row_gids);
         colProvider->getElementGIDs(elements[elmt],col_gids);
         for(std::size_t row=0;row<row_gids.size();row++)
            graph->insertGlobalIndices(row_gids[row],col_gids);
      }
   }

   // finish filling the graph: Make sure the colmap and row maps coincide to
   //                           minimize calls to LID lookups
   graph->fillComplete(getMap(j),getMap(i));

   return graph;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getTpetraMatrix(int i,int j) const
{
   Teuchos::RCP<const MapType> map_i = getMap(i);
   Teuchos::RCP<const MapType> map_j = getMap(j);

   Teuchos::RCP<const CrsGraphType> tGraph = getGraph(i,j);
   Teuchos::RCP<CrsMatrixType> mat = Teuchos::rcp(new CrsMatrixType(tGraph));
   mat->fillComplete(map_j,map_i);

   return mat;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedTpetraMatrix(int i,int j) const
{
   Teuchos::RCP<const MapType> map_i = getGhostedMap(i);
   Teuchos::RCP<const MapType> map_j = getGhostedMap(j);

   Teuchos::RCP<const CrsGraphType> tGraph = getGhostedGraph(i,j);
   Teuchos::RCP<CrsMatrixType> mat = Teuchos::rcp(new CrsMatrixType(tGraph));
   mat->fillComplete(getMap(j),getMap(i));

   return mat;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getTpetraDomainVector(int i) const
{
   Teuchos::RCP<const MapType> tMap = getMap(i);
   return Teuchos::rcp(new VectorType(tMap));
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedTpetraDomainVector(int i) const
{
   Teuchos::RCP<const MapType> tMap = getGhostedMap(i);
   return Teuchos::rcp(new VectorType(tMap));
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getTpetraRangeVector(int i) const
{
   Teuchos::RCP<const MapType> tMap = getMap(i);
   return Teuchos::rcp(new VectorType(tMap));
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Tpetra::Vector<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getGhostedTpetraRangeVector(int i) const
{
   Teuchos::RCP<const MapType> tMap = getGhostedMap(i);
   return Teuchos::rcp(new VectorType(tMap));
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
int
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getBlockRowCount() const
{
   return gidProviders_.size();
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
int
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getBlockColCount() const
{
   return gidProviders_.size();
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
beginFill(LinearObjContainer & loc) const
{
  BTLOC & tloc = Teuchos::dyn_cast<BTLOC>(loc);
  if(tloc.get_A()!=Teuchos::null)
    tloc.beginFill();
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
endFill(LinearObjContainer & loc) const
{
  BTLOC & tloc = Teuchos::dyn_cast<BTLOC>(loc);
  if(tloc.get_A()!=Teuchos::null)
    tloc.endFill();
}

}

#endif // __Panzer_BlockedTpetraLinearObjFactory_impl_hpp__
