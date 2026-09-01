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
#include <utility>
#ifdef PANZER_HAVE_EPETRA_STACK
#include "Panzer_EpetraVector_Write_GlobalEvaluationData.hpp"                    // JMG:  Remove this eventually.
#endif
#include "Panzer_TpetraVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_GlobalIndexer.hpp"

#include "KokkosSparse_SortCrs.hpp"

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
                              const Teuchos::RCP<const BlockedDOFManager> & gidProvider,
                              bool useFEAssembly)
   : blockProvider_(gidProvider), blockedDOFManager_(gidProvider), comm_(comm)
   , useFEAssembly_(useFEAssembly)
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
                              const std::vector<Teuchos::RCP<const panzer::GlobalIndexer>> & gidProviders,
                              bool useFEAssembly)
  : gidProviders_(gidProviders), comm_(comm)
  , useFEAssembly_(useFEAssembly)
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
            //
            // In FE mode the spaces must be stated explicitly rather than deduced from the
            // matrix. A shared FECrsMatrix reports whichever maps its ACTIVE view has, and
            // that view flips between owned+shared and owned across the assembly cycle --
            // a freshly constructed one even starts in owned+shared. Deducing the spaces
            // would stamp this owned operator with ghosted spaces. The classic path is
            // unaffected, since its matrix only ever has the owned maps.
            Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > block;
            if(useFEAssembly_)
               block = Thyra::tpetraLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(
                          getThyraRangeSpace(i),getThyraDomainSpace(j),getTpetraMatrix(i,j));
            else
               block = Thyra::createLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getTpetraMatrix(i,j));
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
            //
            // See getThyraMatrix(): in FE mode the spaces are stated explicitly, since the
            // shared matrix's own maps follow its active view. This ghosted operator keeps
            // the ghosted spaces, matching the ghosted container's vectors and the classic
            // path, even though the matrix underneath is the same object the owned operator
            // wraps.
            Teuchos::RCP<Thyra::LinearOpBase<ScalarT> > block;
            if(useFEAssembly_) {
               Teuchos::RCP<const Thyra::ProductVectorSpaceBase<ScalarT> > gRange
                  = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<ScalarT> >(getGhostedThyraRangeSpace(),true);
               Teuchos::RCP<const Thyra::ProductVectorSpaceBase<ScalarT> > gDomain
                  = Teuchos::rcp_dynamic_cast<const Thyra::ProductVectorSpaceBase<ScalarT> >(getGhostedThyraDomainSpace(),true);
               block = Thyra::tpetraLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(
                          gRange->getBlock(i),gDomain->getBlock(j),getGhostedTpetraMatrix(i,j));
            }
            else
               block = Thyra::createLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(getGhostedTpetraMatrix(i,j));
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

            // In FE mode this block is one shared FECrsMatrix, and endFill()'s endAssembly()
            // has already migrated its ghost rows onto the owned rows. Exporting it onto
            // itself here would double every shared-interface contribution.
            //
            // The test has to be on the extracted Tpetra matrices: even when the blocks are
            // shared, the owned and ghosted containers hold distinct DefaultBlockedLinearOp
            // objects wrapping distinct Thyra::TpetraLinearOps, so comparing get_A() -- or
            // the Thyra blocks -- would never report a match.
            if(tp_in.get()==tp_out.get())
               continue;

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

template <class LocalOrdinalT>
struct entry_type {
  LocalOrdinalT row;
  LocalOrdinalT col;
};

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<const Tpetra::CrsGraph<LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildTpetraGhostedGraph(int i,int j) const
{
  PANZER_FUNC_TIME_MONITOR_DIFF("panzer::BlockedTpetraLinearObjFactory::buildTpetraGhostedGraph",BTLOF);

   using Teuchos::RCP;
   using Teuchos::rcp;
   using exec_space = typename CrsGraphType::execution_space;
   using memory_space = typename NodeT::memory_space;

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

   RCP<CrsGraphType> graph;
   if constexpr (NodeT::is_gpu) {

     // Gather elements from mesh blocks.
     size_t numElements;
     Kokkos::View<LocalOrdinalT *, memory_space> elementsFromBlocks;
     {
       auto numElementBlocks = elementBlockIds.size();

       std::vector<size_t> elementBlockOffsets(numElementBlocks + 1);
       elementBlockOffsets[0] = 0;

       numElements = 0;
       size_t blockNo = 0;
       std::vector<std::string>::const_iterator blockItr;
       for (blockItr = elementBlockIds.begin();
            blockItr != elementBlockIds.end(); ++blockItr) {
         std::string blockId = *blockItr;
         const std::vector<LocalOrdinalT> &elements =
             gidProviders_[0]->getElementBlock(
                 blockId); // each sub provider "should" have the
                           // same elements in each element block
         numElements += elements.size();
         ++blockNo;
         elementBlockOffsets[blockNo] = numElements;
       }
       elementsFromBlocks = Kokkos::View<LocalOrdinalT *, memory_space>(
           "elementsFromBlocks", numElements);
       blockNo = 0;
       for (blockItr = elementBlockIds.begin();
            blockItr != elementBlockIds.end(); ++blockItr) {
         std::string blockId = *blockItr;
         const std::vector<LocalOrdinalT> &elements =
             gidProviders_[0]->getElementBlock(
                 blockId); // each sub provider "should" have the
                           // same elements in each element block
         Kokkos::View<const LocalOrdinalT *, Kokkos::HostSpace,
                      Kokkos::MemoryTraits<Kokkos::Unmanaged>>
             elements_h(elements.data(), elements.size());
         Kokkos::deep_copy(
             Kokkos::subview(
                 elementsFromBlocks,
                 Kokkos::make_pair(elementBlockOffsets[blockNo],
                                   elementBlockOffsets[blockNo + 1])),
             elements_h);
         ++blockNo;
       }
     }

     {

       using local_graph_type = typename CrsGraphType::local_graph_device_type;
       using rowptr_type =
           typename local_graph_type::row_map_type::non_const_type;
       using colidx_type =
           typename local_graph_type::entries_type::non_const_type;

       using entries_map_type =
           Kokkos::UnorderedMap<entry_type<LocalOrdinalT>, void, exec_space>;

       auto numRows = map_i->getLocalNumElements();

       // We are overallocating by 1 here. This simplifies the logic below. But
       // we have to remember to take a subview in the end.
       rowptr_type rowptr("ghostedGraph_rowptr", numRows + 2);

       auto rowLIDs = rowProvider->getLIDs();
       auto colLIDs = colProvider->getLIDs();

       auto numDoFsPerElementRow = rowLIDs.extent(1);
       auto numDoFsPerElementCol = colLIDs.extent(1);

       auto capacity =
           numElements * numDoFsPerElementRow * numDoFsPerElementCol;
       entries_map_type entries(capacity);

       while (true) {

         // Loop over all elements and record the entries that we need in the
         // graph. Also start building the rowptr.
         Kokkos::parallel_for(
             "collect_entries", Kokkos::RangePolicy<exec_space>(0, numElements),
             KOKKOS_LAMBDA(const LocalOrdinalT k) {
               auto elementId = elementsFromBlocks(k);
               entry_type<LocalOrdinalT> entry;
               for (size_t dofNoRow = 0; dofNoRow < numDoFsPerElementRow;
                    ++dofNoRow) {
                 entry.row = rowLIDs(elementId, dofNoRow);
                 for (size_t dofNoCol = 0; dofNoCol < numDoFsPerElementCol;
                      ++dofNoCol) {
                   entry.col = colLIDs(elementId, dofNoCol);
                   auto result = entries.insert(entry);
                   if (result.success()) {
                     // New entry. We offset by 2 here.
                     Kokkos::atomic_inc(&rowptr(entry.row + 2));
                   }
                 }
               }
             });

         if (!entries.failed_insert()) {
           auto numEntries = entries.size();

           // Prefix sum to get offsets.
           // This is not the correct rowptr yet.
           // We have essentially shifted everything by one position.
           // This is useful for when we fill.
           typename rowptr_type::value_type numEntries2;
           Kokkos::parallel_scan(
               "prefix_sum", Kokkos::RangePolicy<exec_space>(0, numRows + 2),
               KOKKOS_LAMBDA(const size_t rlid,
                             typename rowptr_type::value_type &nnz,
                             const bool is_final) {
                 nnz += rowptr(rlid);
                 if (is_final)
                   rowptr(rlid) = nnz;
               },
               numEntries2);
           TEUCHOS_ASSERT_EQUALITY(numEntries, numEntries2);

           // The column indices.
           colidx_type colidx(
               Kokkos::ViewAllocateWithoutInitializing("ghostedGraph_colidx"),
               numEntries);

           // Fill the column indices.
           // We are using the rowptr to figure out offsets.
           // After this step the rowptr is correct.
           Kokkos::parallel_for(
               "fill", Kokkos::RangePolicy<exec_space>(0, entries.capacity()),
               KOKKOS_LAMBDA(const uint32_t c) {
                 if (entries.valid_at(c)) {
                   auto entry = entries.key_at(c);
                   auto offset =
                       Kokkos::atomic_fetch_inc(&rowptr(entry.row + 1));
                   colidx(offset) = entry.col;
                 }
               });

           // Sort the rows.
           KokkosSparse::sort_crs_graph(rowptr, colidx);

           // Create the graph
           graph = rcp(new CrsGraphType(
               map_i, map_j,
               Kokkos::subview(rowptr, Kokkos::make_pair((decltype(numRows))0,
                                                         numRows + 1)),
               colidx));
           graph->fillComplete(getMap(j), getMap(i));

           break;
         } else {
           // We ended up not having enough capacity in the UnorderedMap.
           // Bump it up and try again.
           std::cout << "Insufficient capacity: " << capacity << std::endl;
           capacity *= 2;
           Kokkos::deep_copy(rowptr, 0);
           entries = entries_map_type(capacity);
         }
       }
     }
   } else {
     // Count number of entries in each row of graph; needed for graph
     // constructor
     std::vector<size_t> nEntriesPerRow(map_i->getLocalNumElements(), 0);
     std::vector<std::string>::const_iterator blockItr;
     for (blockItr = elementBlockIds.begin(); blockItr != elementBlockIds.end();
          ++blockItr) {
       std::string blockId = *blockItr;
       // grab elements for this block
       const std::vector<LocalOrdinalT> &elements =
           gidProviders_[0]->getElementBlock(
               blockId); // each sub provider "should" have the
                         // same elements in each element block

       // get information about number of indicies
       std::vector<GlobalOrdinalT> row_gids;
       std::vector<GlobalOrdinalT> col_gids;

       // loop over the elemnts
       for (std::size_t elmt = 0; elmt < elements.size(); elmt++) {

         rowProvider->getElementGIDs(elements[elmt], row_gids);
         colProvider->getElementGIDs(elements[elmt], col_gids);
         for (std::size_t row = 0; row < row_gids.size(); row++) {
           LocalOrdinalT lid = map_i->getLocalElement(row_gids[row]);
           nEntriesPerRow[lid] += col_gids.size();
         }
       }
     }
     Teuchos::ArrayView<const size_t> nEntriesPerRowView(nEntriesPerRow);
     graph = rcp(new CrsGraphType(map_i, map_j, nEntriesPerRowView));

     // graph information about the mesh
     for (blockItr = elementBlockIds.begin(); blockItr != elementBlockIds.end();
          ++blockItr) {
       std::string blockId = *blockItr;

       // grab elements for this block
       const std::vector<LocalOrdinalT> &elements =
           gidProviders_[0]->getElementBlock(
               blockId); // each sub provider "should" have the
                         // same elements in each element block

       // get information about number of indicies
       std::vector<GlobalOrdinalT> row_gids;
       std::vector<GlobalOrdinalT> col_gids;

       // loop over the elemnts
       for (std::size_t elmt = 0; elmt < elements.size(); elmt++) {

         rowProvider->getElementGIDs(elements[elmt], row_gids);
         colProvider->getElementGIDs(elements[elmt], col_gids);
         for (std::size_t row = 0; row < row_gids.size(); row++)
           graph->insertGlobalIndices(row_gids[row], col_gids);
       }
     }

     // finish filling the graph: Make sure the colmap and row maps coincide to
     //                           minimize calls to LID lookups
     graph->fillComplete(getMap(j), getMap(i));
   }

   return graph;
}

// Build the FE graph for block (i,j).
//
// This walks the same element/GID traversal as buildTpetraGhostedGraph(i,j) -- rows come
// from indexer i, columns from indexer j -- but hands the result to the FECrsGraph "V2"
// constructor, which carries both the owned and the owned+shared graph in one object.
//
// There is no getColMap() in this class and none is needed: block (i,j) takes its rows from
// indexer i and its columns from indexer j, so block j's row map IS this block's column map.
// That is the same convention the classic path already uses -- buildTpetraGhostedGraph()
// fill-completes with (getMap(j),getMap(i)).
//
// The V2 constructor requires that the owned row/domain gids appear, in the same order, as a
// leading prefix of the owned+shared row/domain map. Every concrete panzer::GlobalIndexer
// builds getOwnedAndGhostedIndices() as owned_ followed by ghosted_, so getMap(i)/
// getGhostedMap(i) satisfy this directly.
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<typename BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::FECrsGraphType>
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
buildFEGraph(int i,int j) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // NOTE: these must be RCP<const MapType> (not RCP<MapType>). FECrsGraph has both a "V1"
   // ctor (4th positional arg = importer) and a "V2" ctor (4th positional arg =
   // ownedPlusSharedDomainMap); Teuchos::RCP's converting constructor is an unconstrained
   // template, so an RCP<MapType> is equally "convertible" to either as far as overload
   // resolution is concerned. Typing the locals as RCP<const MapType> binds them to V2's
   // domain-map parameter with no conversion at all, resolving the ambiguity in V2's favor.
   RCP<const MapType> ownedRowMap             = getMap(i);
   RCP<const MapType> ownedPlusSharedRowMap   = getGhostedMap(i);
   RCP<const MapType> ownedDomainMap          = getMap(j);
   RCP<const MapType> ownedPlusSharedDomainMap = getGhostedMap(j);

   RCP<const GlobalIndexer> rowProvider = getGlobalIndexer(i);
   RCP<const GlobalIndexer> colProvider = getGlobalIndexer(j);

   std::vector<std::string> elementBlockIds;
   gidProviders_[0]->getElementBlockIds(elementBlockIds); // each sub provider "should" have
                                                          // the same element blocks

   // count entries per owned+shared row, exactly as buildTpetraGhostedGraph() does
   std::vector<size_t> nEntriesPerRow(ownedPlusSharedRowMap->getLocalNumElements(),0);

   std::vector<std::string>::const_iterator blockItr;
   for(blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      const std::vector<LocalOrdinalT> & elements = gidProviders_[0]->getElementBlock(*blockItr);

      std::vector<GlobalOrdinalT> row_gids;
      std::vector<GlobalOrdinalT> col_gids;

      for(std::size_t elmt=0;elmt<elements.size();elmt++) {
         rowProvider->getElementGIDs(elements[elmt],row_gids);
         colProvider->getElementGIDs(elements[elmt],col_gids);
         for(std::size_t row=0;row<row_gids.size();row++) {
            LocalOrdinalT lid = ownedPlusSharedRowMap->getLocalElement(row_gids[row]);
            nEntriesPerRow[lid] += col_gids.size();
         }
      }
   }

   size_t maxNumRowEntries = 0;
   for(std::size_t r=0;r<nEntriesPerRow.size();r++)
      maxNumRowEntries = std::max(maxNumRowEntries,nEntriesPerRow[r]);

   RCP<FECrsGraphType> feGraph = rcp(new FECrsGraphType(
       ownedRowMap, ownedPlusSharedRowMap, maxNumRowEntries,
       ownedPlusSharedDomainMap,
       Teuchos::null,
       ownedDomainMap));

   // Panzer's DOFManager does not guarantee a locally owned element has an owned dof, so
   // Tpetra's debug-only check for that is too strict here; the cost is at most a structurally
   // empty column. Must be set after construction -- the ctor's validator rejects the option.
   {
      Teuchos::RCP<Teuchos::ParameterList> feGraphParams = Teuchos::parameterList();
      feGraphParams->set("Check Col GIDs In At Least One Owned Row",false);
      feGraph->setParameterList(feGraphParams);
   }

   feGraph->beginAssembly();
   for(blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      const std::vector<LocalOrdinalT> & elements = gidProviders_[0]->getElementBlock(*blockItr);

      std::vector<GlobalOrdinalT> row_gids;
      std::vector<GlobalOrdinalT> col_gids;

      for(std::size_t elmt=0;elmt<elements.size();elmt++) {
         rowProvider->getElementGIDs(elements[elmt],row_gids);
         colProvider->getElementGIDs(elements[elmt],col_gids);
         for(std::size_t row=0;row<row_gids.size();row++)
            feGraph->insertGlobalIndices(row_gids[row],col_gids);
      }
   }
   feGraph->endAssembly();

   return feGraph;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<typename BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::FECrsGraphType>
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getFEGraph(int i,int j) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(!useFEAssembly_,std::logic_error,
      "BlockedTpetraLinearObjFactory::getFEGraph: This factory was not constructed with "
      "FE assembly enabled.");

   typedef std::unordered_map<std::pair<int,int>,Teuchos::RCP<FECrsGraphType>,panzer::pair_hash> FEGraphMap;

   typename FEGraphMap::const_iterator itr = feGraphs_.find(std::make_pair(i,j));
   if(itr!=feGraphs_.end())
      return itr->second;

   Teuchos::RCP<FECrsGraphType> graph = buildFEGraph(i,j);
   feGraphs_[std::make_pair(i,j)] = graph;
   return graph;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<typename BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::FECrsMatrixType>
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getFEMatrix(int i,int j) const
{
   TEUCHOS_TEST_FOR_EXCEPTION(!useFEAssembly_,std::logic_error,
      "BlockedTpetraLinearObjFactory::getFEMatrix: This factory was not constructed with "
      "FE assembly enabled.");

   typedef std::unordered_map<std::pair<int,int>,Teuchos::RCP<FECrsMatrixType>,panzer::pair_hash> FEMatrixMap;

   typename FEMatrixMap::const_iterator itr = feMatrices_.find(std::make_pair(i,j));
   if(itr!=feMatrices_.end())
      return itr->second;

   // Cached deliberately. The owned and ghosted containers must end up wrapping the SAME
   // matrix per block, so that endAssembly() can migrate ghost contributions to the owned
   // rows in place. Handing back a fresh matrix here (as the classic getTpetraMatrix does)
   // would give the two containers disjoint storage and silently drop the ghost data.
   Teuchos::RCP<FECrsMatrixType> mat = Teuchos::rcp(new FECrsMatrixType(getFEGraph(i,j)));
   feMatrices_[std::make_pair(i,j)] = mat;
   return mat;
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >
BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
getTpetraMatrix(int i,int j) const
{
   // In FE mode the owned and ghosted views are the same object, so both getters return the
   // shared FE matrix. It is already fill-complete over (getMap(j),getMap(i)) once
   // endAssembly() has run, which is what callers of the "owned" matrix expect.
   if(useFEAssembly_)
     return getFEMatrix(i,j);

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
   // See getTpetraMatrix(): one shared FE matrix serves as both views.
   if(useFEAssembly_)
     return getFEMatrix(i,j);

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

namespace blocked_tpetra_lof_detail {

  /** Pull the Tpetra CrsMatrix out of block (i,j) of a Thyra blocked operator, or return
    * null if that block is excluded. Mirrors the extraction the container already does.
    */
  template <typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
  Teuchos::RCP<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >
  getBlockAsCrsMatrix(Thyra::PhysicallyBlockedLinearOpBase<ScalarT> & Amat,int i,int j)
  {
    using Teuchos::RCP;
    using Teuchos::rcp_dynamic_cast;

    RCP<Thyra::LinearOpBase<ScalarT> > block = Amat.getNonconstBlock(i,j);
    if(block==Teuchos::null)
      return Teuchos::null;

    RCP<Tpetra::Operator<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> > t_block =
        rcp_dynamic_cast<Thyra::TpetraLinearOp<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(block,true)->getTpetraOperator();

    return rcp_dynamic_cast<Tpetra::CrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT> >(t_block,true);
  }

}

// In FE mode the fill lifecycle is driven here rather than delegated to the container.
//
// The container's beginFill()/endFill() call the plain, inherited CrsMatrix::resumeFill() /
// fillComplete() on each block. For an FECrsMatrix that would silently bypass its owned /
// owned+shared state machine, so ghost-row contributions would never be migrated to the
// owned rows. The factory is the right place for the FE version: it is a single object
// shared by both containers and it already owns the per-block matrix cache, whereas the
// container is not templated on Traits and has no back-pointer to the factory.
template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
beginFill(LinearObjContainer & loc) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::PhysicallyBlockedLinearOpBase;

  BTLOC & tloc = Teuchos::dyn_cast<BTLOC>(loc);
  if(tloc.get_A()==Teuchos::null)
    return;

  if(!useFEAssembly_) {
    tloc.beginFill();
    return;
  }

  RCP<PhysicallyBlockedLinearOpBase<ScalarT> > Amat
      = rcp_dynamic_cast<PhysicallyBlockedLinearOpBase<ScalarT> >(tloc.get_A(),true);

  const int blockDim = static_cast<int>(gidProviders_.size());
  for(int i=0;i<blockDim;i++) {
    for(int j=0;j<blockDim;j++) {
      RCP<CrsMatrixType> mat
          = blocked_tpetra_lof_detail::getBlockAsCrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(*Amat,i,j);
      if(mat==Teuchos::null)
        continue;

      RCP<FECrsMatrixType> feMat = rcp_dynamic_cast<FECrsMatrixType>(mat);
      if(feMat==Teuchos::null) {
        // not an FE block (should not happen in FE mode, but stay well behaved)
        mat->resumeFill();
        continue;
      }

      // The owned and ghosted containers hold the SAME matrix per block, and
      // AssemblyEngine::evaluate() calls beginFill() on both. FECrsMatrix::beginAssembly()
      // asserts its fill state is "closed", so the second call would throw; tracking which
      // matrix already has an assembly open per block collapses the pair into the single
      // begin the FE state machine expects. This cannot be asked of the matrix directly --
      // see feAssemblyOpenOn_.
      const std::pair<int,int> key(i,j);
      typename std::unordered_map<std::pair<int,int>,const FECrsMatrixType *,panzer::pair_hash>::const_iterator
          openItr = feAssemblyOpenOn_.find(key);
      if(openItr!=feAssemblyOpenOn_.end() && openItr->second==feMat.get())
        continue;

      feMat->beginAssembly();

      // Clear ONLY the ghost rows of this block.
      //
      // Stale data can only accumulate there. Between assemblies the matrix rests in its
      // OWNED view, whose values alias just the leading chunk of the owned+shared array
      // (see Tpetra_FECrsMatrix_def.hpp, "we'll grab the first chunk of the Owned+Shared
      // matrix's values array"). A caller's setAllToScalar therefore reaches owned rows but
      // never the ghost rows, and endAssembly() does not clear them either, being a
      // combining self-export that leaves its source untouched. Without this the next
      // assembly sums onto the previous one's ghost contributions and inflates every
      // shared-interface dof.
      //
      // Deliberately NOT setAllToScalar(0.0) over the whole block: callers legitimately set
      // matrix values before calling beginFill (adjustForDirichletConditions is exercised
      // exactly that way), and wiping all of it would silently discard their data. Owned
      // rows are the caller's to manage; the ghost rows are scratch space that must start
      // each assembly at zero.
      //
      // Relies on the FECrsGraph invariant that owned rows are a prefix of the owned+shared
      // rows, which is what the V2 constructor guarantees.
      {
        auto rowptrs = feMat->getLocalRowPtrsHost();
        auto lclMat  = feMat->getLocalMatrixDevice();
        const size_t numOwnedRows = getMap(i)->getLocalNumElements();
        const size_t numRows      = static_cast<size_t>(lclMat.numRows());
        if(numOwnedRows < numRows) {
          const size_t firstGhostEntry = rowptrs(numOwnedRows);
          auto vals = lclMat.values;
          if(firstGhostEntry < vals.extent(0))
            Kokkos::deep_copy(Kokkos::subview(vals,Kokkos::make_pair(firstGhostEntry,vals.extent(0))),0.0);
        }
      }

      feAssemblyOpenOn_[key] = feMat.get();
    }
  }
}

template <typename Traits,typename ScalarT,typename LocalOrdinalT,typename GlobalOrdinalT,typename NodeT>
void BlockedTpetraLinearObjFactory<Traits,ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>::
endFill(LinearObjContainer & loc) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::PhysicallyBlockedLinearOpBase;

  BTLOC & tloc = Teuchos::dyn_cast<BTLOC>(loc);
  if(tloc.get_A()==Teuchos::null)
    return;

  if(!useFEAssembly_) {
    tloc.endFill();
    return;
  }

  RCP<PhysicallyBlockedLinearOpBase<ScalarT> > Amat
      = rcp_dynamic_cast<PhysicallyBlockedLinearOpBase<ScalarT> >(tloc.get_A(),true);

  const int blockDim = static_cast<int>(gidProviders_.size());
  for(int i=0;i<blockDim;i++) {
    for(int j=0;j<blockDim;j++) {
      RCP<CrsMatrixType> mat
          = blocked_tpetra_lof_detail::getBlockAsCrsMatrix<ScalarT,LocalOrdinalT,GlobalOrdinalT,NodeT>(*Amat,i,j);
      if(mat==Teuchos::null)
        continue;

      RCP<FECrsMatrixType> feMat = rcp_dynamic_cast<FECrsMatrixType>(mat);
      if(feMat==Teuchos::null) {
        mat->fillComplete(getMap(j),getMap(i));
        continue;
      }

      // See beginFill(): must go through endAssembly(), not the plain fillComplete(), so the
      // owned+shared -> owned cross-rank merge happens. This single endAssembly() IS the
      // ghost->global migration for this block, which is why ghostToGlobalThyraMatrix()
      // skips the export for shared FE blocks.
      //
      // Mirror of the beginFill() guard: endFill() is likewise called on both containers
      // holding the same matrix, and endAssembly() asserts its fill state is "open", so only
      // the first call may run it.
      const std::pair<int,int> key(i,j);
      typename std::unordered_map<std::pair<int,int>,const FECrsMatrixType *,panzer::pair_hash>::const_iterator
          openItr = feAssemblyOpenOn_.find(key);
      if(openItr!=feAssemblyOpenOn_.end() && openItr->second==feMat.get()) {
        feMat->endAssembly();
        feAssemblyOpenOn_.erase(key);
      }
    }
  }
}

}

#endif // __Panzer_BlockedTpetraLinearObjFactory_impl_hpp__
