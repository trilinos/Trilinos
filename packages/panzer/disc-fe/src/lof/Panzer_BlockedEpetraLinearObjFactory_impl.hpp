// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef   __Panzer_BlockedEpetraLinearObjFactory_impl_hpp__
#define   __Panzer_BlockedEpetraLinearObjFactory_impl_hpp__


// Epetra
#include "Epetra_CrsMatrix.h"
#include "Epetra_MpiComm.h"
#include "Epetra_MultiVector.h"
#include "Epetra_Vector.h"

// EpetraExt
#include "EpetraExt_VectorIn.h"
#include "EpetraExt_VectorOut.h"

// Panzer
#include "Panzer_BlockedVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_BlockedVector_Write_GlobalEvaluationData.hpp"
#include "Panzer_EpetraVector_ReadOnly_GlobalEvaluationData.hpp"
#include "Panzer_EpetraVector_Write_GlobalEvaluationData.hpp"
#include "Panzer_Filtered_GlobalIndexer.hpp"
#include "Panzer_HashUtils.hpp"
#include "Panzer_GlobalIndexer.hpp"

// Thyra
#include "Thyra_DefaultBlockedLinearOp.hpp"
#include "Thyra_DefaultProductVector.hpp"
#include "Thyra_DefaultProductVectorSpace.hpp"
#include "Thyra_EpetraLinearOp.hpp"
#include "Thyra_EpetraThyraWrappers.hpp"
#include "Thyra_get_Epetra_Operator.hpp"
#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"

using Teuchos::RCP;

namespace panzer {

// ************************************************************
// class BlockedEpetraLinearObjFactory
// ************************************************************

template <typename Traits,typename LocalOrdinalT>
BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
BlockedEpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                              const Teuchos::RCP<const GlobalIndexer> & gidProvider,
                              bool useDiscreteAdjoint)
   : useColGidProviders_(false), eComm_(Teuchos::null)
   , rawMpiComm_(comm->getRawMpiComm())
   , useDiscreteAdjoint_(useDiscreteAdjoint)
{ 
   rowDOFManagerContainer_ = Teuchos::rcp(new DOFManagerContainer(gidProvider));
   colDOFManagerContainer_ = rowDOFManagerContainer_;

   eComm_ = Teuchos::rcp(new Epetra_MpiComm(*rawMpiComm_));

   makeRoomForBlocks(rowDOFManagerContainer_->getFieldBlocks());

   // build and register the gather/scatter evaluators with 
   // the base class.
   this->buildGatherScatterEvaluators(*this);

   tComm_ = Teuchos::rcp(new Teuchos::MpiComm<int>(rawMpiComm_));
}

template <typename Traits,typename LocalOrdinalT>
BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
BlockedEpetraLinearObjFactory(const Teuchos::RCP<const Teuchos::MpiComm<int> > & comm,
                              const Teuchos::RCP<const GlobalIndexer> & gidProvider,
                              const Teuchos::RCP<const GlobalIndexer> & colGidProvider,
                              bool useDiscreteAdjoint)
   : eComm_(Teuchos::null)
   , rawMpiComm_(comm->getRawMpiComm())
   , useDiscreteAdjoint_(useDiscreteAdjoint)
{ 
   rowDOFManagerContainer_ = Teuchos::rcp(new DOFManagerContainer(gidProvider));
   colDOFManagerContainer_ = Teuchos::rcp(new DOFManagerContainer(colGidProvider));

   eComm_ = Teuchos::rcp(new Epetra_MpiComm(*rawMpiComm_));

   useColGidProviders_ = true;

   makeRoomForBlocks(rowDOFManagerContainer_->getFieldBlocks(),colDOFManagerContainer_->getFieldBlocks());

   // build and register the gather/scatter evaluators with 
   // the base class.
   this->buildGatherScatterEvaluators(*this);

   tComm_ = Teuchos::rcp(new Teuchos::MpiComm<int>(rawMpiComm_));
}

template <typename Traits,typename LocalOrdinalT>
BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::~BlockedEpetraLinearObjFactory()
{ }

// LinearObjectFactory functions 
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename LocalOrdinalT>
void 
BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
readVector(const std::string & identifier,LinearObjContainer & loc,int id) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::dyn_cast;
  using Thyra::ProductVectorBase;

  BlockedEpetraLinearObjContainer & eloc = dyn_cast<BlockedEpetraLinearObjContainer>(loc);

   // extract the vector from linear object container
  RCP<Thyra::VectorBase<double> > vec;
  switch(id) {
  case LinearObjContainer::X:
    vec = eloc.get_x();
    break;
  case LinearObjContainer::DxDt:
    vec = eloc.get_dxdt();
    break;
  case LinearObjContainer::F:
    vec = eloc.get_f();
    break;
  default:
    TEUCHOS_ASSERT(false);
    break;
  };

  int blockRows = this->getBlockRowCount();
  RCP<ProductVectorBase<double> > b_vec = Thyra::nonconstProductVectorBase(vec);

  // convert to Epetra then write out each vector to file
  for(int i=0;i<blockRows;i++) {
    RCP<Thyra::VectorBase<double> > x = b_vec->getNonconstVectorBlock(i); 
    RCP<Epetra_Vector> ex = Thyra::get_Epetra_Vector(*getMap(i),x);

    // build the file name from the identifier
    std::stringstream ss;
    ss << identifier << "-" << i << ".mm";

    // read in vector (wow the MM to Vector is a poorly designed interface!)
    Epetra_Vector * ptr_ex = 0;
    TEUCHOS_ASSERT(0==EpetraExt::MatrixMarketFileToVector(ss.str().c_str(),*getMap(i),ptr_ex));

    *ex = *ptr_ex;
    delete ptr_ex;
  }
}

template <typename Traits,typename LocalOrdinalT>
void 
BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
writeVector(const std::string & identifier,const LinearObjContainer & loc,int id) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::dyn_cast;
  using Thyra::ProductVectorBase;

  const BlockedEpetraLinearObjContainer & eloc = dyn_cast<const BlockedEpetraLinearObjContainer>(loc);

   // extract the vector from linear object container
  RCP<const Thyra::VectorBase<double> > vec;
  switch(id) {
  case LinearObjContainer::X:
    vec = eloc.get_x();
    break;
  case LinearObjContainer::DxDt:
    vec = eloc.get_dxdt();
    break;
  case LinearObjContainer::F:
    vec = eloc.get_f();
    break;
  default:
    TEUCHOS_ASSERT(false);
    break;
  };

  int blockRows = this->getBlockRowCount();
  RCP<const ProductVectorBase<double> > b_vec = Thyra::productVectorBase(vec);

  // convert to Epetra then write out each vector to file
  for(int i=0;i<blockRows;i++) {
    RCP<const Thyra::VectorBase<double> > x = b_vec->getVectorBlock(i); 
    RCP<const Epetra_Vector> ex = Thyra::get_Epetra_Vector(*getMap(i),x);

    // build the file name from the identifier
    std::stringstream ss;
    ss << identifier << "-" << i << ".mm";

    // write out vector
    TEUCHOS_ASSERT(0==EpetraExt::VectorToMatrixMarketFile(ss.str().c_str(),*ex));
  }
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<LinearObjContainer> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::buildLinearObjContainer() const
{
   // if a "single field" DOFManager is used
   if(!rowDOFManagerContainer_->containsBlockedDOFManager() && !colDOFManagerContainer_->containsBlockedDOFManager()) {
     Teuchos::RCP<EpetraLinearObjContainer> container = Teuchos::rcp(new EpetraLinearObjContainer(getColMap(0),getMap(0)));

     return container;
   }

   std::vector<Teuchos::RCP<const Epetra_Map> > blockMaps;
   std::size_t blockDim = getBlockRowCount();
   for(std::size_t i=0;i<blockDim;i++)
      blockMaps.push_back(getMap(i));

   Teuchos::RCP<BlockedEpetraLinearObjContainer > container = Teuchos::rcp(new BlockedEpetraLinearObjContainer);
   container->setMapsForBlocks(blockMaps);

   return container;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<LinearObjContainer> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::buildGhostedLinearObjContainer() const
{
   // if a "single field" DOFManager is used
   if(!rowDOFManagerContainer_->containsBlockedDOFManager() && !colDOFManagerContainer_->containsBlockedDOFManager()) {
     Teuchos::RCP<EpetraLinearObjContainer> container = Teuchos::rcp(new EpetraLinearObjContainer(getGhostedColMap(0),getGhostedMap(0)));

     return container;
   }

   std::vector<Teuchos::RCP<const Epetra_Map> > blockMaps;
   std::size_t blockDim = getBlockRowCount();
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

   if(   !rowDOFManagerContainer_->containsBlockedDOFManager()
      && !colDOFManagerContainer_->containsBlockedDOFManager()) {
     const EpetraLinearObjContainer & e_in = Teuchos::dyn_cast<const EpetraLinearObjContainer>(in); 
     EpetraLinearObjContainer & e_out = Teuchos::dyn_cast<EpetraLinearObjContainer>(out); 
  
     // Operations occur if the GLOBAL container has the correct targets!
     // Users set the GLOBAL continer arguments
     if ( !is_null(e_in.get_x()) && !is_null(e_out.get_x()) && ((mem & LOC::X)==LOC::X))
       globalToGhostEpetraVector(0,*e_in.get_x(),*e_out.get_x(),true);
  
     if ( !is_null(e_in.get_dxdt()) && !is_null(e_out.get_dxdt()) && ((mem & LOC::DxDt)==LOC::DxDt))
       globalToGhostEpetraVector(0,*e_in.get_dxdt(),*e_out.get_dxdt(),true);

     if ( !is_null(e_in.get_f()) && !is_null(e_out.get_f()) && ((mem & LOC::F)==LOC::F))
       globalToGhostEpetraVector(0,*e_in.get_f(),*e_out.get_f(),false);
   }
   else {
     const BLOC & b_in = Teuchos::dyn_cast<const BLOC>(in); 
     BLOC & b_out = Teuchos::dyn_cast<BLOC>(out); 
  
     // Operations occur if the GLOBAL container has the correct targets!
     // Users set the GLOBAL continer arguments
     if ( !is_null(b_in.get_x()) && !is_null(b_out.get_x()) && ((mem & LOC::X)==LOC::X))
       globalToGhostThyraVector(b_in.get_x(),b_out.get_x(),true);
  
     if ( !is_null(b_in.get_dxdt()) && !is_null(b_out.get_dxdt()) && ((mem & LOC::DxDt)==LOC::DxDt))
       globalToGhostThyraVector(b_in.get_dxdt(),b_out.get_dxdt(),true);

     if ( !is_null(b_in.get_f()) && !is_null(b_out.get_f()) && ((mem & LOC::F)==LOC::F))
        globalToGhostThyraVector(b_in.get_f(),b_out.get_f(),false);
   }
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::ghostToGlobalContainer(const LinearObjContainer & in,
                                                                          LinearObjContainer & out,int mem) const
{
   using Teuchos::is_null;

   typedef LinearObjContainer LOC;
   typedef BlockedEpetraLinearObjContainer BLOC;

   if(   !rowDOFManagerContainer_->containsBlockedDOFManager()
      && !colDOFManagerContainer_->containsBlockedDOFManager()) {
     const EpetraLinearObjContainer & e_in = Teuchos::dyn_cast<const EpetraLinearObjContainer>(in); 
     EpetraLinearObjContainer & e_out = Teuchos::dyn_cast<EpetraLinearObjContainer>(out); 

     // Operations occur if the GLOBAL container has the correct targets!
     // Users set the GLOBAL continer arguments
     if ( !is_null(e_in.get_x()) && !is_null(e_out.get_x()) && ((mem & LOC::X)==LOC::X))
       ghostToGlobalEpetraVector(0,*e_in.get_x(),*e_out.get_x(),true);

     if ( !is_null(e_in.get_f()) && !is_null(e_out.get_f()) && ((mem & LOC::F)==LOC::F))
       ghostToGlobalEpetraVector(0,*e_in.get_f(),*e_out.get_f(),false);

     if ( !is_null(e_in.get_A()) && !is_null(e_out.get_A()) && ((mem & LOC::Mat)==LOC::Mat))
       ghostToGlobalEpetraMatrix(0,*e_in.get_A(),*e_out.get_A());
   }
   else {
     const BLOC & b_in = Teuchos::dyn_cast<const BLOC>(in); 
     BLOC & b_out = Teuchos::dyn_cast<BLOC>(out); 

     // Operations occur if the GLOBAL container has the correct targets!
     // Users set the GLOBAL continer arguments
     if ( !is_null(b_in.get_x()) && !is_null(b_out.get_x()) && ((mem & LOC::X)==LOC::X))
       ghostToGlobalThyraVector(b_in.get_x(),b_out.get_x(),true);

     if ( !is_null(b_in.get_f()) && !is_null(b_out.get_f()) && ((mem & LOC::F)==LOC::F))
       ghostToGlobalThyraVector(b_in.get_f(),b_out.get_f(),false);

     if ( !is_null(b_in.get_A()) && !is_null(b_out.get_A()) && ((mem & LOC::Mat)==LOC::Mat))
       ghostToGlobalThyraMatrix(*b_in.get_A(),*b_out.get_A());
   }
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
adjustForDirichletConditions(const LinearObjContainer & localBCRows,
                             const LinearObjContainer & globalBCRows,
                             LinearObjContainer & ghostedObjs,
                             bool zeroVectorRows, bool adjustX) const
{
   typedef ThyraObjContainer<double> TOC;

   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::LinearOpBase;
   using Thyra::PhysicallyBlockedLinearOpBase;
   using Thyra::VectorBase;
   using Thyra::ProductVectorBase;
   using Thyra::get_Epetra_Vector;
   using Thyra::get_Epetra_Operator;

   int rBlockDim = getBlockRowCount();
   int cBlockDim = getBlockColCount();

   // first cast to block LOCs
   const TOC & b_localBCRows = Teuchos::dyn_cast<const TOC>(localBCRows); 
   const TOC & b_globalBCRows = Teuchos::dyn_cast<const TOC>(globalBCRows); 
   TOC & b_ghosted = Teuchos::dyn_cast<TOC>(ghostedObjs); 

   TEUCHOS_ASSERT(b_localBCRows.get_f_th()!=Teuchos::null);
   TEUCHOS_ASSERT(b_globalBCRows.get_f_th()!=Teuchos::null);

   // cast each component as needed to their product form
   RCP<PhysicallyBlockedLinearOpBase<double> > A = rcp_dynamic_cast<PhysicallyBlockedLinearOpBase<double> >(b_ghosted.get_A_th());
   if(A==Teuchos::null && b_ghosted.get_A_th()!=Teuchos::null) {
     // assume it isn't physically blocked, for convenience physically block it
     A = rcp_dynamic_cast<PhysicallyBlockedLinearOpBase<double> >(Thyra::nonconstBlock1x1(b_ghosted.get_A_th()));
   }

   RCP<ProductVectorBase<double> > f          = b_ghosted.get_f_th()==Teuchos::null
                                                            ? Teuchos::null 
                                                            : Thyra::castOrCreateNonconstProductVectorBase(b_ghosted.get_f_th());
   RCP<ProductVectorBase<double> > local_bcs  = b_localBCRows.get_f_th()==Teuchos::null
                                                            ? Teuchos::null 
                                                            : Thyra::castOrCreateNonconstProductVectorBase(b_localBCRows.get_f_th());
   RCP<ProductVectorBase<double> > global_bcs = b_globalBCRows.get_f_th()==Teuchos::null
                                                            ? Teuchos::null 
                                                            : Thyra::castOrCreateNonconstProductVectorBase(b_globalBCRows.get_f_th());

   if(adjustX) f = Thyra::castOrCreateNonconstProductVectorBase(b_ghosted.get_x_th());

   // sanity check!
   if(A!=Teuchos::null) TEUCHOS_ASSERT(A->productRange()->numBlocks()==rBlockDim);
   if(A!=Teuchos::null) TEUCHOS_ASSERT(A->productDomain()->numBlocks()==cBlockDim);
   if(f!=Teuchos::null) TEUCHOS_ASSERT(f->productSpace()->numBlocks()==rBlockDim);
   TEUCHOS_ASSERT(local_bcs->productSpace()->numBlocks()==rBlockDim);
   TEUCHOS_ASSERT(global_bcs->productSpace()->numBlocks()==rBlockDim);

   for(int i=0;i<rBlockDim;i++) {
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

      for(int j=0;j<cBlockDim;j++) {

         // pull out epetra values
         RCP<LinearOpBase<double> > th_A = (A== Teuchos::null)? Teuchos::null : A->getNonconstBlock(i,j);
 
         // don't do anyting if opertor is null
         RCP<Epetra_CrsMatrix> e_A;
         if(th_A==Teuchos::null)
            e_A = Teuchos::null;
         else 
            e_A = rcp_dynamic_cast<Epetra_CrsMatrix>(get_Epetra_Operator(*th_A),true);

         // adjust Block operator
         adjustForDirichletConditions(*e_local_bcs,*e_global_bcs,e_f.ptr(),e_A.ptr(),zeroVectorRows);

         e_f = Teuchos::null; // this is so we only adjust it once on the first pass
      }
   }
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
adjustForDirichletConditions(const Epetra_Vector & local_bcs,
                             const Epetra_Vector & global_bcs,
                             const Teuchos::Ptr<Epetra_Vector> & f,
                             const Teuchos::Ptr<Epetra_CrsMatrix> & A,
                             bool zeroVectorRows) const
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

      if(local_bcs[i]==0.0 || zeroVectorRows) { 
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
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
applyDirichletBCs(const LinearObjContainer & counter,
                  LinearObjContainer & result) const
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Teuchos::dyn_cast;

  typedef Thyra::ProductVectorBase<double> PVector;

  const ThyraObjContainer<double> & th_counter = dyn_cast<const ThyraObjContainer<double> >(counter);
  ThyraObjContainer<double> & th_result  = dyn_cast<ThyraObjContainer<double> >(result);
  
  RCP<const PVector> count = Thyra::castOrCreateProductVectorBase(th_counter.get_f_th().getConst());
  RCP<const PVector> f_in  = Thyra::castOrCreateProductVectorBase(th_counter.get_f_th().getConst());
  RCP<PVector> f_out       = Thyra::castOrCreateNonconstProductVectorBase(th_result.get_f_th());

  int rBlockDim = getBlockRowCount();
  for(int i=0;i<rBlockDim;i++) {

    Teuchos::ArrayRCP<const double> count_array,f_in_array;
    Teuchos::ArrayRCP<double> f_out_array;

    rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(count->getVectorBlock(i),true)->getLocalData(Teuchos::ptrFromRef(count_array));
    rcp_dynamic_cast<const Thyra::SpmdVectorBase<double> >(f_in->getVectorBlock(i),true)->getLocalData(Teuchos::ptrFromRef(f_in_array));
    rcp_dynamic_cast<Thyra::SpmdVectorBase<double> >(f_out->getNonconstVectorBlock(i),true)->getNonconstLocalData(Teuchos::ptrFromRef(f_out_array));

    TEUCHOS_ASSERT(count_array.size()==f_in_array.size());
    TEUCHOS_ASSERT(count_array.size()==f_out_array.size());

    for(Teuchos::ArrayRCP<double>::size_type j=0;j<count_array.size();++j) {
      if(count_array[j]!=0.0)
        f_out_array[j] = f_in_array[j];
    }
  }
}

///////////////////////////////////////////////////////////////////////////////
//
//  buildReadOnlyDomainContainer()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
Teuchos::RCP<ReadOnlyVector_GlobalEvaluationData>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
buildReadOnlyDomainContainer() const
{
  using std::vector;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using BVROGED = panzer::BlockedVector_ReadOnly_GlobalEvaluationData;
  using EVROGED = panzer::EpetraVector_ReadOnly_GlobalEvaluationData;
  using ROVGED  = panzer::ReadOnlyVector_GlobalEvaluationData;

  // If a "single field" DOFManager is used, return a single
  // EpetraVector_ReadOnly_GlobalEvaluationData.
  if (not colDOFManagerContainer_->containsBlockedDOFManager())
  {
    auto ged = rcp(new EVROGED);
    ged->initialize(getGhostedColImport2(0), getGhostedColMap2(0),
      getColMap(0));
    return ged;
  } // end if a "single field" DOFManager is used

  // Otherwise, return a BlockedVector_ReadOnly_GlobalEvaluationData.
  vector<RCP<ROVGED>> gedBlocks;
  for (int i(0); i < getBlockColCount(); ++i)
  {
    auto vecGed = rcp(new EVROGED);
    vecGed->initialize(getGhostedColImport2(i), getGhostedColMap2(i),
      getColMap(i));
    gedBlocks.push_back(vecGed);
  } // end loop over the blocks
  auto ged = rcp(new BVROGED);
  ged->initialize(getGhostedThyraDomainSpace2(), getThyraDomainSpace(),
    gedBlocks);
  return ged;
} // end of buildReadOnlyDomainContainer()

///////////////////////////////////////////////////////////////////////////////
//
//  buildWriteDomainContainer()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
Teuchos::RCP<WriteVector_GlobalEvaluationData>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
buildWriteDomainContainer() const
{
  using std::vector;
  using Teuchos::RCP;
  using Teuchos::rcp;
  using BVWGED = panzer::BlockedVector_Write_GlobalEvaluationData;
  using EVWGED = panzer::EpetraVector_Write_GlobalEvaluationData;
  using WVGED  = panzer::WriteVector_GlobalEvaluationData;

  // If a "single field" DOFManager is used, return a single
  // EpetraVector_Write_GlobalEvaluationData.
  if (not colDOFManagerContainer_->containsBlockedDOFManager())
  {
    auto ged = rcp(new EVWGED);
    ged->initialize(getGhostedColExport2(0), getGhostedColMap2(0),
      getColMap(0));
    return ged;
  } // end if a "single field" DOFManager is used

  // Otherwise, return a BlockedVector_Write_GlobalEvaluationData.
  vector<RCP<WVGED>> gedBlocks;
  for (int i(0); i < getBlockColCount(); ++i)
  {
    auto vecGed = rcp(new EVWGED);
    vecGed->initialize(getGhostedColExport2(i), getGhostedColMap2(i),
      getColMap(i));
    gedBlocks.push_back(vecGed);
  } // end loop over the blocks
  auto ged = rcp(new BVWGED);
  ged->initialize(getGhostedThyraDomainSpace2(), getThyraDomainSpace(),
    gedBlocks);
  return ged;
} // end of buildWriteDomainContainer()

template <typename Traits,typename LocalOrdinalT>
Teuchos::MpiComm<int> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getComm() const
{
   return *tComm_;
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
initializeContainer(int mem,LinearObjContainer & loc) const
{
   typedef ThyraObjContainer<double> TOC;

   TOC & toc = Teuchos::dyn_cast<TOC>(loc);
   initializeContainer_internal(mem,toc);
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
initializeGhostedContainer(int mem,LinearObjContainer & loc) const
{
   typedef LinearObjContainer LOC;
   typedef ThyraObjContainer<double> TOC;

   TOC & toc = Teuchos::dyn_cast<TOC>(loc);
   initializeGhostedContainer_internal(mem,toc);

   if(rowDOFManagerContainer_->containsBlockedDOFManager()) {
     typedef BlockedEpetraLinearObjContainer BLOC;

     BLOC & bloc = Teuchos::dyn_cast<BLOC>(loc);

     if((mem & LOC::F) == LOC::F)
       bloc.setRequiresDirichletAdjustment(true);

     if((mem & LOC::Mat) == LOC::Mat) 
       bloc.setRequiresDirichletAdjustment(true);
   }
   else {
     EpetraLinearObjContainer & eloc = Teuchos::dyn_cast<EpetraLinearObjContainer>(loc);

     if((mem & LOC::F) == LOC::F)
       eloc.setRequiresDirichletAdjustment(true);

     if((mem & LOC::Mat) == LOC::Mat) 
       eloc.setRequiresDirichletAdjustment(true);
   }
}

// Generic methods 
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
initializeContainer_internal(int mem,ThyraObjContainer<double> & loc) const
{
   typedef LinearObjContainer LOC;

   loc.clear();

   if((mem & LOC::X) == LOC::X)
      loc.set_x_th(getThyraDomainVector());

   if((mem & LOC::DxDt) == LOC::DxDt)
      loc.set_dxdt_th(getThyraDomainVector());
    
   if((mem & LOC::F) == LOC::F)
      loc.set_f_th(getThyraRangeVector());

   if((mem & LOC::Mat) == LOC::Mat)
      loc.set_A_th(getThyraMatrix());
}

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
initializeGhostedContainer_internal(int mem,ThyraObjContainer<double> & loc) const
{
   typedef LinearObjContainer LOC;

   loc.clear();

   if((mem & LOC::X) == LOC::X)
      loc.set_x_th(getGhostedThyraDomainVector());

   if((mem & LOC::DxDt) == LOC::DxDt)
      loc.set_dxdt_th(getGhostedThyraDomainVector());
    
   if((mem & LOC::F) == LOC::F)
      loc.set_f_th(getGhostedThyraRangeVector());

   if((mem & LOC::Mat) == LOC::Mat)
      loc.set_A_th(getGhostedThyraMatrix());
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
Teuchos::RCP<const GlobalIndexer> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGlobalIndexer(int i) const
{
   return rowDOFManagerContainer_->getFieldDOFManagers()[i];
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<const GlobalIndexer> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getColGlobalIndexer(int i) const
{
   return colDOFManagerContainer_->getFieldDOFManagers()[i];
}

///////////////////////////////////////////////////////////////////////////////
//
//  makeRoomForBlocks()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
void
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
makeRoomForBlocks(
  std::size_t blockCnt,
  std::size_t colBlockCnt)
{
  maps_.resize(blockCnt); 
  ghostedMaps_.resize(blockCnt); 
  ghostedMaps2_.resize(blockCnt); 
  importers_.resize(blockCnt); 
  importers2_.resize(blockCnt); 
  exporters_.resize(blockCnt); 
  if (colBlockCnt > 0)
  {
    colMaps_.resize(colBlockCnt); 
    colGhostedMaps_.resize(colBlockCnt); 
    colGhostedMaps2_.resize(colBlockCnt); 
    colImporters_.resize(colBlockCnt); 
    colImporters2_.resize(colBlockCnt); 
    colExporters_.resize(colBlockCnt); 
  } // end if (colBlockCnt > 0)
} // end of makeRoomForBlocks()

// Thyra methods 
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getThyraDomainSpace() const
{
   if(domainSpace_==Teuchos::null) {
     if(colDOFManagerContainer_->containsBlockedDOFManager()) {
       // loop over all vectors and build the vector space
       std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double> > > vsArray;
       for(int i=0;i<getBlockColCount();i++)  
         vsArray.push_back(Thyra::create_VectorSpace(getColMap(i)));

       domainSpace_ = Thyra::productVectorSpace<double>(vsArray);
     }
     else {
       // the domain space is not blocked (just an SPMD vector), build it from
       // the zeroth column
       domainSpace_ = Thyra::create_VectorSpace(getColMap(0));
     }
   }
   
   return domainSpace_;
}

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getThyraRangeSpace() const
{
   if(rangeSpace_==Teuchos::null) {
     if(rowDOFManagerContainer_->containsBlockedDOFManager()) {
       // loop over all vectors and build the vector space
       std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double> > > vsArray;
       for(int i=0;i<getBlockRowCount();i++)  
          vsArray.push_back(Thyra::create_VectorSpace(getMap(i)));
 
       rangeSpace_ = Thyra::productVectorSpace<double>(vsArray);
     }
     else {
       // the range space is not blocked (just an SPMD vector), build it from
       // the zeroth row
       rangeSpace_ = Thyra::create_VectorSpace(getMap(0));
     }
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
   // return a flat epetra matrix
   if(!rowDOFManagerContainer_->containsBlockedDOFManager() && 
      !colDOFManagerContainer_->containsBlockedDOFManager()) {
     return Thyra::nonconstEpetraLinearOp(getEpetraMatrix(0,0));
   }

   Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > blockedOp = Thyra::defaultBlockedLinearOp<double>();

   // get the block dimension
   std::size_t rBlockDim = getBlockRowCount();
   std::size_t cBlockDim = getBlockColCount();

   // this operator will be square
   blockedOp->beginBlockFill(rBlockDim,cBlockDim);

   // loop over each block
   for(std::size_t i=0;i<rBlockDim;i++) { 
      for(std::size_t j=0;j<cBlockDim;j++) {
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

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedThyraDomainSpace()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedThyraDomainSpace() const
{
  using std::vector;
  using Teuchos::RCP;
  using Thyra::create_VectorSpace;
  using Thyra::productVectorSpace;
  using Thyra::VectorSpaceBase;
  if (ghostedDomainSpace_.is_null())
  {
    if (colDOFManagerContainer_->containsBlockedDOFManager())
    {
      // Loop over all vectors and build the vector space.
      vector<RCP<const VectorSpaceBase<double>>> vsArray;
      for (int i(0); i < getBlockColCount(); ++i)
        vsArray.push_back(create_VectorSpace(getGhostedColMap(i)));
      ghostedDomainSpace_ = productVectorSpace<double>(vsArray);
    }
    else // if (not colDOFManagerContainer_->containsBlockedDOFManager())
    {
      // The domain space is not blocked (that is, we're just dealing with a
      // SPMD vector), so build it from the zeroth column.
      ghostedDomainSpace_ = create_VectorSpace(getGhostedColMap(0));
    } // end if (colDOFManagerContainer_->containsBlockedDOFManager()) or not
  } // end if (ghostedDomainSpace_.is_null())
  return ghostedDomainSpace_;
} // end of getGhostedThyraDomainSpace()

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedThyraDomainSpace2()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
Teuchos::RCP<const Thyra::VectorSpaceBase<double>>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedThyraDomainSpace2() const
{
  using std::vector;
  using Teuchos::RCP;
  using Thyra::create_VectorSpace;
  using Thyra::productVectorSpace;
  using Thyra::VectorSpaceBase;
  if (ghostedDomainSpace_.is_null())
  {
    if (colDOFManagerContainer_->containsBlockedDOFManager())
    {
      // Loop over all vectors and build the vector space.
      vector<RCP<const VectorSpaceBase<double>>> vsArray;
      for (int i(0); i < getBlockColCount(); ++i)
        vsArray.push_back(create_VectorSpace(getGhostedColMap2(i)));
      ghostedDomainSpace_ = productVectorSpace<double>(vsArray);
    }
    else // if (not colDOFManagerContainer_->containsBlockedDOFManager())
    {
      // The domain space is not blocked (that is, we're just dealing with a
      // SPMD vector), so build it from the zeroth column.
      ghostedDomainSpace_ = create_VectorSpace(getGhostedColMap2(0));
    } // end if (colDOFManagerContainer_->containsBlockedDOFManager()) or not
  } // end if (ghostedDomainSpace_.is_null())
  return ghostedDomainSpace_;
} // end of getGhostedThyraDomainSpace2()

template <typename Traits,typename LocalOrdinalT>
Teuchos::RCP<const Thyra::VectorSpaceBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGhostedThyraRangeSpace() const
{
   if(ghostedRangeSpace_==Teuchos::null) {
     if(rowDOFManagerContainer_->containsBlockedDOFManager()) {
       // loop over all vectors and build the vector space
       std::vector<Teuchos::RCP<const Thyra::VectorSpaceBase<double> > > vsArray;
       for(int i=0;i<getBlockRowCount();i++)  
         vsArray.push_back(Thyra::create_VectorSpace(getGhostedMap(i)));

       ghostedRangeSpace_ = Thyra::productVectorSpace<double>(vsArray);
     }
     else {
       // the range space is not blocked (just an SPMD vector), build it from
       // the zeroth row
       ghostedRangeSpace_ = Thyra::create_VectorSpace(getGhostedMap(0));
     }
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
Teuchos::RCP<Thyra::LinearOpBase<double> > BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGhostedThyraMatrix() const
{
   // return a flat epetra matrix
   if(!rowDOFManagerContainer_->containsBlockedDOFManager() && 
      !colDOFManagerContainer_->containsBlockedDOFManager()) {
     return Thyra::nonconstEpetraLinearOp(getGhostedEpetraMatrix(0,0));
   }

   Teuchos::RCP<Thyra::PhysicallyBlockedLinearOpBase<double> > blockedOp = Thyra::defaultBlockedLinearOp<double>();

   // get the block dimension
   std::size_t rBlockDim = getBlockRowCount();
   std::size_t cBlockDim = getBlockColCount();

   // this operator will be square
   blockedOp->beginBlockFill(rBlockDim,cBlockDim);

   // loop over each block
   for(std::size_t i=0;i<rBlockDim;i++) { 
      for(std::size_t j=0;j<cBlockDim;j++) {
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
                         const Teuchos::RCP<Thyra::VectorBase<double> > & out,bool col) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;
   using Thyra::get_Epetra_Vector;

   std::size_t blockDim = col ? getBlockColCount() : getBlockRowCount();

   // get product vectors
   RCP<const ProductVectorBase<double> > prod_in = Thyra::castOrCreateProductVectorBase(in);
   RCP<ProductVectorBase<double> > prod_out      = Thyra::castOrCreateNonconstProductVectorBase(out);

   TEUCHOS_ASSERT(prod_in->productSpace()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(prod_out->productSpace()->numBlocks()==(int) blockDim);

   for(std::size_t i=0;i<blockDim;i++) {
      // first get each Epetra vector
      RCP<const Epetra_Vector> ep_in;
      RCP<Epetra_Vector> ep_out;

      if(not col) {
        ep_in = get_Epetra_Vector(*getGhostedMap(i),prod_in->getVectorBlock(i));
        ep_out = get_Epetra_Vector(*getMap(i),prod_out->getNonconstVectorBlock(i));
      } else {
        ep_in = get_Epetra_Vector(*getGhostedColMap(i),prod_in->getVectorBlock(i));
        ep_out = get_Epetra_Vector(*getColMap(i),prod_out->getNonconstVectorBlock(i));
      }

      // use Epetra to do global communication
      ghostToGlobalEpetraVector(i,*ep_in,*ep_out,col);
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

   // get the block dimension
   std::size_t rBlockDim = getBlockRowCount();
   std::size_t cBlockDim = getBlockColCount();

   // get product vectors
   const PhysicallyBlockedLinearOpBase<double> & prod_in = dyn_cast<const PhysicallyBlockedLinearOpBase<double> >(in);
   PhysicallyBlockedLinearOpBase<double> & prod_out      = dyn_cast<PhysicallyBlockedLinearOpBase<double> >(out);

   TEUCHOS_ASSERT(prod_in.productRange()->numBlocks()==(int) rBlockDim);
   TEUCHOS_ASSERT(prod_in.productDomain()->numBlocks()==(int) cBlockDim);
   TEUCHOS_ASSERT(prod_out.productRange()->numBlocks()==(int) rBlockDim);
   TEUCHOS_ASSERT(prod_out.productDomain()->numBlocks()==(int) cBlockDim);

   for(std::size_t i=0;i<rBlockDim;i++) {
      for(std::size_t j=0;j<cBlockDim;j++) {
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
                         const Teuchos::RCP<Thyra::VectorBase<double> > & out,bool col) const
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::ProductVectorBase;
   using Thyra::get_Epetra_Vector;

   std::size_t blockDim = col ? getBlockColCount() : getBlockRowCount();

   // get product vectors
   RCP<const ProductVectorBase<double> > prod_in = Thyra::castOrCreateProductVectorBase(in);
   RCP<ProductVectorBase<double> > prod_out      = Thyra::castOrCreateNonconstProductVectorBase(out);

   TEUCHOS_ASSERT(prod_in->productSpace()->numBlocks()==(int) blockDim);
   TEUCHOS_ASSERT(prod_out->productSpace()->numBlocks()==(int) blockDim);

   for(std::size_t i=0;i<blockDim;i++) {
      // first get each Epetra vector
      RCP<const Epetra_Vector> ep_in; 
      RCP<Epetra_Vector> ep_out;      

      if(not col) {
        ep_in = get_Epetra_Vector(*getMap(i),prod_in->getVectorBlock(i));
        ep_out = get_Epetra_Vector(*getGhostedMap(i),prod_out->getNonconstVectorBlock(i));
      }
      else {
        ep_in = get_Epetra_Vector(*getColMap(i),prod_in->getVectorBlock(i));
        ep_out = get_Epetra_Vector(*getGhostedColMap(i),prod_out->getNonconstVectorBlock(i));
      }

      // use Epetra to do global communication
      globalToGhostEpetraVector(i,*ep_in,*ep_out,col);
   }
}

// Epetra methods 
/////////////////////////////////////////////////////////////////////

template <typename Traits,typename LocalOrdinalT>
void BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
ghostToGlobalEpetraVector(int i,const Epetra_Vector & in,Epetra_Vector & out,bool col) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<Epetra_Export> exporter = col ? getGhostedColExport(i) : getGhostedExport(i);
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
globalToGhostEpetraVector(int i,const Epetra_Vector & in,Epetra_Vector & out,bool col) const
{
   using Teuchos::RCP;

   // do the global distribution
   RCP<Epetra_Import> importer = col ? getGhostedColImport(i) : getGhostedImport(i);
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
      maps_[i] = buildMap(i);

   return maps_[i];
}

// get the map from the matrix
template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getColMap(int i) const
{
   if(not useColGidProviders_)
     return getMap(i); 

   if(colMaps_[i]==Teuchos::null) 
      colMaps_[i] = buildColMap(i);

   return colMaps_[i];
}

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedMap()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedMap(
  int i) const
{
  if (ghostedMaps_[i].is_null())
    ghostedMaps_[i] = buildGhostedMap(i);
  return ghostedMaps_[i];
} // end of getGhostedMap()

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedMap2()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedMap2(
  int i) const
{
  if (ghostedMaps2_[i].is_null())
    ghostedMaps2_[i] = buildGhostedMap2(i);
  return ghostedMaps2_[i];
} // end of getGhostedMap2()

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedColMap()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedColMap(
  int i) const
{
  if (not useColGidProviders_)
    return getGhostedMap(i); 
  if (colGhostedMaps_[i].is_null())
    colGhostedMaps_[i] = buildColGhostedMap(i);
  return colGhostedMaps_[i];
} // end of getGhostedColMap()

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedColMap2()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedColMap2(
  int i) const
{
  if (not useColGidProviders_)
    return getGhostedMap2(i); 
  if (colGhostedMaps2_[i].is_null())
    colGhostedMaps2_[i] = buildColGhostedMap2(i);
  return colGhostedMaps2_[i];
} // end of getGhostedColMap2()

// get the graph of the crs matrix
template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getGraph(int i,int j) const
{
  typedef std::unordered_map<std::pair<int,int>,Teuchos::RCP<Epetra_CrsGraph>,panzer::pair_hash> GraphMap;
   
   GraphMap::const_iterator itr = graphs_.find(std::make_pair(i,j));
   Teuchos::RCP<Epetra_CrsGraph> graph;
   if(itr==graphs_.end()) {
      graph = buildGraph(i,j);
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
   typedef std::unordered_map<std::pair<int,int>,Teuchos::RCP<Epetra_CrsGraph>,panzer::pair_hash> GraphMap;
   
   GraphMap::const_iterator itr = ghostedGraphs_.find(std::make_pair(i,j));
   Teuchos::RCP<Epetra_CrsGraph> ghostedGraph;
   if(itr==ghostedGraphs_.end()) {
      ghostedGraph = buildGhostedGraph(i,j,true);
      ghostedGraphs_[std::make_pair(i,j)] = ghostedGraph;
   }
   else
      ghostedGraph = itr->second;

   TEUCHOS_ASSERT(ghostedGraph!=Teuchos::null);
   return ghostedGraph;
}

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedImport()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Import>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedImport(
  int i) const
{
  using Teuchos::rcp;
  if (importers_[i].is_null())
    importers_[i] = rcp(new Epetra_Import(*getGhostedMap(i), *getMap(i)));
  return importers_[i];
} // end of getGhostedImport()

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedImport2()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Import>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedImport2(
  int i) const
{
  using Teuchos::rcp;
  if (importers2_[i].is_null())
    importers2_[i] = rcp(new Epetra_Import(*getGhostedMap2(i), *getMap(i)));
  return importers2_[i];
} // end of getGhostedImport2()

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedColImport()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Import>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedColImport(
  int i) const
{
  using Teuchos::rcp;
  if (not useColGidProviders_)
    return getGhostedImport(i);
  if (colImporters_[i].is_null())
    colImporters_[i] =
      rcp(new Epetra_Import(*getGhostedColMap(i), *getColMap(i)));
  return colImporters_[i];
} // end of getGhostedColImport()

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedColImport2()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Import>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedColImport2(
  int i) const
{
  using Teuchos::rcp;
  if (not useColGidProviders_)
    return getGhostedImport2(i);
  if (colImporters2_[i].is_null())
    colImporters2_[i] =
      rcp(new Epetra_Import(*getGhostedColMap2(i), *getColMap(i)));
  return colImporters2_[i];
} // end of getGhostedColImport2()

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedExport()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Export>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedExport(
  int i) const
{
  using Teuchos::rcp;
  if (exporters_[i].is_null())
    exporters_[i] = rcp(new Epetra_Export(*getGhostedMap(i), *getMap(i)));
  return exporters_[i];
} // end of getGhostedExport()

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedExport2()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Export>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedExport2(
  int i) const
{
  using Teuchos::rcp;
  if (exporters_[i].is_null())
    exporters_[i] = rcp(new Epetra_Export(*getGhostedMap2(i), *getMap(i)));
  return exporters_[i];
} // end of getGhostedExport2()

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedColExport()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Export>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedColExport(
  int i) const
{
  using Teuchos::rcp;
  if (not useColGidProviders_)
    return getGhostedExport(i);
  if (colExporters_[i].is_null())
    colExporters_[i] = rcp(new Epetra_Export(*getGhostedColMap(i),
      *getColMap(i)));
  return colExporters_[i];
} // end of getGhostedColExport()

///////////////////////////////////////////////////////////////////////////////
//
//  getGhostedColExport2()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Export>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
getGhostedColExport2(
  int i) const
{
  using Teuchos::rcp;
  if (not useColGidProviders_)
    return getGhostedExport2(i);
  if (colExporters_[i].is_null())
    colExporters_[i] = rcp(new Epetra_Export(*getGhostedColMap2(i),
      *getColMap(i)));
  return colExporters_[i];
} // end of getGhostedColExport2()

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
buildMap(int i) const
{
   std::vector<int> indices;

   // get the global indices
   getGlobalIndexer(i)->getOwnedIndicesAsInt(indices);

   return Teuchos::rcp(new Epetra_Map(-1,indices.size(),&indices[0],0,*eComm_));
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
buildColMap(int i) const
{
   if(not useColGidProviders_)
     return buildMap(i);

   std::vector<int> indices;

   // get the global indices
   getColGlobalIndexer(i)->getOwnedIndicesAsInt(indices);

   return Teuchos::rcp(new Epetra_Map(-1,indices.size(),&indices[0],0,*eComm_));
}

///////////////////////////////////////////////////////////////////////////////
//
//  buildGhostedMap()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
buildGhostedMap(
  int i) const
{
  using std::vector;
  using Teuchos::rcp;
  vector<int> indices;
  getGlobalIndexer(i)->getOwnedAndGhostedIndicesAsInt(indices);
  return rcp(new Epetra_Map(-1, indices.size(), &indices[0], 0, *eComm_));
} // end of buildGhostedMap()

///////////////////////////////////////////////////////////////////////////////
//
//  buildGhostedMap2()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
buildGhostedMap2(
  int i) const
{
  using std::vector;
  using Teuchos::rcp;
  vector<int> indices;
  getGlobalIndexer(i)->getGhostedIndicesAsInt(indices);
  return rcp(new Epetra_Map(-1, indices.size(), &indices[0], 0, *eComm_));
} // end of buildGhostedMap2()

///////////////////////////////////////////////////////////////////////////////
//
//  buildColGhostedMap()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
buildColGhostedMap(
  int i) const
{
  using std::vector;
  using Teuchos::rcp;
  if (not useColGidProviders_)
    return buildGhostedMap(i);
  vector<int> indices;
  getColGlobalIndexer(i)->getOwnedAndGhostedIndicesAsInt(indices);
  return rcp(new Epetra_Map(-1, indices.size(), &indices[0], 0, *eComm_));
} // end of buildColGhostedMap()

///////////////////////////////////////////////////////////////////////////////
//
//  buildColGhostedMap2()
//
///////////////////////////////////////////////////////////////////////////////
template<typename Traits, typename LocalOrdinalT>
const Teuchos::RCP<Epetra_Map>
BlockedEpetraLinearObjFactory<Traits, LocalOrdinalT>::
buildColGhostedMap2(
  int i) const
{
  using std::vector;
  using Teuchos::rcp;
  if (not useColGidProviders_)
    return buildGhostedMap2(i);
  vector<int> indices;
  getColGlobalIndexer(i)->getGhostedIndicesAsInt(indices);
  return rcp(new Epetra_Map(-1, indices.size(), &indices[0], 0, *eComm_));
} // end of buildColGhostedMap2()

// get the graph of the crs matrix
template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
buildGraph(int i,int j) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;

   // build the map and allocate the space for the graph and
   // grab the ghosted graph
   RCP<Epetra_Map> map_i = getMap(i);
   RCP<Epetra_Map> map_j = getColMap(j);

   TEUCHOS_ASSERT(map_i!=Teuchos::null);
   TEUCHOS_ASSERT(map_j!=Teuchos::null);

   RCP<Epetra_CrsGraph> graph  = rcp(new Epetra_CrsGraph(Copy,*map_i,0));
   RCP<Epetra_CrsGraph> oGraph = buildFilteredGhostedGraph(i,j);
     // this is the only place buildFilteredGhostedGraph is called. That is because
     // only the unghosted graph should reflect any of the filtering.

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
buildGhostedGraph(int i,int j,bool optimizeStorage) const
{
   // build the map and allocate the space for the graph
   Teuchos::RCP<Epetra_Map> rowMap = getGhostedMap(i);
   Teuchos::RCP<Epetra_Map> colMap = getGhostedColMap(j);
   Teuchos::RCP<Epetra_CrsGraph> graph = Teuchos::rcp(new Epetra_CrsGraph(Copy,*rowMap,*colMap,0));

   std::vector<std::string> elementBlockIds;
   
   Teuchos::RCP<const GlobalIndexer> rowProvider, colProvider;
 
   rowProvider = getGlobalIndexer(i);
   colProvider = getColGlobalIndexer(j);

   rowProvider->getElementBlockIds(elementBlockIds); // each sub provider "should" have the
                                                     // same element blocks
                                                        
   const Teuchos::RCP<const ConnManager> conn_mgr = colProvider->getConnManager();
   const bool han = conn_mgr.is_null() ? false : conn_mgr->hasAssociatedNeighbors();

   // graph information about the mesh
   std::vector<std::string>::const_iterator blockItr;
   for(blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;

      // grab elements for this block
      const std::vector<LocalOrdinalT> & elements = rowProvider->getElementBlock(blockId); // each sub provider "should" have the
                                                                                           // same elements in each element block

      // get information about number of indicies
      std::vector<int> row_gids;
      std::vector<int> col_gids;

      // loop over the elemnts
      for(std::size_t elmt=0;elmt<elements.size();elmt++) {
         rowProvider->getElementGIDsAsInt(elements[elmt],row_gids);
         colProvider->getElementGIDsAsInt(elements[elmt],col_gids);

         if (han) {
           const std::vector<LocalOrdinalT>& aes = conn_mgr->getAssociatedNeighbors(elements[elmt]);
           for (typename std::vector<LocalOrdinalT>::const_iterator eit = aes.begin();
                eit != aes.end(); ++eit) {
             std::vector<int> other_col_gids;
             colProvider->getElementGIDsAsInt(*eit, other_col_gids);
             col_gids.insert(col_gids.end(), other_col_gids.begin(), other_col_gids.end());
           }
         }

         for(std::size_t row=0;row<row_gids.size();row++)
            graph->InsertGlobalIndices(row_gids[row],col_gids.size(),&col_gids[0]);
      }
   }

   // finish filling the graph: Make sure the colmap and row maps coincide to 
   //                           minimize calls to LID lookups
   graph->FillComplete(*colMap,*rowMap);
   if(optimizeStorage)
     graph->OptimizeStorage();

   return graph;
}

template <typename Traits,typename LocalOrdinalT>
const Teuchos::RCP<Epetra_CrsGraph> BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
buildFilteredGhostedGraph(int i,int j) const
{
   using Teuchos::RCP;
   using Teuchos::rcp;
   using Teuchos::rcp_dynamic_cast;

   // figure out if the domain is filtered
   RCP<const Filtered_GlobalIndexer> filtered_ugi 
       = rcp_dynamic_cast<const Filtered_GlobalIndexer>(getColGlobalIndexer(j));

   // domain is unfiltered, a filtered graph is just the original ghosted graph
   if(filtered_ugi==Teuchos::null)
     return buildGhostedGraph(i,j,true);

   // get all local indices that are active (i.e. unfiltered)
   std::vector<int> ghostedActive;
   filtered_ugi->getOwnedAndGhostedNotFilteredIndicator(ghostedActive);

   // This will build a new ghosted graph without optimized storage so entries can be removed.
   Teuchos::RCP<Epetra_CrsGraph> filteredGraph = buildGhostedGraph(i,j,false); 
       // false implies that storage is not optimzied 

   // remove filtered column entries
   for(int k=0;k<filteredGraph->NumMyRows();++k) {
     std::vector<int> removedIndices;
     int numIndices = 0;
     int * indices = 0;
     TEUCHOS_ASSERT(filteredGraph->ExtractMyRowView(k,numIndices,indices)==0);

     for(int m=0;m<numIndices;++m) {
       if(ghostedActive[indices[m]]==0)
         removedIndices.push_back(indices[m]);
     }

     TEUCHOS_ASSERT(filteredGraph->RemoveMyIndices(k,Teuchos::as<int>(removedIndices.size()),&removedIndices[0])==0);
   }

   // finish filling the graph
   Teuchos::RCP<Epetra_Map> rowMap = getGhostedMap(i);
   Teuchos::RCP<Epetra_Map> colMap = getGhostedColMap(j);

   TEUCHOS_ASSERT(filteredGraph->FillComplete(*colMap,*rowMap)==0);
   TEUCHOS_ASSERT(filteredGraph->OptimizeStorage()==0);

   return filteredGraph;
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
   return eComm_;
}

template <typename Traits,typename LocalOrdinalT>
int BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getBlockRowCount() const
{
  return rowDOFManagerContainer_->getFieldBlocks();
}

template <typename Traits,typename LocalOrdinalT>
int BlockedEpetraLinearObjFactory<Traits,LocalOrdinalT>::
getBlockColCount() const
{
  return colDOFManagerContainer_->getFieldBlocks();
}

}

#endif // __Panzer_BlockedEpetraLinearObjFactory_impl_hpp__
