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

#ifndef PANZER_SCATTER_DIRICHLET_RESIDUAL_BLOCEDEPETRA_IMPL_HPP
#define PANZER_SCATTER_DIRICHLET_RESIDUAL_BLOCEDEPETRA_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Epetra_Map.h"
#include "Epetra_Vector.h"
#include "Epetra_CrsMatrix.h"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BlockedTpetraLinearObjContainer.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Thyra_SpmdVectorBase.hpp"
#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_get_Epetra_Operator.hpp"

#include "Teuchos_FancyOStream.hpp"

#include <unordered_map>

template <typename EvalT,typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterDirichletResidual_BlockedTpetra<EvalT,TRAITS,LO,GO,NodeT>::
ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & /* indexer */,
                                       const Teuchos::ParameterList& p)
{
  std::string scatterName = p.get<std::string>("Scatter Name");
  Teuchos::RCP<PHX::FieldTag> scatterHolder =
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  Teuchos::RCP<PHX::DataLayout> dl =
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis")->functional;

  // build the vector of fields that this is dependent on
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    PHX::MDField<const ScalarT,Cell,NODE> scatterField = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterField.fieldTag());
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder);

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
// Specialization: Residual
// **********************************************************************


template <typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer,
                                const Teuchos::ParameterList& p)
   : globalIndexer_(indexer)
   , globalDataKey_("Residual Scatter Container")
{
  std::string scatterName = p.get<std::string>("Scatter Name");
  scatterHolder_ =
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  // grab map from evaluated names to field names
  fieldMap_ = p.get< Teuchos::RCP< std::map<std::string,std::string> > >("Dependent Map");

  // determine if we are scattering an initial condition
  scatterIC_ = p.isParameter("Scatter Initial Condition") ? p.get<bool>("Scatter Initial Condition") : false;

  Teuchos::RCP<PHX::DataLayout> dl = (!scatterIC_) ?
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis")->functional :
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional;
  if (!scatterIC_) {
    side_subcell_dim_ = p.get<int>("Side Subcell Dimension");
    local_side_id_ = p.get<int>("Local Side ID");
  }

  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  checkApplyBC_ = p.isParameter("Check Apply BC") ? p.get<bool>("Check Apply BC") : false;
  applyBC_.resize(names.size()); // must allocate (even if not used) to support lambda capture
  if (checkApplyBC_) {
    for (std::size_t eq = 0; eq < names.size(); ++eq) {
      applyBC_[eq] = PHX::MDField<const bool,Cell,NODE>(std::string("APPLY_BC_")+fieldMap_->find(names[eq])->second,dl);
      this->addDependentField(applyBC_[eq]);
    }
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  const Workset & workset_0 = (*d.worksets_)[0];
  const std::string blockId = this->wda(workset_0).block_id;

  fieldIds_.resize(scatterFields_.size());
  fieldOffsets_.resize(scatterFields_.size());
  basisIndexForMDFieldOffsets_.resize(scatterFields_.size());
  fieldGlobalIndexers_.resize(scatterFields_.size());
  productVectorBlockIndex_.resize(scatterFields_.size());
  int maxElementBlockGIDCount = -1;
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {
    // get field ID from DOF manager
    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;

    const int globalFieldNum = globalIndexer_->getFieldNum(fieldName); // Field number in the aggregate BlockDOFManager
    productVectorBlockIndex_[fd] = globalIndexer_->getFieldBlock(globalFieldNum);
    fieldGlobalIndexers_[fd] = globalIndexer_->getFieldDOFManagers()[productVectorBlockIndex_[fd]];
    fieldIds_[fd] = fieldGlobalIndexers_[fd]->getFieldNum(fieldName); // Field number in the sub-global-indexer

    // Offsets and basisIndex depend on whether scattering IC or Dirichlet BC
    if (!scatterIC_) {
      const auto& offsetPair = fieldGlobalIndexers_[fd]->getGIDFieldOffsets_closure(blockId,fieldIds_[fd],side_subcell_dim_,local_side_id_);
      {
        const auto& offsets =  offsetPair.first;
        fieldOffsets_[fd] = Kokkos::View<int*,PHX::Device>("ScatterDirichletResidual_BlockedTpetra(Residual):fieldOffsets",offsets.size());
        auto hostOffsets = Kokkos::create_mirror_view(fieldOffsets_[fd]);
        for (std::size_t i=0; i < offsets.size(); ++i)
          hostOffsets(i) = offsets[i];
        Kokkos::deep_copy(fieldOffsets_[fd], hostOffsets);
      }
      {
        const auto& basisIndex =  offsetPair.second;
        basisIndexForMDFieldOffsets_[fd] = Kokkos::View<int*,PHX::Device>("ScatterDirichletResidual_BlockedTpetra(Residual):basisIndexForMDFieldOffsets",basisIndex.size());
        auto hostBasisIndex = Kokkos::create_mirror_view(basisIndexForMDFieldOffsets_[fd]);
        for (std::size_t i=0; i < basisIndex.size(); ++i)
          hostBasisIndex(i) = basisIndex[i];
        Kokkos::deep_copy(basisIndexForMDFieldOffsets_[fd], hostBasisIndex);
      }
    }
    else {
      // For ICs, only need offsets, not basisIndex
      const std::vector<int>& offsets = fieldGlobalIndexers_[fd]->getGIDFieldOffsets(blockId,fieldIds_[fd]);
      fieldOffsets_[fd] = Kokkos::View<int*,PHX::Device>("ScatterDirichletResidual_BlockedTpetra(Residual):fieldOffsets",offsets.size());
      auto hostOffsets = Kokkos::create_mirror_view(fieldOffsets_[fd]);
      for (std::size_t i=0; i < offsets.size(); ++i)
        hostOffsets(i) = offsets[i];
      Kokkos::deep_copy(fieldOffsets_[fd], hostOffsets);
    }

    maxElementBlockGIDCount = std::max(fieldGlobalIndexers_[fd]->getElementBlockGIDCount(blockId),maxElementBlockGIDCount);
    typename PHX::Device().fence();
  }

  // We will use one workset lid view for all fields, but has to be
  // sized big enough to hold the largest elementBlockGIDCount in the
  // ProductVector.
  worksetLIDs_ = Kokkos::View<LO**,PHX::Device>("ScatterResidual_BlockedTpetra(Residual):worksetLIDs",
                                                scatterFields_[0].extent(0),
                                                maxElementBlockGIDCount);
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   // extract dirichlet counter from container
   Teuchos::RCP<const ContainerType> blockContainer
         = Teuchos::rcp_dynamic_cast<ContainerType>(d.gedc->getDataObject("Dirichlet Counter"),true);

   dirichletCounter_ = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockContainer->get_f(),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));

   // extract linear object container
   blockedContainer_ = Teuchos::rcp_dynamic_cast<const ContainerType>(d.gedc->getDataObject(globalDataKey_),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(blockedContainer_));
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;
   using Thyra::VectorBase;
   using Thyra::ProductVectorBase;

   const auto& localCellIds = this->wda(workset).cell_local_ids_k;

   RCP<ProductVectorBase<double> > thyraScatterTarget = (!scatterIC_) ?
     rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer_->get_f(),true) :
     rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer_->get_x(),true);

   // Loop over scattered fields
   int currentWorksetLIDSubBlock = -1;
   for (std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
     // workset LIDs only change for different sub blocks
     if (productVectorBlockIndex_[fieldIndex] != currentWorksetLIDSubBlock) {
       fieldGlobalIndexers_[fieldIndex]->getElementLIDs(localCellIds,worksetLIDs_);
       currentWorksetLIDSubBlock = productVectorBlockIndex_[fieldIndex];
     }

     // Get Scatter target block
     const auto& tpetraScatterTarget = *((rcp_dynamic_cast<Thyra::TpetraVector<RealType,LO,GO,NodeT>>(thyraScatterTarget->getNonconstVectorBlock(productVectorBlockIndex_[fieldIndex]),true))->getTpetraVector());
     const auto& kokkosScatterTarget = tpetraScatterTarget.template getLocalView<PHX::mem_space>();

     // Get dirichlet counter block
     const auto& tpetraDirichletCounter = *((rcp_dynamic_cast<Thyra::TpetraVector<RealType,LO,GO,NodeT>>(dirichletCounter_->getNonconstVectorBlock(productVectorBlockIndex_[fieldIndex]),true))->getTpetraVector());
     const auto& kokkosDirichletCounter = tpetraDirichletCounter.template getLocalView<PHX::mem_space>();

     // Class data fields for lambda capture
     const auto fieldOffsets = fieldOffsets_[fieldIndex];
     const auto basisIndices = basisIndexForMDFieldOffsets_[fieldIndex];
     const auto worksetLIDs = worksetLIDs_;
     const auto fieldValues = scatterFields_[fieldIndex];
     const auto applyBC = applyBC_[fieldIndex].get_static_view();
     const bool checkApplyBC = checkApplyBC_;

     if (!scatterIC_) {

       Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device>(0,workset.num_cells), KOKKOS_LAMBDA (const int& cell) {
         for (int basis=0; basis < static_cast<int>(fieldOffsets.size()); ++basis) {
           const int lid = worksetLIDs(cell,fieldOffsets(basis));
           if (lid < 0) // not on this processor!
             continue;
           const int basisIndex = basisIndices(basis);

           // Possible warp divergence for hierarchic
           if (checkApplyBC)
             if (!applyBC(cell,basisIndex))
               continue;

           kokkosScatterTarget(lid,0) = fieldValues(cell,basisIndex);
           kokkosDirichletCounter(lid,0) = 1.0;
         }
       });

     } else {

       Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device>(0,workset.num_cells), KOKKOS_LAMBDA (const int& cell) {
         for (int basis=0; basis < static_cast<int>(fieldOffsets.size()); ++basis) {
           const int lid = worksetLIDs(cell,fieldOffsets(basis));
           if (lid < 0) // not on this processor!
             continue;
           kokkosScatterTarget(lid,0) = fieldValues(cell,basis);
           kokkosDirichletCounter(lid,0) = 1.0;
         }
       });

     }
   }

}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template <typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
ScatterDirichletResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer,
                                const Teuchos::ParameterList& p)
   : globalIndexer_(indexer)
   , globalDataKey_("Residual Scatter Container")
{
  std::string scatterName = p.get<std::string>("Scatter Name");
  scatterHolder_ =
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  // grab map from evaluated names to field names
  fieldMap_ = p.get< Teuchos::RCP< std::map<std::string,std::string> > >("Dependent Map");

  Teuchos::RCP<PHX::DataLayout> dl =
    p.get< Teuchos::RCP<panzer::PureBasis> >("Basis")->functional;

  side_subcell_dim_ = p.get<int>("Side Subcell Dimension");
  local_side_id_ = p.get<int>("Local Side ID");

  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  checkApplyBC_ = p.get<bool>("Check Apply BC");
  applyBC_.resize(names.size()); // must allocate (even if not used) to support lambda capture
  if (checkApplyBC_) {
    for (std::size_t eq = 0; eq < names.size(); ++eq) {
      applyBC_[eq] = PHX::MDField<const bool,Cell,NODE>(std::string("APPLY_BC_")+fieldMap_->find(names[eq])->second,dl);
      this->addDependentField(applyBC_[eq]);
    }
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual (Jacobian)");
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d,
                      PHX::FieldManager<TRAITS>& /* fm */)
{
  const Workset & workset_0 = (*d.worksets_)[0];
  const std::string blockId = this->wda(workset_0).block_id;

  fieldIds_.resize(scatterFields_.size());
  fieldOffsets_.resize(scatterFields_.size());
  basisIndexForMDFieldOffsets_.resize(scatterFields_.size());
  productVectorBlockIndex_.resize(scatterFields_.size());
  for(std::size_t fd=0;fd<scatterFields_.size();++fd) {

    std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    const int globalFieldNum = globalIndexer_->getFieldNum(fieldName); // Field number in the aggregate BlockDOFManager
    productVectorBlockIndex_[fd] = globalIndexer_->getFieldBlock(globalFieldNum);
    const auto& fieldGlobalIndexer = globalIndexer_->getFieldDOFManagers()[productVectorBlockIndex_[fd]];
    fieldIds_[fd] = fieldGlobalIndexer->getFieldNum(fieldName); // Field number in the sub-global-indexer

    const auto& offsetPair = fieldGlobalIndexer->getGIDFieldOffsets_closure(blockId,fieldIds_[fd],side_subcell_dim_,local_side_id_);
    {
      const auto& offsets =  offsetPair.first;
      fieldOffsets_[fd] = Kokkos::View<int*,PHX::Device>("ScatterDirichletResidual_BlockedTpetra(Jacobian):fieldOffsets",offsets.size());
      auto hostOffsets = Kokkos::create_mirror_view(fieldOffsets_[fd]);
      for (std::size_t i=0; i < offsets.size(); ++i)
        hostOffsets(i) = offsets[i];
      Kokkos::deep_copy(fieldOffsets_[fd], hostOffsets);
    }
    {
      const auto& basisIndex =  offsetPair.second;
      basisIndexForMDFieldOffsets_[fd] = Kokkos::View<int*,PHX::Device>("ScatterDirichletResidual_BlockedTpetra(Jacobian):basisIndexForMDFieldOffsets",basisIndex.size());
      auto hostBasisIndex = Kokkos::create_mirror_view(basisIndexForMDFieldOffsets_[fd]);
      for (std::size_t i=0; i < basisIndex.size(); ++i)
        hostBasisIndex(i) = basisIndex[i];
      Kokkos::deep_copy(basisIndexForMDFieldOffsets_[fd], hostBasisIndex);
    }
  }

  // This is sized differently than the Residual implementation since
  // we need the LIDs for all sub-blocks, not just the single
  // sub-block for the field residual scatter.
  int elementBlockGIDCount = 0;
  for (const auto blockDOFMgr : globalIndexer_->getFieldDOFManagers())
    elementBlockGIDCount += blockDOFMgr->getElementBlockGIDCount(blockId);

  worksetLIDs_ = Kokkos::View<LO**,PHX::Device>("ScatterDirichletResidual_BlockedTpetra(Jacobian):worksetLIDs",
                                                scatterFields_[0].extent(0),
                                                elementBlockGIDCount);

  // Compute the block offsets
  const auto& blockGlobalIndexers = globalIndexer_->getFieldDOFManagers();
  const int numBlocks = static_cast<int>(globalIndexer_->getFieldDOFManagers().size());
  blockOffsets_ = Kokkos::View<LO*,PHX::Device>("ScatterDirichletResidual_BlockedTpetra(Jacobian):blockOffsets_",
                                                numBlocks+1); // Number of fields, plus a sentinel
  const auto hostBlockOffsets = Kokkos::create_mirror_view(blockOffsets_);
  for (int blk=0;blk<numBlocks;blk++) {
    int blockOffset = globalIndexer_->getBlockGIDOffset(blockId,blk);
    hostBlockOffsets(blk) = blockOffset;
  }
  blockOffsets_(numBlocks) = blockOffsets_(numBlocks-1) + blockGlobalIndexers[blockGlobalIndexers.size()-1]->getElementBlockGIDCount(blockId);
  Kokkos::deep_copy(blockOffsets_,hostBlockOffsets);

  // Make sure the that hard coded derivative dimension in the
  // evaluate call is large enough to hold all derivatives for each
  // sub block load
  for (int blk=0;blk<numBlocks;blk++) {
    const int blockDerivativeSize = hostBlockOffsets(blk+1) - hostBlockOffsets(blk);
    TEUCHOS_TEST_FOR_EXCEPTION(blockDerivativeSize > maxDerivativeArraySize_, std::runtime_error,
                               "ERROR: the derivative dimension for sub block "
                               << blk << "with a value of " << blockDerivativeSize
                               << "is larger than the size allocated for cLIDs and vals "
                               << "in the evaluate call! You must manually increase the "
                               << "size and recompile!");
  }

  typename PHX::Device().fence();
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   // extract dirichlet counter from container
   Teuchos::RCP<const ContainerType> blockContainer
         = Teuchos::rcp_dynamic_cast<const ContainerType>(d.gedc->getDataObject("Dirichlet Counter"),true);

   dirichletCounter_ = Teuchos::rcp_dynamic_cast<Thyra::ProductVectorBase<double> >(blockContainer->get_f(),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(dirichletCounter_));

   // extract linear object container
   blockedContainer_ = Teuchos::rcp_dynamic_cast<const ContainerType>(d.gedc->getDataObject(globalDataKey_),true);
   TEUCHOS_ASSERT(!Teuchos::is_null(blockedContainer_));
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterDirichletResidual_BlockedTpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::VectorBase;
  using Thyra::ProductVectorBase;
  using Thyra::BlockedLinearOpBase;

  const auto& localCellIds = this->wda(workset).cell_local_ids_k;

  const int numFieldBlocks = globalIndexer_->getNumFieldBlocks();
  const RCP<const ContainerType> blockedContainer = blockedContainer_;
  const RCP<ProductVectorBase<double>> thyraBlockResidual = rcp_dynamic_cast<ProductVectorBase<double>>(blockedContainer_->get_f());
  const bool haveResidual = Teuchos::nonnull(thyraBlockResidual);
  const RCP<BlockedLinearOpBase<double>> Jac = rcp_dynamic_cast<BlockedLinearOpBase<double>>(blockedContainer_->get_A(),true);

  // Get the local data for the sub-block crs matrices. First allocate
  // on host and then deep_copy to device. The sub-blocks are
  // unmanaged since they are allocated and ref counted separately on
  // host.
  using LocalMatrixType = KokkosSparse::CrsMatrix<double,LO,PHX::Device,Kokkos::MemoryTraits<Kokkos::Unmanaged>, size_t>;
  typename Kokkos::View<LocalMatrixType**,PHX::Device>::HostMirror
    hostJacTpetraBlocks("panzer::ScatterResidual_BlockTpetra<Jacobian>::hostJacTpetraBlocks", numFieldBlocks,numFieldBlocks);

  Kokkos::View<int**,PHX::Device> blockExistsInJac =   Kokkos::View<int**,PHX::Device>("blockExistsInJac_",numFieldBlocks,numFieldBlocks);
  auto hostBlockExistsInJac = Kokkos::create_mirror_view(blockExistsInJac);

  for (int row=0; row < numFieldBlocks; ++row) {
    for (int col=0; col < numFieldBlocks; ++col) {
      const auto thyraTpetraOperator = rcp_dynamic_cast<Thyra::TpetraLinearOp<double,LO,GO,NodeT>>(Jac->getNonconstBlock(row,col),false);
      if (nonnull(thyraTpetraOperator)) {

        // HACK to enforce views in the CrsGraph to be
        // Unmanaged. Passing in the MemoryTrait<Unmanaged> doesn't
        // work as the CrsGraph in the CrsMatrix ignores the
        // MemoryTrait. Need to use the runtime constructor by passing
        // in points to ensure Unmanaged.  See:
        // https://github.com/kokkos/kokkos/issues/1581

        // These two lines are the original code we can revert to when #1581 is fixed.
        // const auto crsMatrix = rcp_dynamic_cast<Tpetra::CrsMatrix<double,LO,GO,NodeT>>(thyraTpetraOperator->getTpetraOperator(),true);
        // new (&hostJacTpetraBlocks(row,col)) KokkosSparse::CrsMatrix<double,LO,PHX::Device,Kokkos::MemoryTraits<Kokkos::Unmanaged>> (crsMatrix->getLocalMatrix());

        // Instead do this
        {
          // Grab the local managed matrix and graph
          const auto tpetraCrsMatrix = rcp_dynamic_cast<Tpetra::CrsMatrix<double,LO,GO,NodeT>>(thyraTpetraOperator->getTpetraOperator(),true);
          const auto managedMatrix = tpetraCrsMatrix->getLocalMatrix();
          const auto managedGraph = managedMatrix.graph;

          // Create runtime unmanaged versions
          using StaticCrsGraphType = typename LocalMatrixType::StaticCrsGraphType;
          StaticCrsGraphType unmanagedGraph;
          unmanagedGraph.entries = typename StaticCrsGraphType::entries_type(managedGraph.entries.data(),managedGraph.entries.extent(0));
          unmanagedGraph.row_map = typename StaticCrsGraphType::row_map_type(managedGraph.row_map.data(),managedGraph.row_map.extent(0));
          unmanagedGraph.row_block_offsets = typename StaticCrsGraphType::row_block_type(managedGraph.row_block_offsets.data(),managedGraph.row_block_offsets.extent(0));

          typename LocalMatrixType::values_type unmanagedValues(managedMatrix.values.data(),managedMatrix.values.extent(0));
          LocalMatrixType unmanagedMatrix(managedMatrix.values.label(), managedMatrix.numCols(), unmanagedValues, unmanagedGraph);
          new (&hostJacTpetraBlocks(row,col)) LocalMatrixType(unmanagedMatrix);
        }

        hostBlockExistsInJac(row,col) = 1;
      }
      else {
        hostBlockExistsInJac(row,col) = 0;
      }
    }
  }
  typename Kokkos::View<LocalMatrixType**,PHX::Device>
    jacTpetraBlocks("panzer::ScatterResidual_BlockedTpetra<Jacobian>::jacTpetraBlocks",numFieldBlocks,numFieldBlocks);
  Kokkos::deep_copy(jacTpetraBlocks,hostJacTpetraBlocks);
  Kokkos::deep_copy(blockExistsInJac,hostBlockExistsInJac);
  typename PHX::Device().fence();

  // worksetLIDs is larger for Jacobian than Residual fill. Need the
  // entire set of field offsets for derivative indexing no matter
  // which block row you are scattering. The residual only needs the
  // lids for the sub-block that it is scattering to. The subviews
  // below are to offset the LID blocks correctly.
  const auto& globalIndexers = globalIndexer_->getFieldDOFManagers();
  for (size_t block=0; block < globalIndexers.size(); ++block) {
    const auto subviewOfBlockLIDs = Kokkos::subview(worksetLIDs_,Kokkos::ALL(), std::make_pair(blockOffsets_(block),blockOffsets_(block+1)));
    globalIndexers[block]->getElementLIDs(localCellIds,subviewOfBlockLIDs);
  }

  // Loop over scattered fields
  for (std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {

    const int blockRowIndex = productVectorBlockIndex_[fieldIndex];
    typename Tpetra::Vector<double,LO,GO,PHX::Device>::dual_view_type::t_dev kokkosResidual;
    if (haveResidual) {
      const auto& tpetraResidual = *((rcp_dynamic_cast<Thyra::TpetraVector<RealType,LO,GO,NodeT>>(thyraBlockResidual->getNonconstVectorBlock(blockRowIndex),true))->getTpetraVector());
      kokkosResidual = tpetraResidual.template getLocalView<PHX::mem_space>();
    }

    // Get dirichlet counter block
    const auto& tpetraDirichletCounter = *((rcp_dynamic_cast<Thyra::TpetraVector<RealType,LO,GO,NodeT>>(dirichletCounter_->getNonconstVectorBlock(blockRowIndex),true))->getTpetraVector());
    const auto& kokkosDirichletCounter = tpetraDirichletCounter.template getLocalView<PHX::mem_space>();

    // Class data fields for lambda capture
    const Kokkos::View<const int*,PHX::Device> fieldOffsets = fieldOffsets_[fieldIndex];
    const auto& basisIndices = basisIndexForMDFieldOffsets_[fieldIndex];
    const Kokkos::View<const LO**,PHX::Device> worksetLIDs = worksetLIDs_;
    const PHX::View<const ScalarT**> fieldValues = scatterFields_[fieldIndex].get_static_view();
    const Kokkos::View<const LO*,PHX::Device> blockOffsets = blockOffsets_;
    const auto& applyBC = applyBC_[fieldIndex].get_static_view();
    const bool checkApplyBC = checkApplyBC_;

    Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device>(0,workset.num_cells), KOKKOS_LAMBDA (const int& cell) {
      LO cLIDs[maxDerivativeArraySize_];
      typename Sacado::ScalarType<ScalarT>::type vals[maxDerivativeArraySize_];

      for(int basis=0; basis < static_cast<int>(fieldOffsets.size()); ++basis) {
        const int rowLID = worksetLIDs(cell,fieldOffsets(basis));

        if (rowLID < 0) // not on this processor!
          continue;

        const int basisIndex = basisIndices(basis);

        // Possible warp divergence for hierarchic
        if (checkApplyBC)
          if (!applyBC(cell,basisIndex))
            continue;

        typedef PHX::MDField<const ScalarT,Cell,NODE> FieldType;
        typename FieldType::array_type::reference_type tmpFieldVal = fieldValues(cell,basisIndex);

        if (haveResidual)
          kokkosResidual(rowLID,0) = tmpFieldVal.val();

        kokkosDirichletCounter(rowLID,0) = 1.0;

        // Zero out entire matrix row
        for (int blockColIndex=0; blockColIndex < numFieldBlocks; ++blockColIndex) {
          if (blockExistsInJac(blockRowIndex,blockColIndex)) {
            const auto& rowEntries = jacTpetraBlocks(blockRowIndex,blockColIndex).row(rowLID);
            for (int i=0; i < rowEntries.length; ++i)
              rowEntries.value(i) = Traits::RealType(0.0);
          }
        }

        // Set values
        for (int blockColIndex=0; blockColIndex < numFieldBlocks; ++blockColIndex) {
          if (blockExistsInJac(blockRowIndex,blockColIndex)) {
            const int start = blockOffsets(blockColIndex);
            const int stop = blockOffsets(blockColIndex+1);
            const int sensSize = stop-start;
            // Views may be padded. Use contiguous arrays here
            for (int i=0; i < sensSize; ++i) {
              cLIDs[i] = worksetLIDs(cell,start+i);
              vals[i] = tmpFieldVal.fastAccessDx(start+i);
            }
            jacTpetraBlocks(blockRowIndex,blockColIndex).replaceValues(rowLID,cLIDs,sensSize,vals,true,true);
          }
        }
      }
    });

  }

  // Placement delete on view of matrices
  for (int row=0; row < numFieldBlocks; ++row) {
    for (int col=0; col < numFieldBlocks; ++col) {
      if (hostBlockExistsInJac(row,col)) {
        hostJacTpetraBlocks(row,col).~CrsMatrix();
      }
    }
  }

}

// **********************************************************************

#endif
