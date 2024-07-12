// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_SCATTER_RESIDUAL_BLOCKEDEPETRA_IMPL_HPP
#define PANZER_SCATTER_RESIDUAL_BLOCKEDEPETRA_IMPL_HPP

#include "Teuchos_RCP.hpp"
#include "Teuchos_Assert.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_BlockedDOFManager.hpp"
#include "Panzer_PureBasis.hpp"
#include "Panzer_BlockedTpetraLinearObjContainer.hpp"
#include "Panzer_LOCPair_GlobalEvaluationData.hpp"
#include "Panzer_HashUtils.hpp"
#include "Panzer_GlobalEvaluationDataContainer.hpp"

#include "Thyra_ProductVectorBase.hpp"
#include "Thyra_BlockedLinearOpBase.hpp"
#include "Thyra_TpetraVector.hpp"
#include "Thyra_TpetraLinearOp.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "KokkosSparse_CrsMatrix.hpp"

#include "Phalanx_DataLayout_MDALayout.hpp"

#include "Teuchos_FancyOStream.hpp"

template <typename EvalT,typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterResidual_BlockedTpetra<EvalT,TRAITS,LO,GO,NodeT>::
ScatterResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & /* indexer */,
                              const Teuchos::ParameterList& p)
{
  std::string scatterName = p.get<std::string>("Scatter Name");
  Teuchos::RCP<PHX::FieldTag> scatterHolder =
    Teuchos::rcp(new PHX::Tag<ScalarT>(scatterName,Teuchos::rcp(new PHX::MDALayout<Dummy>(0))));

  // get names to be evaluated
  const std::vector<std::string>& names =
    *(p.get< Teuchos::RCP< std::vector<std::string> > >("Dependent Names"));

  Teuchos::RCP<PHX::DataLayout> dl =
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional;

  // build the vector of fields that this is dependent on
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    PHX::MDField<const ScalarT,Cell,NODE> field = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(field.fieldTag());
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder);

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
// Specialization: Residual
// **********************************************************************

template <typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
ScatterResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer,
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
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional;

  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual");
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Residual,TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d,
		      PHX::FieldManager<TRAITS>& /* fm */)
{
  const Workset & workset_0 = (*d.worksets_)[0];
  const std::string blockId = this->wda(workset_0).block_id;

  fieldIds_.resize(scatterFields_.size());
  fieldOffsets_.resize(scatterFields_.size());
  fieldGlobalIndexers_.resize(scatterFields_.size());
  productVectorBlockIndex_.resize(scatterFields_.size());
  int maxElementBlockGIDCount = -1;
  for(std::size_t fd=0; fd < scatterFields_.size(); ++fd) {
    const std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    const int globalFieldNum = globalIndexer_->getFieldNum(fieldName); // Field number in the aggregate BlockDOFManager
    productVectorBlockIndex_[fd] = globalIndexer_->getFieldBlock(globalFieldNum);
    fieldGlobalIndexers_[fd] = globalIndexer_->getFieldDOFManagers()[productVectorBlockIndex_[fd]];
    fieldIds_[fd] = fieldGlobalIndexers_[fd]->getFieldNum(fieldName); // Field number in the sub-global-indexer

    const std::vector<int>& offsets = fieldGlobalIndexers_[fd]->getGIDFieldOffsets(blockId,fieldIds_[fd]);
    fieldOffsets_[fd] = PHX::View<int*>("ScatterResidual_BlockedTpetra(Residual):fieldOffsets",offsets.size());
    auto hostOffsets = Kokkos::create_mirror_view(fieldOffsets_[fd]);
    for (std::size_t i=0; i < offsets.size(); ++i)
      hostOffsets(i) = offsets[i];
    Kokkos::deep_copy(fieldOffsets_[fd], hostOffsets);

    maxElementBlockGIDCount = std::max(fieldGlobalIndexers_[fd]->getElementBlockGIDCount(blockId),maxElementBlockGIDCount);
  }

  // We will use one workset lid view for all fields, but has to be
  // sized big enough to hold the largest elementBlockGIDCount in the
  // ProductVector.
  worksetLIDs_ = PHX::View<LO**>("ScatterResidual_BlockedTpetra(Residual):worksetLIDs",
                                                scatterFields_[0].extent(0),
						maxElementBlockGIDCount);
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Residual, TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   // extract linear object container
   blockedContainer_ = rcp_dynamic_cast<const ContainerType>(d.gedc->getDataObject(globalDataKey_));

   if(blockedContainer_==Teuchos::null) {
     RCP<const LOCPair_GlobalEvaluationData> gdata = rcp_dynamic_cast<const LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true);
     blockedContainer_ = rcp_dynamic_cast<const ContainerType>(gdata->getGhostedLOC());
   }
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Residual,TRAITS,LO,GO,NodeT>::
evaluateFields(typename TRAITS::EvalData workset)
{
  using Teuchos::RCP;
  using Teuchos::rcp_dynamic_cast;
  using Thyra::VectorBase;
  using Thyra::ProductVectorBase;

  const auto& localCellIds = this->wda(workset).cell_local_ids_k;
  const RCP<ProductVectorBase<double>> thyraBlockResidual = rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer_->get_f(),true);

  // Loop over scattered fields
  int currentWorksetLIDSubBlock = -1;
  for (std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {
    // workset LIDs only change for different sub blocks
    if (productVectorBlockIndex_[fieldIndex] != currentWorksetLIDSubBlock) {
      fieldGlobalIndexers_[fieldIndex]->getElementLIDs(localCellIds,worksetLIDs_);
      currentWorksetLIDSubBlock = productVectorBlockIndex_[fieldIndex];
    }

    auto& tpetraResidual = *((rcp_dynamic_cast<Thyra::TpetraVector<RealType,LO,GO,NodeT>>(thyraBlockResidual->getNonconstVectorBlock(productVectorBlockIndex_[fieldIndex]),true))->getTpetraVector());
    const auto& kokkosResidual = tpetraResidual.getLocalViewDevice(Tpetra::Access::ReadWrite);

    // Class data fields for lambda capture
    const auto& fieldOffsets = fieldOffsets_[fieldIndex];
    const auto& worksetLIDs = worksetLIDs_;
    const auto& fieldValues = scatterFields_[fieldIndex].get_static_view();

    Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device>(0,workset.num_cells), KOKKOS_LAMBDA (const int& cell) {
      for(int basis=0; basis < static_cast<int>(fieldOffsets.size()); ++basis) {
	      const int lid = worksetLIDs(cell,fieldOffsets(basis));
	      Kokkos::atomic_add(&kokkosResidual(lid,0), fieldValues(cell,basis));
      }
    });
  }
}

// **********************************************************************
// Specialization: Jacobian
// **********************************************************************

template <typename TRAITS,typename LO,typename GO,typename NodeT>
panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Jacobian, TRAITS,LO,GO,NodeT>::
ScatterResidual_BlockedTpetra(const Teuchos::RCP<const BlockedDOFManager> & indexer,
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
    p.get< Teuchos::RCP<const panzer::PureBasis> >("Basis")->functional;

  // build the vector of fields that this is dependent on
  scatterFields_.resize(names.size());
  for (std::size_t eq = 0; eq < names.size(); ++eq) {
    scatterFields_[eq] = PHX::MDField<const ScalarT,Cell,NODE>(names[eq],dl);

    // tell the field manager that we depend on this field
    this->addDependentField(scatterFields_[eq]);
  }

  // this is what this evaluator provides
  this->addEvaluatedField(*scatterHolder_);

  if (p.isType<std::string>("Global Data Key"))
     globalDataKey_ = p.get<std::string>("Global Data Key");

  this->setName(scatterName+" Scatter Residual (Jacobian)");
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Jacobian,TRAITS,LO,GO,NodeT>::
postRegistrationSetup(typename TRAITS::SetupData d,
		      PHX::FieldManager<TRAITS>& /* fm */)
{
  const Workset & workset_0 = (*d.worksets_)[0];
  const std::string blockId = this->wda(workset_0).block_id;

  fieldIds_.resize(scatterFields_.size());
  fieldOffsets_.resize(scatterFields_.size());
  productVectorBlockIndex_.resize(scatterFields_.size());
  for (std::size_t fd=0; fd < scatterFields_.size(); ++fd) {
    const std::string fieldName = fieldMap_->find(scatterFields_[fd].fieldTag().name())->second;
    const int globalFieldNum = globalIndexer_->getFieldNum(fieldName); // Field number in the aggregate BlockDOFManager
    productVectorBlockIndex_[fd] = globalIndexer_->getFieldBlock(globalFieldNum);
    const auto& fieldGlobalIndexer = globalIndexer_->getFieldDOFManagers()[productVectorBlockIndex_[fd]];
    fieldIds_[fd] = fieldGlobalIndexer->getFieldNum(fieldName); // Field number in the sub-global-indexer

    const std::vector<int>& offsets = globalIndexer_->getGIDFieldOffsets(blockId,globalFieldNum);
    fieldOffsets_[fd] = PHX::View<int*>("ScatterResidual_BlockedTpetra(Jacobian):fieldOffsets",offsets.size());
    auto hostOffsets = Kokkos::create_mirror_view(fieldOffsets_[fd]);
    for (std::size_t i=0; i < offsets.size(); ++i)
      hostOffsets(i) = offsets[i];
    Kokkos::deep_copy(fieldOffsets_[fd], hostOffsets);
  }

  // This is sized differently than the Residual implementation since
  // we need the LIDs for all sub-blocks, not just the single
  // sub-block for the field residual scatter.
  int elementBlockGIDCount = 0;
  for (const auto& blockDOFMgr : globalIndexer_->getFieldDOFManagers())
    elementBlockGIDCount += blockDOFMgr->getElementBlockGIDCount(blockId);

  worksetLIDs_ = Kokkos::View<LO**, Kokkos::LayoutRight, PHX::Device>(
                    "ScatterResidual_BlockedTpetra(Jacobian):worksetLIDs_",
                    scatterFields_[0].extent(0), elementBlockGIDCount );

  // Compute the block offsets
  const auto& blockGlobalIndexers = globalIndexer_->getFieldDOFManagers();
  const int numBlocks = static_cast<int>(globalIndexer_->getFieldDOFManagers().size());
  blockOffsets_ = PHX::View<LO*>("ScatterResidual_BlockedTpetra(Jacobian):blockOffsets_",
                                                numBlocks+1); // Number of fields, plus a sentinel
  const auto hostBlockOffsets = Kokkos::create_mirror_view(blockOffsets_);
  for (int blk=0;blk<numBlocks;blk++) {
    int blockOffset = globalIndexer_->getBlockGIDOffset(blockId,blk);
    hostBlockOffsets(blk) = blockOffset;
  }
  hostBlockOffsets(numBlocks) = hostBlockOffsets(numBlocks-1) + blockGlobalIndexers[blockGlobalIndexers.size()-1]->getElementBlockGIDCount(blockId);
  Kokkos::deep_copy(blockOffsets_,hostBlockOffsets);

  // Make sure the that derivative dimension in the evaluate call is large
  // enough to hold all derivatives for each sub block load
  int max_blockDerivativeSize = 0;
  for (int blk=0;blk<numBlocks;blk++) {
    const int blockDerivativeSize = hostBlockOffsets(blk+1) - hostBlockOffsets(blk);
    if ( blockDerivativeSize > max_blockDerivativeSize )
      max_blockDerivativeSize = blockDerivativeSize;
  }
  workset_vals_ = Kokkos::View<typename Sacado::ScalarType<ScalarT>::type**, Kokkos::LayoutRight, PHX::Device>(
                    "ScatterResidual_BlockedTpetra(Jacobian):workset_vals_",
                    scatterFields_[0].extent(0), max_blockDerivativeSize );
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Jacobian,TRAITS,LO,GO,NodeT>::
preEvaluate(typename TRAITS::PreEvalData d)
{
   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   // extract linear object container
   blockedContainer_ = rcp_dynamic_cast<const ContainerType>(d.gedc->getDataObject(globalDataKey_));

   if(blockedContainer_==Teuchos::null) {
     RCP<const LOCPair_GlobalEvaluationData> gdata = rcp_dynamic_cast<const LOCPair_GlobalEvaluationData>(d.gedc->getDataObject(globalDataKey_),true);
     blockedContainer_ = rcp_dynamic_cast<const ContainerType>(gdata->getGhostedLOC());
   }
}

// **********************************************************************
template <typename TRAITS,typename LO,typename GO,typename NodeT>
void panzer::ScatterResidual_BlockedTpetra<panzer::Traits::Jacobian,TRAITS,LO,GO,NodeT>::
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
  const RCP<ProductVectorBase<double>> thyraBlockResidual = rcp_dynamic_cast<ProductVectorBase<double> >(blockedContainer_->get_f());
  const bool haveResidual = Teuchos::nonnull(thyraBlockResidual);
  const RCP<BlockedLinearOpBase<double>> Jac = rcp_dynamic_cast<BlockedLinearOpBase<double> >(blockedContainer_->get_A(),true);

  // Get the local data for the sub-block crs matrices. First allocate
  // on host and then deep_copy to device. The sub-blocks are
  // unmanaged since they are allocated and ref counted separately on
  // host.
  using LocalMatrixType = KokkosSparse::CrsMatrix<double,LO,PHX::Device,Kokkos::MemoryTraits<Kokkos::Unmanaged>, size_t>;
  typename PHX::View<LocalMatrixType**>::HostMirror
    hostJacTpetraBlocks("panzer::ScatterResidual_BlockTpetra<Jacobian>::hostJacTpetraBlocks", numFieldBlocks,numFieldBlocks);

  PHX::View<int**> blockExistsInJac =   PHX::View<int**>("blockExistsInJac_",numFieldBlocks,numFieldBlocks);
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
          const auto managedMatrix = tpetraCrsMatrix->getLocalMatrixDevice();
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
  typename PHX::View<LocalMatrixType**>
    jacTpetraBlocks("panzer::ScatterResidual_BlockedTpetra<Jacobian>::jacTpetraBlocks",numFieldBlocks,numFieldBlocks);
  Kokkos::deep_copy(jacTpetraBlocks,hostJacTpetraBlocks);
  Kokkos::deep_copy(blockExistsInJac,hostBlockExistsInJac);

  // worksetLIDs is larger for Jacobian than Residual fill. Need the
  // entire set of field offsets for derivative indexing no matter
  // which block row you are scattering. The residual only needs the
  // lids for the sub-block that it is scattering to. The subviews
  // below are to offset the LID blocks correctly.
  const auto& globalIndexers = globalIndexer_->getFieldDOFManagers();
  auto blockOffsets_h = Kokkos::create_mirror_view(blockOffsets_);
  Kokkos::deep_copy(blockOffsets_h, blockOffsets_);
  for (size_t block=0; block < globalIndexers.size(); ++block) {
    const auto subviewOfBlockLIDs = Kokkos::subview(worksetLIDs_,Kokkos::ALL(), std::make_pair(blockOffsets_h(block),blockOffsets_h(block+1)));
    globalIndexers[block]->getElementLIDs(localCellIds,subviewOfBlockLIDs);
  }

  // Loop over scattered fields
  for (std::size_t fieldIndex = 0; fieldIndex < scatterFields_.size(); fieldIndex++) {

    const int blockRowIndex = productVectorBlockIndex_[fieldIndex];
    typename Tpetra::Vector<double,LO,GO,PHX::Device>::dual_view_type::t_dev kokkosResidual;
    if (haveResidual) {
      auto& tpetraResidual = *((rcp_dynamic_cast<Thyra::TpetraVector<RealType,LO,GO,NodeT>>(thyraBlockResidual->getNonconstVectorBlock(blockRowIndex),true))->getTpetraVector());
      kokkosResidual = tpetraResidual.getLocalViewDevice(Tpetra::Access::ReadWrite);
    }

    // Class data fields for lambda capture
    const PHX::View<const int*> fieldOffsets = fieldOffsets_[fieldIndex];
    const PHX::View<const ScalarT**> fieldValues = scatterFields_[fieldIndex].get_static_view();
    const PHX::View<const LO*> blockOffsets = blockOffsets_;

    const Kokkos::View<const LO**, Kokkos::LayoutRight, PHX::Device> worksetLIDs = worksetLIDs_;
    Kokkos::View<typename Sacado::ScalarType<ScalarT>::type**, Kokkos::LayoutRight, PHX::Device> workset_vals = workset_vals_;

    Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device>(0,workset.num_cells), KOKKOS_LAMBDA (const int& cell) {
      for(int basis=0; basis < static_cast<int>(fieldOffsets.size()); ++basis) {
        typedef PHX::MDField<const ScalarT,Cell,NODE> FieldType;
        typename FieldType::array_type::reference_type tmpFieldVal = fieldValues(cell,basis);
        const int rowLID = worksetLIDs(cell,fieldOffsets(basis));

        if (haveResidual)
          Kokkos::atomic_add(&kokkosResidual(rowLID,0), tmpFieldVal.val());

        for (int blockColIndex=0; blockColIndex < numFieldBlocks; ++blockColIndex) {
          if (blockExistsInJac(blockRowIndex,blockColIndex)) {
            const int start = blockOffsets(blockColIndex);
            const int stop = blockOffsets(blockColIndex+1);
            const int sensSize = stop-start;

            for (int i=0; i < sensSize; ++i)
              workset_vals(cell,i) = tmpFieldVal.fastAccessDx(start+i);

            jacTpetraBlocks(blockRowIndex,blockColIndex).sumIntoValues(rowLID, &worksetLIDs(cell,start), sensSize, &workset_vals(cell,0), true,true);
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
