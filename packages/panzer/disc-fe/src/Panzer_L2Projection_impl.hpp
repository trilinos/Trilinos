// @HEADER
// @HEADER

#ifndef PANZER_L2_PROJECTION_IMPL_HPP
#define PANZER_L2_PROJECTION_IMPL_HPP

#include "Teuchos_DefaultMpiComm.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Panzer_BasisDescriptor.hpp"
#include "Panzer_TpetraLinearObjContainer.hpp"
#include "Panzer_TpetraLinearObjFactory.hpp"
#include "Panzer_BlockedTpetraLinearObjFactory.hpp"
#include "Panzer_ElementBlockIdToPhysicsIdMap.hpp"
#include "Panzer_DOFManagerFactory.hpp"
#include "Panzer_BlockedDOFManagerFactory.hpp"
#include "Panzer_IntrepidFieldPattern.hpp"
#include "Panzer_WorksetDescriptor.hpp"
#include "Panzer_WorksetContainer.hpp"
#include "Panzer_Workset.hpp"

namespace panzer {

  template<typename LO, typename GO>
  void panzer::L2Projection<LO,GO>::setup(const panzer::BasisDescriptor& targetBasis,
                                          const panzer::IntegrationDescriptor& integrationDescriptor,
                                          const Teuchos::RCP<const Teuchos::MpiComm<int>>& comm,
                                          const Teuchos::RCP<const panzer::ConnManager<LO,GO>>& connManager,
                                          const std::vector<std::string>& elementBlockNames,
                                          const Teuchos::RCP<panzer::WorksetContainer> worksetContainer)
  {
    targetBasisDescriptor_ = targetBasis;
    integrationDescriptor_ = integrationDescriptor;
    comm_ = comm;
    connManager_ = connManager;
    elementBlockNames_ = elementBlockNames;
    worksetContainer_ = worksetContainer;
    setupCalled_ = true;

    // Build target DOF Manager
    targetGlobalIndexer_ =
      Teuchos::rcp(new panzer::DOFManager<LO,GO>(Teuchos::rcp_const_cast<panzer::ConnManager<LO,GO>>(connManager),*(comm->getRawMpiComm())));

    // For hybrid mesh, blocks could have different topology
    for (const auto& eBlock : elementBlockNames_) {
      std::vector<shards::CellTopology> topologies;
      connManager_->getElementBlockTopologies(topologies);
      std::vector<std::string> ebNames;
      connManager_->getElementBlockIds(ebNames);
      const auto search = std::find(ebNames.cbegin(),ebNames.cend(),eBlock);
      TEUCHOS_ASSERT(search != ebNames.cend());
      const int index = std::distance(ebNames.cbegin(),search);
      const auto& cellTopology = topologies[index];

      auto intrepidBasis = panzer::createIntrepid2Basis<PHX::Device,double,double>(targetBasisDescriptor_.getType(),
                                                                                   targetBasisDescriptor_.getOrder(),
                                                                                   cellTopology);
      Teuchos::RCP<const panzer::FieldPattern> fieldPattern(new panzer::Intrepid2FieldPattern(intrepidBasis));
      // Field name is the basis type
      targetGlobalIndexer_->addField(eBlock,targetBasisDescriptor_.getType(),fieldPattern);
    }

    targetGlobalIndexer_->buildGlobalUnknowns();

    // Check workset needs are correct
  }

  template<typename LO, typename GO>
  Teuchos::RCP<Tpetra::CrsMatrix<double,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>
  panzer::L2Projection<LO,GO>::buildMassMatrix()
  {
    TEUCHOS_ASSERT(setupCalled_);

    // Allocate the owned matrix
    std::vector<Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO>>> indexers;
    indexers.push_back(targetGlobalIndexer_);

    panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,LO,GO,panzer::TpetraNodeType> factory(comm_,indexers);

    auto ownedMatrix = factory.getTpetraMatrix(0,0);
    auto ghostedMatrix = factory.getGhostedTpetraMatrix(0,0);
    ownedMatrix->resumeFill();
    ownedMatrix->setAllToScalar(0.0);
    ghostedMatrix->resumeFill();
    ghostedMatrix->setAllToScalar(0.0);
    PHX::Device::fence();

    auto M = ghostedMatrix->getLocalMatrix();
    const int fieldIndex = targetGlobalIndexer_->getFieldNum(targetBasisDescriptor_.getType());

    // Loop over element blocks and fill mass matrix
    for (const auto& block : elementBlockNames_) {

      // Based on descriptor, currently assumes there should only be one workset
      panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,true);
      const auto& worksets = worksetContainer_->getWorksets(wd);

      for (const auto& workset : *worksets) {

        const auto& basisValues = workset.getBasisValues(targetBasisDescriptor_,integrationDescriptor_);
        
        const auto& unweightedBasis = basisValues.basis_scalar;
        const auto& weightedBasis = basisValues.weighted_basis_scalar;
        
        // Offsets (this assumes UVM, need to fix)
        const std::vector<LO>& offsets = targetGlobalIndexer_->getGIDFieldOffsets(block,fieldIndex);
        PHX::View<LO*> kOffsets("MassMatrix: Offsets",offsets.size());
        for(const auto& i : offsets)
          kOffsets(i) = offsets[i];
        
        PHX::Device::fence();
        
        // Local Ids
        Kokkos::View<LO**,PHX::Device> localIds("MassMatrix: LocalIds", workset.numOwnedCells()+workset.numGhostCells()+workset.numVirtualCells(),
                                                targetGlobalIndexer_->getElementBlockGIDCount(block));
        
        // Remove the ghosted cell ids or the call to getElementLocalIds will spill array bounds
        const auto cellLocalIdsNoGhost = Kokkos::subview(workset.cell_local_ids_k,std::make_pair(0,workset.numOwnedCells()));
        
        targetGlobalIndexer_->getElementLIDs(cellLocalIdsNoGhost,localIds);
        
        const int numBasisPoints = static_cast<int>(weightedBasis.extent(1));
        Kokkos::parallel_for(workset.numOwnedCells(),KOKKOS_LAMBDA (const int& cell) {
            
          LO cLIDs[256];
          const int numIds = static_cast<int>(localIds.extent(1));
          for(int i=0;i<numIds;++i)
            cLIDs[i] = localIds(cell,i);
          
          double vals[256];
          const int numQP = static_cast<int>(unweightedBasis.extent(2));
          
          for (int qp=0; qp < numQP; ++qp) {
            for (int row=0; row < numBasisPoints; ++row) {
              int offset = kOffsets(row);
              LO lid = localIds(cell,offset);
              
              for (int col=0; col < numIds; ++col)
                vals[col] = unweightedBasis(cell,row,qp) * weightedBasis(cell,col,qp);
              
              M.sumIntoValues(lid,cLIDs,numIds,vals,true,true);
            }
          }
   
        });
      }
    }
    PHX::exec_space::fence();

    ghostedMatrix->fillComplete();
    const auto exporter = factory.getGhostedExport(0);
    ownedMatrix->doExport(*ghostedMatrix, *exporter, Tpetra::ADD);
    ownedMatrix->fillComplete();

    return ownedMatrix;
  }

  template<typename LO, typename GO>
  Teuchos::RCP<Tpetra::MultiVector<double,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>
  panzer::L2Projection<LO,GO>::buildInverseLumpedMassMatrix()
  {
    using Teuchos::rcp;
    const auto massMatrix = this->buildMassMatrix();
    const auto lumpedMassMatrix = rcp(new Tpetra::MultiVector<double,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>(massMatrix->getDomainMap(),1,true));
    const auto tmp = rcp(new Tpetra::MultiVector<double,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>(massMatrix->getRangeMap(),1,false));
    tmp->putScalar(1.0);
    massMatrix->apply(*tmp,*lumpedMassMatrix);
    lumpedMassMatrix->reciprocal(*lumpedMassMatrix);
    return lumpedMassMatrix;
  }

  template<typename LO, typename GO>
  Teuchos::RCP<Tpetra::CrsMatrix<double,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>
  panzer::L2Projection<LO,GO>::buildRHSMatrix(const panzer::UniqueGlobalIndexer<LO,GO>& sourceGlobalIndexer,
                                              const Teuchos::RCP<const Tpetra::Map<LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>& inputOwnedSourceMap,
                                              const std::string& sourceFieldName,
                                              const panzer::BasisDescriptor& sourceBasisDescriptor,
                                              const int directionIndex)
  {
    // *******************
    // Create Retangular matrix (both ghosted and owned).
    // *******************
    using Teuchos::RCP;
    using Teuchos::rcp;
    using MapType = Tpetra::Map<LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>;
    using GraphType = Tpetra::CrsGraph<LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>;
    using ExportType = Tpetra::Export<LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>;
    using MatrixType = Tpetra::CrsMatrix<double,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>;

    // *******************
    // Build ghosted graph
    // *******************

    RCP<MapType> ghostedTargetMap;
    {
      std::vector<GO> indices;
      targetGlobalIndexer_->getOwnedAndGhostedIndices(indices);
      ghostedTargetMap = rcp(new MapType(Teuchos::OrdinalTraits<GO>::invalid(),indices,0,comm_));
    }

    RCP<MapType> ghostedSourceMap;
    {
      std::vector<GO> indices;
      sourceGlobalIndexer.getOwnedAndGhostedIndices(indices);
      ghostedSourceMap = rcp(new MapType(Teuchos::OrdinalTraits<GO>::invalid(),indices,0,comm_));
    }

    RCP<GraphType> ghostedGraph = rcp(new GraphType(ghostedTargetMap,ghostedSourceMap,0));

    // Now insert the non-zero pattern per row
    std::vector<std::string> elementBlockIds;
    targetGlobalIndexer_->getElementBlockIds(elementBlockIds);
    std::vector<std::string>::const_iterator blockItr;
    for (blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;
      const std::vector<LO> & elements = targetGlobalIndexer_->getElementBlock(blockId);

      std::vector<GO> row_gids;
      std::vector<GO> col_gids;

      for(std::size_t elmt=0;elmt<elements.size();elmt++) {
        targetGlobalIndexer_->getElementGIDs(elements[elmt],row_gids);
        sourceGlobalIndexer.getElementGIDs(elements[elmt],col_gids);
        for(std::size_t row=0;row<row_gids.size();row++)
          ghostedGraph->insertGlobalIndices(row_gids[row],col_gids);
      }
    }

    // Fill complete with range and domain map
    ghostedGraph->fillComplete(ghostedSourceMap,ghostedTargetMap);

    // *****************
    // Build owned graph
    // *****************
    RCP<MapType> ownedTargetMap;
    {
      std::vector<GO> indices;
      targetGlobalIndexer_->getOwnedIndices(indices);
      ownedTargetMap = rcp(new MapType(Teuchos::OrdinalTraits<GO>::invalid(),indices,0,comm_));
    }

    RCP<const MapType> ownedSourceMap = inputOwnedSourceMap;
    if (is_null(ownedSourceMap)) {
      std::vector<GO> indices;
      sourceGlobalIndexer.getOwnedIndices(indices);
      ownedSourceMap = rcp(new MapType(Teuchos::OrdinalTraits<GO>::invalid(),indices,0,comm_));
    }

    RCP<GraphType> ownedGraph = rcp(new GraphType(ownedTargetMap,0));
    RCP<const ExportType> exporter = rcp(new ExportType(ghostedTargetMap,ownedTargetMap));
    ownedGraph->doExport(*ghostedGraph, *exporter, Tpetra::INSERT);
    ownedGraph->fillComplete(ownedSourceMap,ownedTargetMap);

    RCP<MatrixType> ghostedMatrix = rcp(new MatrixType(ghostedGraph));
    RCP<MatrixType> ownedMatrix = rcp(new MatrixType(ownedGraph));
    // ghostedMatrix.fillComplete();
    // ghostedMatrix.resumeFill();

    ghostedMatrix->setAllToScalar(0.0);
    ownedMatrix->setAllToScalar(0.0);
    PHX::Device::fence();

    // *******************
    // Fill ghosted matrix
    // *******************
    for (const auto& block : elementBlockNames_) {

      panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,true);
      const auto& worksets = worksetContainer_->getWorksets(wd);
      for (const auto& workset : *worksets) {

        // Get target basis values: current implementation assumes target basis is HGrad
        const auto& targetBasisValues = workset.getBasisValues(targetBasisDescriptor_,integrationDescriptor_);
        const auto& targetWeightedBasis = targetBasisValues.weighted_basis_scalar.get_static_view();

        // Sources can be any basis
        const auto& sourceBasisValues = workset.getBasisValues(sourceBasisDescriptor,integrationDescriptor_);
        Kokkos::View<const double***,PHX::Device> sourceUnweightedScalarBasis;
        Kokkos::View<const double****,PHX::Device> sourceUnweightedVectorBasis;
        bool useRankThreeBasis = false; // default to gradient or vector basis
        if ( (sourceBasisDescriptor.getType() == "HGrad") || (sourceBasisDescriptor.getType() == "Const") ) {
          if (directionIndex == -1) { // Project dof value
            sourceUnweightedScalarBasis = sourceBasisValues.basis_scalar.get_static_view();
            useRankThreeBasis = true;
          }
          else { // Project dof gradient of scalar basis
            sourceUnweightedVectorBasis = sourceBasisValues.grad_basis.get_static_view();
          }
        }
        else { // Project vector value
          sourceUnweightedVectorBasis = sourceBasisValues.basis_vector.get_static_view();
        }

        // Get the element local ids
        Kokkos::View<LO**,PHX::Device> targetLocalIds("buildRHSMatrix: targetLocalIds", workset.numOwnedCells(),
                                                      targetGlobalIndexer_->getElementBlockGIDCount(block));
        Kokkos::View<LO**,PHX::Device> sourceLocalIds("buildRHSMatrix: sourceLocalIds", workset.numOwnedCells(),
                                                      sourceGlobalIndexer.getElementBlockGIDCount(block));
        {
          // Remove the ghosted cell ids or the call to getElementLocalIds will spill array bounds
          const auto cellLocalIdsNoGhost = Kokkos::subview(workset.cell_local_ids_k,std::make_pair(0,workset.numOwnedCells()));
          targetGlobalIndexer_->getElementLIDs(cellLocalIdsNoGhost,targetLocalIds);
          sourceGlobalIndexer.getElementLIDs(cellLocalIdsNoGhost,sourceLocalIds);
        }

        // Get the offsets
        Kokkos::View<LO*,PHX::Device> targetFieldOffsets;
        {
          const auto fieldIndex = targetGlobalIndexer_->getFieldNum(targetBasisDescriptor_.getType());
          const std::vector<LO>& offsets = targetGlobalIndexer_->getGIDFieldOffsets(block,fieldIndex);
          targetFieldOffsets = Kokkos::View<LO*,PHX::Device>("L2Projection:buildRHS:targetFieldOffsets",offsets.size());
          const auto hostOffsets = Kokkos::create_mirror_view(targetFieldOffsets);
          for(size_t i=0; i < offsets.size(); ++i)
            hostOffsets(i) = offsets[i];
          Kokkos::deep_copy(targetFieldOffsets,hostOffsets);
          PHX::Device::fence();
        }

        Kokkos::View<LO*,PHX::Device> sourceFieldOffsets;
        {
          const auto fieldIndex = sourceGlobalIndexer.getFieldNum(sourceFieldName);
          const std::vector<LO>& offsets = sourceGlobalIndexer.getGIDFieldOffsets(block,fieldIndex);
          sourceFieldOffsets = Kokkos::View<LO*,PHX::Device>("L2Projection:buildRHS:sourceFieldOffsets",offsets.size());
          const auto hostOffsets = Kokkos::create_mirror_view(sourceFieldOffsets);
          for(size_t i=0; i <offsets.size(); ++i)
            hostOffsets(i) = offsets[i];
          Kokkos::deep_copy(sourceFieldOffsets,hostOffsets);
          PHX::Device::fence();
        }

        const auto localMatrix = ghostedMatrix->getLocalMatrix();
        const int numRows = static_cast<int>(targetWeightedBasis.extent(1));
        int tmpNumCols = -1;
        int tmpNumQP = -1;
        if (useRankThreeBasis) {
          tmpNumCols = static_cast<int>(sourceUnweightedScalarBasis.extent(1));
          tmpNumQP = static_cast<int>(sourceUnweightedScalarBasis.extent(2));
        }
        else {
          tmpNumCols = static_cast<int>(sourceUnweightedVectorBasis.extent(1));
          tmpNumQP = static_cast<int>(sourceUnweightedVectorBasis.extent(2));
        }
        const int numCols = tmpNumCols;
        const int numQP = tmpNumQP;

        Kokkos::parallel_for(Kokkos::RangePolicy<PHX::Device>(0,workset.numOwnedCells()),KOKKOS_LAMBDA (const int& cell) {
          LO cLIDs[256];
          double vals[256];
          for (int row = 0; row < numRows; ++row) {
            const int rowOffset = targetFieldOffsets(row);
            const int rowLID = targetLocalIds(cell,rowOffset);
            for (int col = 0; col < numCols; ++col)
              vals[col] = 0.0;

            for (int col = 0; col < numCols; ++col) {
              for (int qp = 0; qp < numQP; ++qp) {
                const int colOffset = sourceFieldOffsets(col);
                const int colLID = sourceLocalIds(cell,colOffset);
                cLIDs[col] = colLID;
                if (useRankThreeBasis)
                  vals[col] += sourceUnweightedScalarBasis(cell,col,qp) * targetWeightedBasis(cell,row,qp);
                else
                  vals[col] += sourceUnweightedVectorBasis(cell,col,qp,directionIndex) * targetWeightedBasis(cell,row,qp);
              }
            }
            localMatrix.sumIntoValues(rowLID,cLIDs,numCols,vals,true,true);
          }
        });
      }

    }
    ghostedMatrix->fillComplete(ghostedSourceMap,ghostedTargetMap);
    ownedMatrix->resumeFill();
    ownedMatrix->doExport(*ghostedMatrix,*exporter,Tpetra::ADD);
    ownedMatrix->fillComplete(ownedSourceMap,ownedTargetMap);

    return ownedMatrix;
  }

}

#endif
