// @HEADER
// @HEADER

#ifndef PANZER_L2_PROJECTION_IMPL_HPP
#define PANZER_L2_PROJECTION_IMPL_HPP

#include "Teuchos_DefaultMpiComm.hpp"
#include "Tpetra_CrsGraph.hpp"
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
                                          const Teuchos::RCP<Teuchos::MpiComm<int>>& comm,
                                          const Teuchos::RCP<const panzer::ConnManager<LO,GO>>& connManager,
                                          const std::vector<std::string>& elementBlockNames,
                                          const Teuchos::RCP<panzer::WorksetContainer>& worksetContainer)
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
      panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,false);
      const auto& worksets = worksetContainer_->getWorksets(wd);
      TEUCHOS_ASSERT(worksets->size() == 1);
      const auto& workset = (*worksets)[0];
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

      std::cout << "workset.owned=" << workset.numOwnedCells() << ", ghosted=" << workset.numGhostCells() << ", virtual=" << workset.numVirtualCells() << std::endl;

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
    PHX::exec_space::fence();

    ghostedMatrix->fillComplete();
    const auto exporter = factory.getGhostedExport(0);
    ownedMatrix->doExport(*ghostedMatrix, *exporter, Tpetra::ADD);
    ownedMatrix->fillComplete();

    return ownedMatrix;
  }

  template<typename LO, typename GO>
  Teuchos::RCP<Tpetra::CrsMatrix<double,LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>
  panzer::L2Projection<LO,GO>::buildRHSMatrix(const Teuchos::RCP<panzer::DOFManager<LO,GO>>& sourceGlobalIndexer,
                                              const Teuchos::RCP<Tpetra::Map<LO,GO,Kokkos::Compat::KokkosDeviceWrapperNode<PHX::Device>>>& inputOwnedSourceMap,
                                              const std::string& sourceFieldName,
                                              const panzer::BasisDescriptor& sourceBasisDescriptor,
                                              const bool sourceIsVectorBasis,
                                              const int vectorBasisIndex)
  {
    TEUCHOS_ASSERT(nonnull(sourceGlobalIndexer));

    /*
      NOTE: This will be the minimal memory implementation. It is not
      yet complete. In the meantime, will use a memory hog
      implementaiton.

    TEUCHOS_ASSERT(setupCalled_);

    // Create local maps with lid->gid lookups
    std::vector<Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO>>> indexers;
    indexers.push_back(targetGlobalIndexer_);
    panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,LO,GO,panzer::TpetraNodeType> targetFactory(comm_,indexers);
    auto ghostedTargetMap = targetFactory.getGhostedMap(0);
    auto localTargetMap = ghostedTargetMap->getLocalMap();

    indexers.clear();
    indexers.push_back(sourceGlobalIndexer);
    panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,LO,GO,panzer::TpetraNodeType> sourceFactory(comm_,indexers);
    auto ghostedSourceMap = sourceFactory.getGhostedMap(0);
    auto localSourceMap = ghostedSourceMap->getLocalMap();
    
    // Count the unknown entries per row on the ghosted map, indexing is storing GO in std::pair<Row,Col>
    Kokkos::DualView<LO*, PHX::Device> entriesPerRow("L2Projection:entriesPerRow",ghostedTargetMap->getNodeNumElements());

    // Stores row LID and col GID pairs to generate graph
    Kokkos::UnorderedMap<Kokkos::pair<LO,GO>,void,PHX::mem_space> globalNodeSet;
    {
      size_t setCapacity = (28ull * entriesPerRow.extent(0));
      unsigned failedInsertCount = 0;
        
      do {
        Kokkos::deep_copy(entriesPerRow.template view<PHX::mem_space>(),0);
        globalNodeSet = Kokkos::UnorderedMap<Kokkos::pair<LO,GO>,void,PHX::mem_space>( ( setCapacity += failedInsertCount ) );
        setCapacity = globalNodeSet.capacity(); // Actaul size may be larger than requested:
        failedInsertCount = 0;

        for (const auto& block : elementBlockNames_) {

          // Offsets (this assumes UVM, need to fix)
          const int targetFieldIndex = targetGlobalIndexer_->getFieldNum(targetBasisDescriptor_.getType());
          const std::vector<LO>& targetHostOffsets = targetGlobalIndexer_->getGIDFieldOffsets(block,targetFieldIndex);
          PHX::View<LO*> targetFieldOffsets("L2Projection:buildRHSMatrix: targetFieldOffsets",targetHostOffsets.size());
          for(const auto& i : targetHostOffsets)
            targetFieldOffsets(i) = targetHostOffsets[i];

          // Offsets (this assumes UVM, need to fix)
          const int sourceFieldIndex = sourceGlobalIndexer->getFieldNum(sourceFieldName);
          const std::vector<LO>& sourceHostOffsets = sourceGlobalIndexer->getGIDFieldOffsets(block,sourceFieldIndex);
          PHX::View<LO*> sourceFieldOffsets("L2Projection:buildRHSMatrix: sourceFieldOffsets",sourceHostOffsets.size());
          for(const auto& i : sourceHostOffsets)
            sourceFieldOffsets(i) = sourceHostOffsets[i];

          PHX::Device::fence();
          
          panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,false);
          const auto& worksets = worksetContainer_->getWorksets(wd);
          for (const auto& workset : *worksets) {
            
            const auto& targetBasisValues = workset.getBasisValues(targetBasisDescriptor_,integrationDescriptor_);
            const auto& sourceBasisValues = workset.getBasisValues(sourceBasisDescriptor,integrationDescriptor_);

            // Local Ids
            Kokkos::View<LO**,PHX::Device> targetLocalIds("buildRHSMatrix: targetLocalIds", workset.numOwnedCells(),
                                                          targetGlobalIndexer_->getElementBlockGIDCount(block));        
            // Remove the ghosted cell ids or the call to getElementLocalIds will spill array bounds
            const auto targetCellLocalIdsNoGhost = Kokkos::subview(workset.cell_local_ids_k,std::make_pair(0,workset.numOwnedCells()));
            targetGlobalIndexer_->getElementLIDs(targetCellLocalIdsNoGhost,targetLocalIds);

            Kokkos::View<LO**,PHX::Device> sourceLocalIds("buildRHSMatrix: sourceLocalIds", workset.numOwnedCells(),
                                                          sourceGlobalIndexer->getElementBlockGIDCount(block));        
            // Remove the ghosted cell ids or the call to getElementLocalIds will spill array bounds
            const auto sourceCellLocalIdsNoGhost = Kokkos::subview(workset.cell_local_ids_k,std::make_pair(0,workset.numOwnedCells()));
            sourceGlobalIndexer->getElementLIDs(sourceCellLocalIdsNoGhost,sourceLocalIds);

            // Launch parallel kernel
            Kokkos::parallel_reduce(Kokkos::RangePolicy<PHX::exec_space>(0,10), KOKKOS_LAMBDA (const int& cell, unsigned& failedInsertCount)
            {
              for (int tb = 0; tb < static_cast<int>(targetBasis.extent(1)); ++tb) {
                LO targetOffset = targetOffsets(tb); 
                LO targetLID = targetLocalIds(cell,targetFieldOffset);
                
                for (int sb = 0; sb < static_cast<int>(sourceBasis.extent(1)); ++sb) {
                  LO sourceOffset = sourceOffsets(sb); 
                  LO sourceLID = sourceLocalIds(cell,sourceFieldOffset);
                  GO sourceGID = localSourceMap.getGlobalElement(sourceLID);

                  const auto result = global_node_set_.insert( Kokkos::make_pair(targetLID,sourceGID) );
                  
                  if (!result.failed()) {
                    Kokkos::atomic_fetch_add(&row_count_( row_gid * num_equations_ + row_eq_offset ), num_equations_);
                  }
                  else {
                    ++failed_count;
                  }
                }
              }
            }, failedInsertCount);

            PHX::exec_space::fence();
          }
        }  
      } while (failedInsertCount); 
    }

    // Maybe don't need these?
    entriesPerRow.template modify<PHX::mem_space>();
    entriesPerRow.template sync<typename PHX::mem_space::HostSpace>();
    */

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
      sourceGlobalIndexer->getOwnedAndGhostedIndices(indices);
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
        sourceGlobalIndexer->getElementGIDs(elements[elmt],col_gids);
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
 
    RCP<MapType> ownedSourceMap = inputOwnedSourceMap;
    if (is_null(ownedSourceMap)) {
      std::vector<GO> indices;
      sourceGlobalIndexer->getOwnedIndices(indices);
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

      panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,false);
      const auto& worksets = worksetContainer_->getWorksets(wd);
      for (const auto& workset : *worksets) {
        const auto& targetBasisValues = workset.getBasisValues(targetBasisDescriptor_,integrationDescriptor_);
        const auto& sourceBasisValues = workset.getBasisValues(sourceBasisDescriptor,integrationDescriptor_);
        
        // Assumes target is in HDiv
        const auto& targetUnweightedBasis = targetBasisValues.basis_scalar.get_static_view();
        const auto& targetWeightedBasis = targetBasisValues.weighted_basis_scalar.get_static_view();

        // One component of vector
        Kokkos::View<const double***,PHX::Device> sourceUnweightedScalarBasis;
        Kokkos::View<const double***,PHX::Device> sourceWeightedScalarBasis;
        Kokkos::View<const double****,PHX::Device> sourceUnweightedVectorBasis;
        Kokkos::View<const double****,PHX::Device> sourceWeightedVectorBasis;

        if (sourceIsVectorBasis) {
          sourceUnweightedVectorBasis = sourceBasisValues.basis_vector.get_static_view();
          sourceWeightedVectorBasis = sourceBasisValues.weighted_basis_vector.get_static_view();
        }
        else {
          sourceUnweightedScalarBasis = sourceBasisValues.basis_scalar.get_static_view();
          sourceWeightedScalarBasis = sourceBasisValues.weighted_basis_scalar.get_static_view();
        }

        // Get the local ids
        Kokkos::View<LO**,PHX::Device> targetLocalIds("buildRHSMatrix: targetLocalIds", workset.numOwnedCells(),
                                                      targetGlobalIndexer_->getElementBlockGIDCount(block));
        Kokkos::View<LO**,PHX::Device> sourceLocalIds("buildRHSMatrix: sourceLocalIds", workset.numOwnedCells(),
                                                      sourceGlobalIndexer->getElementBlockGIDCount(block));
        {
          // Remove the ghosted cell ids or the call to getElementLocalIds will spill array bounds
          const auto cellLocalIdsNoGhost = Kokkos::subview(workset.cell_local_ids_k,std::make_pair(0,workset.numOwnedCells()));
          targetGlobalIndexer_->getElementLIDs(cellLocalIdsNoGhost,targetLocalIds);
          sourceGlobalIndexer->getElementLIDs(cellLocalIdsNoGhost,sourceLocalIds);
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
          const auto fieldIndex = sourceGlobalIndexer->getFieldNum(sourceFieldName);
          const std::vector<LO>& offsets = sourceGlobalIndexer->getGIDFieldOffsets(block,fieldIndex);
          sourceFieldOffsets = Kokkos::View<LO*,PHX::Device>("L2Projection:buildRHS:sourceFieldOffsets",offsets.size());
          const auto hostOffsets = Kokkos::create_mirror_view(sourceFieldOffsets);
          for(size_t i=0; i <offsets.size(); ++i)
            hostOffsets(i) = offsets[i];
          Kokkos::deep_copy(sourceFieldOffsets,hostOffsets);
          PHX::Device::fence();
        }

        const auto localMatrix = ghostedMatrix->getLocalMatrix();
        const int numRows = static_cast<int>(targetWeightedBasis.extent(1));
        const int numCols = static_cast<int>(sourceWeightedVectorBasis.extent(1));
        const int numQP = static_cast<int>(sourceWeightedVectorBasis.extent(2));

        std::cout << "ROGER numRows=" << numRows << ", numCols=" << numCols << ", numQP=" << numQP << std::endl;


        if (sourceIsVectorBasis) {
          // loop over cells
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
                  vals[col] += sourceUnweightedVectorBasis(cell,col,qp,vectorBasisIndex)*targetWeightedBasis(cell,row,qp);
                }
              }
              localMatrix.sumIntoValues(rowLID,cLIDs,numCols,vals,true,true);
            }
          });
        }
        else {
          // loop over cells
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
                  vals[col] += sourceUnweightedScalarBasis(cell,col,qp)*targetWeightedBasis(cell,row,qp);
                }
              }
              localMatrix.sumIntoValues(rowLID,cLIDs,numCols,vals,true,true);
            }
          });
          
        }
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
