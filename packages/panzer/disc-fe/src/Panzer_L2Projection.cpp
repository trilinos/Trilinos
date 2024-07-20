// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_L2_PROJECTION_IMPL_HPP
#define PANZER_L2_PROJECTION_IMPL_HPP

#include "Teuchos_DefaultMpiComm.hpp"
#include "Tpetra_CrsGraph.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_CrsMatrix.hpp"
#include "Panzer_L2Projection.hpp"
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

  void panzer::L2Projection::setup(const panzer::BasisDescriptor& targetBasis,
                                   const panzer::IntegrationDescriptor& integrationDescriptor,
                                   const Teuchos::RCP<const Teuchos::MpiComm<int>>& comm,
                                   const Teuchos::RCP<const panzer::ConnManager>& connManager,
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
    useUserSuppliedBasisValues_ = false;

    // Build target DOF Manager
    targetGlobalIndexer_ =
      Teuchos::rcp(new panzer::DOFManager(Teuchos::rcp_const_cast<panzer::ConnManager>(connManager),*(comm->getRawMpiComm())));

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

  void panzer::L2Projection::useBasisValues(const std::map<std::string,Teuchos::RCP<panzer::BasisValues2<double>>>& bv)
  {
    useUserSuppliedBasisValues_ = true;
    basisValues_ = bv;

    // Check that the basis and integration descriptor match what was
    // supplied in setup.
    for (const auto& eblock : elementBlockNames_) {
      TEUCHOS_ASSERT(basisValues_[eblock]->getBasisDescriptor()==targetBasisDescriptor_);
    }
  }

  Teuchos::RCP<panzer::GlobalIndexer>
  panzer::L2Projection::getTargetGlobalIndexer() const
  {return targetGlobalIndexer_;}

  Teuchos::RCP<Tpetra::CrsMatrix<double,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>>
  panzer::L2Projection::buildMassMatrix(bool use_lumping,
                                        const std::unordered_map<std::string,double>* elementBlockMultipliers)
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("L2Projection::Build Mass Matrix",TopBuildMassMatrix);
    TEUCHOS_ASSERT(setupCalled_);

    if (elementBlockMultipliers != nullptr) {
      TEUCHOS_ASSERT(elementBlockMultipliers->size() == elementBlockNames_.size());
    }

    // Allocate the owned matrix
    std::vector<Teuchos::RCP<const panzer::GlobalIndexer>> indexers;
    indexers.push_back(targetGlobalIndexer_);

    panzer::BlockedTpetraLinearObjFactory<panzer::Traits,double,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType> factory(comm_,indexers);

    auto ownedMatrix = factory.getTpetraMatrix(0,0);
    auto ghostedMatrix = factory.getGhostedTpetraMatrix(0,0);
    ownedMatrix->resumeFill();
    ownedMatrix->setAllToScalar(0.0);
    ghostedMatrix->resumeFill();
    ghostedMatrix->setAllToScalar(0.0);

    const int fieldIndex = targetGlobalIndexer_->getFieldNum(targetBasisDescriptor_.getType());

    const bool is_scalar = targetBasisDescriptor_.getType()=="HGrad" || targetBasisDescriptor_.getType()=="Const" || targetBasisDescriptor_.getType()=="HVol";

    // Loop over element blocks and fill mass matrix
    if(is_scalar){
      auto M = ghostedMatrix->getLocalMatrixDevice();
      for (const auto& block : elementBlockNames_) {

        double ebMultiplier = 1.0;
        if (elementBlockMultipliers != nullptr)
          ebMultiplier = elementBlockMultipliers->find(block)->second;

        // Get the cell local ids and set the BasisValues object (BV
        // can come from a user defined set or from WorksetContainer).
        const panzer::BasisValues2<double>* bv_ptr;
        int num_cells_owned_ghosted_virtual = 0;
        int num_cells_owned = 0;
        Kokkos::View<const panzer::LocalOrdinal*> cell_local_ids;
        if (useUserSuppliedBasisValues_) {
          // Skip this block if there are no elements in this block on this mpi process
          auto tmp = connManager_->getElementBlock(block); // no ghosting or virtual in this case
          if (tmp.size() == 0)
            continue;

          Kokkos::View<panzer::LocalOrdinal*>::HostMirror cell_local_ids_host(tmp.data(),tmp.size());
          Kokkos::View<panzer::LocalOrdinal*> cell_local_ids_nonconst("cell_local_ids",tmp.size());
          Kokkos::deep_copy(cell_local_ids_nonconst,cell_local_ids_host);
          cell_local_ids = cell_local_ids_nonconst;

          bv_ptr = basisValues_[block].get();
          num_cells_owned_ghosted_virtual = cell_local_ids.extent(0);
          num_cells_owned = cell_local_ids.extent(0);
        }
        else {
          // Based on descriptor, currently assumes there should only
          // be one workset (partitioned path assumes a single
          // workset).
          panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,true);
          const auto worksets = worksetContainer_->getWorksets(wd);
          // Skip this block if there are no elements in this block on this mpi process
          if (worksets->size() == 0)
            continue;

          TEUCHOS_ASSERT(worksets->size() == 1);

          const auto& workset = (*worksets)[0];
          bv_ptr = &(workset.getBasisValues(targetBasisDescriptor_,integrationDescriptor_));
          num_cells_owned_ghosted_virtual = workset.numOwnedCells()+workset.numGhostCells()+workset.numVirtualCells();
          num_cells_owned = workset.numOwnedCells();
          cell_local_ids = workset.getLocalCellIDs();
        }
        const auto& basisValues = *bv_ptr;

        const auto unweightedBasis = basisValues.getBasisValues(false).get_static_view();
        const auto weightedBasis = basisValues.getBasisValues(true).get_static_view();

        const std::vector<panzer::LocalOrdinal>& offsets = targetGlobalIndexer_->getGIDFieldOffsets(block,fieldIndex);
        PHX::View<panzer::LocalOrdinal*> kOffsets("MassMatrix: Offsets",offsets.size());
        auto kOffsets_h = Kokkos::create_mirror_view(kOffsets);

        for(const auto& i : offsets)
          kOffsets_h(i) = offsets[i];

        Kokkos::deep_copy(kOffsets, kOffsets_h);

        // Local Ids
        PHX::View<panzer::LocalOrdinal**> localIds("MassMatrix: LocalIds", num_cells_owned_ghosted_virtual,
                                                   targetGlobalIndexer_->getElementBlockGIDCount(block));

        // Remove the ghosted cell ids or the call to getElementLocalIds will spill array bounds
        const auto cellLocalIdsNoGhost = Kokkos::subview(cell_local_ids,std::make_pair(0,num_cells_owned));

        targetGlobalIndexer_->getElementLIDs(cellLocalIdsNoGhost,localIds);

        const int numBasisPoints = static_cast<int>(weightedBasis.extent(1));
        if ( use_lumping ) {
          Kokkos::parallel_for(num_cells_owned,KOKKOS_LAMBDA (const int& cell) {
            double total_mass = 0.0, trace = 0.0;

            panzer::LocalOrdinal cLIDs[256];
            const int numIds = static_cast<int>(localIds.extent(1));
            for(int i=0;i<numIds;++i)
              cLIDs[i] = localIds(cell,i);

            double vals[256]={0.0};
            const int numQP = static_cast<int>(unweightedBasis.extent(2));

            for (int row=0; row < numBasisPoints; ++row) {
              for (int col=0; col < numIds; ++col) {
                for (int qp=0; qp < numQP; ++qp) {
                  auto tmp = unweightedBasis(cell,row,qp) * weightedBasis(cell,col,qp) * ebMultiplier;
                  total_mass += tmp;
                  if (col == row )
                    trace += tmp;
                }
              }
            }

            for (int row=0; row < numBasisPoints; ++row) {
              for (int col=0; col < numBasisPoints; ++col)
                vals[col] = 0;

              int offset = kOffsets(row);
              panzer::LocalOrdinal lid = localIds(cell,offset);
              int col = row;
              vals[col] = 0.0;
              for (int qp=0; qp < numQP; ++qp)
                vals[col] += unweightedBasis(cell,row,qp) * weightedBasis(cell,col,qp) * ebMultiplier * total_mass / trace;

              M.sumIntoValues(lid,cLIDs,numIds,vals,true,true);
            }
          });

        } else {
          Kokkos::parallel_for(num_cells_owned,KOKKOS_LAMBDA (const int& cell) {
            panzer::LocalOrdinal cLIDs[256];
            const int numIds = static_cast<int>(localIds.extent(1));
            for(int i=0;i<numIds;++i)
              cLIDs[i] = localIds(cell,i);

            double vals[256];
            const int numQP = static_cast<int>(unweightedBasis.extent(2));

            for (int row=0; row < numBasisPoints; ++row) {
              int offset = kOffsets(row);
              panzer::LocalOrdinal lid = localIds(cell,offset);

              for (int col=0; col < numIds; ++col) {
                vals[col] = 0.0;
                for (int qp=0; qp < numQP; ++qp)
                  vals[col] += unweightedBasis(cell,row,qp) * weightedBasis(cell,col,qp) * ebMultiplier;
              }
              M.sumIntoValues(lid,cLIDs,numIds,vals,true,true);

            }
          });
        }
      }
    } else {
      auto M = ghostedMatrix->getLocalMatrixDevice();
      for (const auto& block : elementBlockNames_) {

        double ebMultiplier = 1.0;
        if (elementBlockMultipliers != nullptr)
          ebMultiplier = elementBlockMultipliers->find(block)->second;

        // Get the cell local ids and set the BasisValues object (BV
        // can come from a user defined set or from WorksetContainer).
        const panzer::BasisValues2<double>* bv_ptr;
        int num_cells_owned_ghosted_virtual = 0;
        int num_cells_owned = 0;
        Kokkos::View<const panzer::LocalOrdinal*> cell_local_ids;
        if (useUserSuppliedBasisValues_) {
          // Skip this block if there are no elements in this block on this mpi process
          auto tmp = connManager_->getElementBlock(block); // no ghosting or virtual in this case
          if (tmp.size() == 0)
            continue;

          Kokkos::View<panzer::LocalOrdinal*>::HostMirror cell_local_ids_host(tmp.data(),tmp.size());
          Kokkos::View<panzer::LocalOrdinal*> cell_local_ids_nonconst("cell_local_ids",tmp.size());
          Kokkos::deep_copy(cell_local_ids_nonconst,cell_local_ids_host);
          cell_local_ids = cell_local_ids_nonconst;

          bv_ptr = basisValues_[block].get();
          num_cells_owned_ghosted_virtual = cell_local_ids.extent(0);
          num_cells_owned = cell_local_ids.extent(0);
        }
        else {
          // Based on descriptor, currently assumes there should only
          // be one workset (partitioned path assumes a single
          // workset).
          panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,true);
          const auto worksets = worksetContainer_->getWorksets(wd);
          // Skip this block if there are no elements in this block on this mpi process
          if (worksets->size() == 0)
            continue;

          TEUCHOS_ASSERT(worksets->size() == 1);

          const auto& workset = (*worksets)[0];
          bv_ptr = &(workset.getBasisValues(targetBasisDescriptor_,integrationDescriptor_));
          num_cells_owned_ghosted_virtual = workset.numOwnedCells()+workset.numGhostCells()+workset.numVirtualCells();
          num_cells_owned = workset.numOwnedCells();
          cell_local_ids = workset.getLocalCellIDs();
        }
        const auto& basisValues = *bv_ptr;

        const auto unweightedBasis = basisValues.getVectorBasisValues(false).get_static_view();
        const auto weightedBasis = basisValues.getVectorBasisValues(true).get_static_view();

        const std::vector<panzer::LocalOrdinal>& offsets = targetGlobalIndexer_->getGIDFieldOffsets(block,fieldIndex);
        PHX::View<panzer::LocalOrdinal*> kOffsets("MassMatrix: Offsets",offsets.size());
        auto kOffsets_h = Kokkos::create_mirror_view(kOffsets);

        for(const auto& i : offsets)
          kOffsets_h(i) = offsets[i];

        Kokkos::deep_copy(kOffsets, kOffsets_h);

        // Local Ids
        PHX::View<panzer::LocalOrdinal**> localIds("MassMatrix: LocalIds", num_cells_owned_ghosted_virtual,
                                                   targetGlobalIndexer_->getElementBlockGIDCount(block));

        // Remove the ghosted cell ids or the call to getElementLocalIds will spill array bounds
        const PHX::View<const int*> cellLocalIdsNoGhost = Kokkos::subview(cell_local_ids,std::make_pair(0,num_cells_owned));

        targetGlobalIndexer_->getElementLIDs(cellLocalIdsNoGhost,localIds);

        const int numBasisPoints = static_cast<int>(weightedBasis.extent(1));
        Kokkos::parallel_for(num_cells_owned,KOKKOS_LAMBDA (const int& cell) {

          panzer::LocalOrdinal cLIDs[256];
          const int numIds = static_cast<int>(localIds.extent(1));
          for(int i=0;i<numIds;++i)
            cLIDs[i] = localIds(cell,i);

          double vals[256];
          const int numQP = static_cast<int>(unweightedBasis.extent(2));

          for (int qp=0; qp < numQP; ++qp) {
            for (int row=0; row < numBasisPoints; ++row) {
              int offset = kOffsets(row);
              panzer::LocalOrdinal lid = localIds(cell,offset);

              for (int col=0; col < numIds; ++col){
                vals[col] = 0.0;
                for(int dim=0; dim < static_cast<int>(weightedBasis.extent(3)); ++dim)
                  vals[col] += unweightedBasis(cell,row,qp,dim) * weightedBasis(cell,col,qp,dim) * ebMultiplier;
              }

              M.sumIntoValues(lid,cLIDs,numIds,vals,true,true);
            }
          }
        });
      }
    }

    {
      PANZER_FUNC_TIME_MONITOR_DIFF("Exporting of mass matrix",ExportMM);
      auto map = factory.getMap(0);
      ghostedMatrix->fillComplete(map,map);
      const auto exporter = factory.getGhostedExport(0);
      ownedMatrix->doExport(*ghostedMatrix, *exporter, Tpetra::ADD);
      ownedMatrix->fillComplete();
    }
    return ownedMatrix;
  }

  Teuchos::RCP<Tpetra::MultiVector<double,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>>
  panzer::L2Projection::buildInverseLumpedMassMatrix()
  {
    PANZER_FUNC_TIME_MONITOR_DIFF("L2Projection<panzer::LocalOrdinal,panzer::GlobalOrdinal>::buildInverseLumpedMassMatrix",buildInvLMM);
    using Teuchos::rcp;
    const auto massMatrix = this->buildMassMatrix(true);
    const auto lumpedMassMatrix = rcp(new Tpetra::MultiVector<double,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>(massMatrix->getDomainMap(),1,true));
    const auto tmp = rcp(new Tpetra::MultiVector<double,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>(massMatrix->getRangeMap(),1,false));
    tmp->putScalar(1.0);
    {
      PANZER_FUNC_TIME_MONITOR_DIFF("Apply",Apply);
      massMatrix->apply(*tmp,*lumpedMassMatrix);
    }
    {
      PANZER_FUNC_TIME_MONITOR_DIFF("Reciprocal",Reciprocal);
      lumpedMassMatrix->reciprocal(*lumpedMassMatrix);
    }
    return lumpedMassMatrix;
  }

  Teuchos::RCP<Tpetra::CrsMatrix<double,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>>
  panzer::L2Projection::buildRHSMatrix(const panzer::GlobalIndexer& sourceGlobalIndexer,
                                       const Teuchos::RCP<const Tpetra::Map<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>>& inputOwnedSourceMap,
                                       const std::string& sourceFieldName,
                                       const panzer::BasisDescriptor& sourceBasisDescriptor,
                                       const int directionIndex)
  {
    // *******************
    // Create Retangular matrix (both ghosted and owned).
    // *******************
    using Teuchos::RCP;
    using Teuchos::rcp;
    using MapType = Tpetra::Map<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>;
    using GraphType = Tpetra::CrsGraph<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>;
    using ExportType = Tpetra::Export<panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>;
    using MatrixType = Tpetra::CrsMatrix<double,panzer::LocalOrdinal,panzer::GlobalOrdinal,panzer::TpetraNodeType>;

    // *******************
    // Build ghosted graph
    // *******************

    RCP<MapType> ghostedTargetMap;
    {
      std::vector<panzer::GlobalOrdinal> indices;
      targetGlobalIndexer_->getOwnedAndGhostedIndices(indices);
      ghostedTargetMap = rcp(new MapType(Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::invalid(),indices,0,comm_));
    }

    RCP<MapType> ghostedSourceMap;
    {
      std::vector<panzer::GlobalOrdinal> indices;
      sourceGlobalIndexer.getOwnedAndGhostedIndices(indices);
      ghostedSourceMap = rcp(new MapType(Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::invalid(),indices,0,comm_));
    }


    // Now insert the non-zero pattern per row
    // count number of entries per row; required by CrsGraph constructor
    std::vector<size_t> nEntriesPerRow(ghostedTargetMap->getLocalNumElements(),0);
    std::vector<std::string> elementBlockIds;
    targetGlobalIndexer_->getElementBlockIds(elementBlockIds);
    std::vector<std::string>::const_iterator blockItr;
    for (blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;
      const std::vector<panzer::LocalOrdinal> & elements = targetGlobalIndexer_->getElementBlock(blockId);

      std::vector<panzer::GlobalOrdinal> row_gids;
      std::vector<panzer::GlobalOrdinal> col_gids;

      for(std::size_t elmt=0;elmt<elements.size();elmt++) {
        targetGlobalIndexer_->getElementGIDs(elements[elmt],row_gids);
        sourceGlobalIndexer.getElementGIDs(elements[elmt],col_gids);
        for(std::size_t row=0;row<row_gids.size();row++) {
          panzer::LocalOrdinal lid =
                               ghostedTargetMap->getLocalElement(row_gids[row]);
          nEntriesPerRow[lid] += col_gids.size();
        }
      }
    }

    Teuchos::ArrayView<const size_t> nEntriesPerRowView(nEntriesPerRow);
    RCP<GraphType> ghostedGraph = rcp(new GraphType(ghostedTargetMap,ghostedSourceMap,nEntriesPerRowView));

    for (blockItr=elementBlockIds.begin();blockItr!=elementBlockIds.end();++blockItr) {
      std::string blockId = *blockItr;
      const std::vector<panzer::LocalOrdinal> & elements = targetGlobalIndexer_->getElementBlock(blockId);

      std::vector<panzer::GlobalOrdinal> row_gids;
      std::vector<panzer::GlobalOrdinal> col_gids;

      for(std::size_t elmt=0;elmt<elements.size();elmt++) {
        targetGlobalIndexer_->getElementGIDs(elements[elmt],row_gids);
        sourceGlobalIndexer.getElementGIDs(elements[elmt],col_gids);
        for(std::size_t row=0;row<row_gids.size();row++)
          ghostedGraph->insertGlobalIndices(row_gids[row],col_gids);
      }
    }

    RCP<MapType> ownedTargetMap;
    {
      std::vector<panzer::GlobalOrdinal> indices;
      targetGlobalIndexer_->getOwnedIndices(indices);
      ownedTargetMap = rcp(new MapType(Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::invalid(),indices,0,comm_));
    }

    RCP<const MapType> ownedSourceMap = inputOwnedSourceMap;
    if (is_null(ownedSourceMap)) {
      std::vector<panzer::GlobalOrdinal> indices;
      sourceGlobalIndexer.getOwnedIndices(indices);
      ownedSourceMap = rcp(new MapType(Teuchos::OrdinalTraits<panzer::GlobalOrdinal>::invalid(),indices,0,comm_));
    }

    // Fill complete with owned range and domain map
    ghostedGraph->fillComplete(ownedSourceMap,ownedTargetMap);

    // *****************
    // Build owned graph
    // *****************

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

    // *******************
    // Fill ghosted matrix
    // *******************
    for (const auto& block : elementBlockNames_) {

      panzer::WorksetDescriptor wd(block,panzer::WorksetSizeType::ALL_ELEMENTS,true,true);
      const auto& worksets = worksetContainer_->getWorksets(wd);
      for (const auto& workset : *worksets) {

        // Get target basis values: current implementation assumes target basis is HGrad
        const auto& targetBasisValues = workset.getBasisValues(targetBasisDescriptor_,integrationDescriptor_);
        const auto& targetWeightedBasis = targetBasisValues.getBasisValues(true).get_static_view();

        // Sources can be any basis
        const auto& sourceBasisValues = workset.getBasisValues(sourceBasisDescriptor,integrationDescriptor_);
        PHX::View<const double***> sourceUnweightedScalarBasis;
        PHX::View<const double****> sourceUnweightedVectorBasis;
        bool useRankThreeBasis = false; // default to gradient or vector basis
        if ( (sourceBasisDescriptor.getType() == "HGrad") || (sourceBasisDescriptor.getType() == "Const") || (sourceBasisDescriptor.getType() == "HVol") ) {
          if (directionIndex == -1) { // Project dof value
            sourceUnweightedScalarBasis = sourceBasisValues.getBasisValues(false).get_static_view();
            useRankThreeBasis = true;
          }
          else { // Project dof gradient of scalar basis
            sourceUnweightedVectorBasis = sourceBasisValues.getGradBasisValues(false).get_static_view();
          }
        }
        else { // Project vector value
          sourceUnweightedVectorBasis = sourceBasisValues.getVectorBasisValues(false).get_static_view();
        }

        // Get the element local ids
        PHX::View<panzer::LocalOrdinal**> targetLocalIds("buildRHSMatrix: targetLocalIds", workset.numOwnedCells(),
                                                      targetGlobalIndexer_->getElementBlockGIDCount(block));
        PHX::View<panzer::LocalOrdinal**> sourceLocalIds("buildRHSMatrix: sourceLocalIds", workset.numOwnedCells(),
                                                      sourceGlobalIndexer.getElementBlockGIDCount(block));
        {
          // Remove the ghosted cell ids or the call to getElementLocalIds will spill array bounds
          const auto cellLocalIdsNoGhost = Kokkos::subview(workset.cell_local_ids_k,std::make_pair(0,workset.numOwnedCells()));
          targetGlobalIndexer_->getElementLIDs(cellLocalIdsNoGhost,targetLocalIds);
          sourceGlobalIndexer.getElementLIDs(cellLocalIdsNoGhost,sourceLocalIds);
        }

        // Get the offsets
        PHX::View<panzer::LocalOrdinal*> targetFieldOffsets;
        {
          const auto fieldIndex = targetGlobalIndexer_->getFieldNum(targetBasisDescriptor_.getType());
          const std::vector<panzer::LocalOrdinal>& offsets = targetGlobalIndexer_->getGIDFieldOffsets(block,fieldIndex);
          targetFieldOffsets = PHX::View<panzer::LocalOrdinal*>("L2Projection:buildRHS:targetFieldOffsets",offsets.size());
          const auto hostOffsets = Kokkos::create_mirror_view(targetFieldOffsets);
          for(size_t i=0; i < offsets.size(); ++i)
            hostOffsets(i) = offsets[i];
          Kokkos::deep_copy(targetFieldOffsets,hostOffsets);
        }

        PHX::View<panzer::LocalOrdinal*> sourceFieldOffsets;
        {
          const auto fieldIndex = sourceGlobalIndexer.getFieldNum(sourceFieldName);
          const std::vector<panzer::LocalOrdinal>& offsets = sourceGlobalIndexer.getGIDFieldOffsets(block,fieldIndex);
          TEUCHOS_TEST_FOR_EXCEPTION(offsets.size() == 0, std::runtime_error,
                                     "ERROR: panzer::L2Projection::buildRHSMatrix() - The source field, \""
                                     << sourceFieldName << "\", does not exist in element block \""
                                     << block << "\"!");
          sourceFieldOffsets = PHX::View<panzer::LocalOrdinal*>("L2Projection:buildRHS:sourceFieldOffsets",offsets.size());
          const auto hostOffsets = Kokkos::create_mirror_view(sourceFieldOffsets);
          for(size_t i=0; i <offsets.size(); ++i)
            hostOffsets(i) = offsets[i];
          Kokkos::deep_copy(sourceFieldOffsets,hostOffsets);
        }

        const auto localMatrix = ghostedMatrix->getLocalMatrixDevice();
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
          panzer::LocalOrdinal cLIDs[256];
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
    ghostedMatrix->fillComplete(ownedSourceMap,ownedTargetMap);
    ownedMatrix->resumeFill();
    ownedMatrix->doExport(*ghostedMatrix,*exporter,Tpetra::ADD);
    ownedMatrix->fillComplete(ownedSourceMap,ownedTargetMap);

    return ownedMatrix;
  }

}

#endif
