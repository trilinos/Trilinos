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

#include "PanzerDiscFE_config.hpp"
#include "Panzer_GlobalIndexer.hpp"
#include "Panzer_IntrepidOrientation.hpp"

namespace panzer {

  void buildIntrepidOrientation(std::vector<Intrepid2::Orientation> & orientation,
                                panzer::ConnManager & connMgr)
  {
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::RCP;
    using Teuchos::rcp;

    orientation.clear();

    // Retrive element blocks and its meta data
    const int numElementBlocks = connMgr.numElementBlocks();

    std::vector<std::string> elementBlockIds;
    std::vector<shards::CellTopology> elementBlockTopologies;

    connMgr.getElementBlockIds(elementBlockIds);
    connMgr.getElementBlockTopologies(elementBlockTopologies);

    TEUCHOS_TEST_FOR_EXCEPTION(numElementBlocks <= 0 &&
                               numElementBlocks != static_cast<int>(elementBlockIds.size()) &&
                               numElementBlocks != static_cast<int>(elementBlockTopologies.size()),
                               std::logic_error,
                               "panzer::buildIntrepidOrientation: Number of element blocks does not match to element block meta data");

    // Currently panzer support only one type of elements for whole mesh (use the first cell topology)
    const auto cellTopo = elementBlockTopologies.at(0);
    const int numVerticesPerCell = cellTopo.getVertexCount();

    const auto fp = NodalFieldPattern(cellTopo);
    connMgr.buildConnectivity(fp);

    // Count and pre-alloc orientations
    int total_elems = 0;
    for (int i=0;i<numElementBlocks;++i) {
      total_elems += connMgr.getElementBlock(elementBlockIds.at(i)).size();
    }

    orientation.resize(total_elems);
    // Loop over element blocks
    for (int i=0;i<numElementBlocks;++i) {
      // get elements in a block
      const auto &elementBlock = connMgr.getElementBlock(elementBlockIds.at(i));

      const int numElementsPerBlock = elementBlock.size();

      // construct orientation information
      for (int c=0;c<numElementsPerBlock;++c) {
        const int localCellId = elementBlock.at(c);
        Kokkos::View<const panzer::GlobalOrdinal*,Kokkos::HostSpace>
          vertices(connMgr.getConnectivity(localCellId), numVerticesPerCell);
          // This function call expects a view for the vertices, not the nodes
        orientation[localCellId] = (Intrepid2::Orientation::getOrientation(cellTopo, vertices));
      }
    }
  }

  Teuchos::RCP<std::vector<Intrepid2::Orientation> >
  buildIntrepidOrientation(const Teuchos::RCP<const GlobalIndexer> globalIndexer)
  {
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::RCP;
    using Teuchos::rcp;

    auto orientation = rcp(new std::vector<Intrepid2::Orientation>);

    {
      RCP<const GlobalIndexer> ugi
        = rcp_dynamic_cast<const GlobalIndexer>(globalIndexer);

      if (ugi!=Teuchos::null) {
        const auto connMgr = ugi->getConnManager()->noConnectivityClone();

        TEUCHOS_TEST_FOR_EXCEPTION(connMgr == Teuchos::null,std::logic_error,
                                   "panzer::buildIntrepidOrientation: ConnManager is null!");

        buildIntrepidOrientation(*orientation, *connMgr);
        return orientation;
      }
    }

    TEUCHOS_TEST_FOR_EXCEPTION(true,std::logic_error,
                               "panzer::buildIntrepidOrientation: Could not cast GlobalIndexer");
  }

  void buildIntrepidOrientations(const std::vector<std::string>& eBlockNames,
                                 const panzer::ConnManager & connMgrInput,
                                 std::map<std::string,std::vector<Intrepid2::Orientation>> & orientations)
  {
    using Teuchos::rcp_dynamic_cast;
    using Teuchos::RCP;
    using Teuchos::rcp;

    auto connMgrPtr = connMgrInput.noConnectivityClone();
    auto& connMgr = *connMgrPtr;

    // Map element block name to topology
    std::map<std::string,shards::CellTopology> eb_name_to_topo;
    {
      std::vector<std::string> elementBlockIds;
      connMgr.getElementBlockIds(elementBlockIds);

      std::vector<shards::CellTopology> elementBlockTopologies;
      connMgr.getElementBlockTopologies(elementBlockTopologies);

      for (size_t i=0; i < elementBlockIds.size(); ++i) {
        eb_name_to_topo[elementBlockIds[i]] = elementBlockTopologies[i];
      }
    }

    // Currently panzer supports only one element topology for whole mesh (use the first cell topology)
    const auto cellTopo = eb_name_to_topo[eBlockNames[0]];
    const int numVerticesPerCell = cellTopo.getVertexCount();

    const auto fp = NodalFieldPattern(cellTopo);
    connMgr.buildConnectivity(fp);

    // Loop over requested element blocks
    for (size_t i=0;i<eBlockNames.size();++i) {

      const auto &lids = connMgr.getElementBlock(eBlockNames[i]);
      const size_t num_elems = lids.size();
      auto& orts = orientations[eBlockNames[i]];
      orts.resize(num_elems);

      for (size_t c=0;c<num_elems;++c) {
        const int lid = lids[c];
        Kokkos::View<const panzer::GlobalOrdinal*,Kokkos::HostSpace> 
          vertices(connMgr.getConnectivity(lid),numVerticesPerCell);

        // This function call expects a view for the vertices, not the nodes
        orts[c] = Intrepid2::Orientation::getOrientation(cellTopo, vertices);
      }
    }
  }

} // end namespace panzer
