// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_CellTopologyInfo.hpp"
#include "Panzer_IntegrationRule.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"

#include "Teuchos_Assert.hpp"
#include "Phalanx_DataLayout_MDALayout.hpp"

panzer::CellTopologyInfo::
CellTopologyInfo(int numCells, const Teuchos::RCP<const shards::CellTopology>& cellTopo)
{
  num_cells = numCells; 

  dimension = cellTopo->getDimension();
  num_edges = cellTopo->getEdgeCount();
  cell_topo_name = cellTopo->getName();
  
  topology = cellTopo;

  initializeDataLayouts();
}


void panzer::CellTopologyInfo::initializeDataLayouts()
{
  using Teuchos::rcp;
  using PHX::MDALayout;

  edge_scalar = rcp(new MDALayout<Cell,Edge>(num_cells, num_edges));
  edge_vector = rcp(new MDALayout<Cell,Edge,Dim>(num_cells, num_edges, dimension));

}
