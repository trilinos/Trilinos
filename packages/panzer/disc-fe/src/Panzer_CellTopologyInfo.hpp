// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CELLTOPOLOGYINFO_HPP
#define PANZER_CELLTOPOLOGYINFO_HPP

#include <string>
#include "Teuchos_RCP.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Kokkos_DynRankView.hpp"
#include "Intrepid2_Basis.hpp"

#include "Panzer_Dimension.hpp"
#include "Panzer_CellData.hpp"
#include "Panzer_IntrepidBasisFactory.hpp"


// Added by Suzey: 06/18/2012, to obtain the edge information of a cell topology
// and to create the edge data layouts

namespace panzer {

  class CellTopologyInfo { 

  public:
    
    CellTopologyInfo(int numCells, const Teuchos::RCP<const shards::CellTopology>& cellTopo);

    int getNumCells() const
    { return num_cells; }
    
    int getDimension() const
    { return dimension; }
    
    int getNumEdges() const
    { return num_edges; }
    
    std::string getCellName() const
    { return cell_topo_name; }
    
    Teuchos::RCP<const shards::CellTopology> getCellTopology() const
    { return topology; }

  public:

    //!<Cell,Edge>
    Teuchos::RCP<PHX::DataLayout> edge_scalar;  

    //!<Cell,Edge,Dim>
    Teuchos::RCP<PHX::DataLayout> edge_vector; 

  private:

    //! Initialize data layouts
    void initializeDataLayouts();

    Teuchos::RCP<const shards::CellTopology> topology;

    std::string cell_topo_name;

    int num_cells;
    int dimension;
    int num_edges; 

  };

}

#endif
