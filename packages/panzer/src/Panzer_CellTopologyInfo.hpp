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

#ifndef PANZER_CELLTOPOLOGYINFO_HPP
#define PANZER_CELLTOPOLOGYINFO_HPP

#include <string>
#include "Teuchos_RCP.hpp"

#include "Phalanx_DataLayout.hpp"

#include "Intrepid_FieldContainer.hpp"
#include "Intrepid_Basis.hpp"

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
