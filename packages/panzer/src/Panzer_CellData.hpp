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

#ifndef PANZER_CELL_DATA_HPP
#define PANZER_CELL_DATA_HPP

#include "Panzer_config.hpp"

#include "Teuchos_Assert.hpp"
#include "Teuchos_RCP.hpp"

#include "Shards_CellTopology.hpp"
#include "Shards_BasicTopologies.hpp"

namespace panzer {
  
  /** \brief Data for determining cell topology and dimensionality

      This class provides information about the finite elements/cells
      and allows us to reuse all the evaluation machinery for both
      volume and surface integrations.  When performing neumann
      conditions, we integrate over a side but must still use all
      basis points from the volume if we need a gradient that contains
      a normal component to the surface being integrated over.
  */ 
  class CellData {
    
  public:

    /** for testing purposes only */
    explicit CellData() 
    { std::cout << "WARNING: Default constructor for panzer::CellData is for testing purposes only!" << std::endl; } 
    
    /** Build cell data that uses volume data.  CellTopology is on the
      * volume cells!
      */
    explicit CellData(std::size_t num_cells, int base_cell_dimension, 
                      const Teuchos::RCP<const shards::CellTopology> & ct) :
      m_num_cells(num_cells),
      m_dimension(base_cell_dimension),
      m_is_side(false),
      m_side(-1),
      m_cell_topo(ct)
    { }

    /** Build cell data that uses side cells.  CellTopology is on the
      * volume cells!
      */
    explicit CellData(std::size_t num_cells, int base_cell_dimension,
		      int local_side_id,const Teuchos::RCP<const shards::CellTopology> & ct) :
      m_num_cells(num_cells),
      m_dimension(base_cell_dimension),
      m_is_side(true),
      m_side(local_side_id),
      m_cell_topo(ct)
    { }
	
    bool isSide() const 
    {return m_is_side;}
    
    int side() const 
    {
      TEUCHOS_TEST_FOR_EXCEPTION(!m_is_side, std::logic_error,
			 "Cannot return side index, CellData is not a side!");
      return m_side;
    }
    	
    std::size_t numCells() const 
    {return m_num_cells;}
    	
    //! Dimension of the base cell.  NOT the dimension of the local side, even if the side() method returns true.
    int baseCellDimension() const 
    {return m_dimension;}
    
    //! Get CellTopology for the base cell.
    Teuchos::RCP<const shards::CellTopology> getCellTopology() const
    { return m_cell_topo; }

  private:
    std::size_t m_num_cells;
    int m_dimension;
    bool m_is_side;
    int m_side;
      
 
    Teuchos::RCP<const shards::CellTopology> m_cell_topo;
  };

}

#endif
