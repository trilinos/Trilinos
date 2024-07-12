// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef PANZER_CELL_DATA_HPP
#define PANZER_CELL_DATA_HPP

#include "PanzerDiscFE_config.hpp"

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
    { }
    
    /** Build cell data that uses volume data.  CellTopology is on the
      * volume cells!
      */
    explicit CellData(std::size_t num_cells,
                      const Teuchos::RCP<const shards::CellTopology> & ct) :
      m_num_cells(num_cells),
      m_is_side(false),
      m_side(-1),
      m_cell_topo(ct)
    { }

    /** Build cell data that uses side cells.  CellTopology is on the
      * volume cells!
      */
    explicit CellData(std::size_t num_cells,
		      int local_side_id,const Teuchos::RCP<const shards::CellTopology> & ct) :
      m_num_cells(num_cells),
      m_is_side(local_side_id >= 0),
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
    {return m_cell_topo->getDimension();}
    
    //! Get CellTopology for the base cell.
    // TODO BWR this terminology is a bit confusing... this is not the base cell topology
    Teuchos::RCP<const shards::CellTopology> getCellTopology() const
    { return m_cell_topo; }

  private:
    std::size_t m_num_cells;
    bool m_is_side;
    int m_side;
      
 
    Teuchos::RCP<const shards::CellTopology> m_cell_topo;
  };

}

#endif
