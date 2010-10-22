
#ifndef PANZER_CELL_DATA_HPP
#define PANZER_CELL_DATA_HPP

#include "Teuchos_TestForException.hpp"

namespace panzer {
  
  class CellData {
    
  public:
    
    explicit CellData(std::size_t num_cells, int base_cell_dimension) :
      m_num_cells(num_cells),
      m_dimension(base_cell_dimension),
      m_is_side(false),
      m_side(-1) 
    { }
    
    explicit CellData(std::size_t num_cells, int base_cell_dimension,
		      int local_side_id) : 
      m_num_cells(num_cells),
      m_dimension(base_cell_dimension),
      m_is_side(true),
      m_side(local_side_id) 
    { }
	
    bool isSide() const 
    {return m_is_side;}
    
    int side() const 
    {
      TEST_FOR_EXCEPTION(!m_is_side, std::logic_error,
			 "Cannot return side index, CellData is not a side!");
      return m_side;
    }
    	
    std::size_t numCells() const 
    {return m_num_cells;}
    	
    //! Dimension of the base cell.  NOT the dimension of the local side, even if the side() method returns true.
    int baseCellDimension() const 
    {return m_dimension;}
    
  private:
    std::size_t m_num_cells;
    int m_dimension;
    bool m_is_side;
    int m_side;
      
  };

}

#endif
