//#define HAVE_SHARDS_DEBUG

#include <stdexcept>
#include <sstream>

#include <stk_mesh/fem/CellTopology.hpp>

#include <stk_util/environment/ReportHandler.hpp>

namespace stk {
namespace mesh {
namespace fem {

typedef CellTopologyData_Subcell Subcell ;

void CellTopology::requireCell() const
{
  ThrowErrorMsgIf( m_cell       == NULL, "m_cell is NULL" );
  ThrowErrorMsgIf( m_cell->base == NULL, m_cell->name << " has NULL base" );
}

void CellTopology::requireDimension( const unsigned subcellDim ) const
{
  ThrowInvalidArgMsgIf( subcellDim > 3,
                        "dim = " << subcellDim << " > 3 )" );
}

void CellTopology::requireSubcell( const unsigned subcellDim ,
                                   const unsigned subcellOrd ) const
{
  ThrowInvalidArgMsgIf( m_cell->subcell_count[ subcellDim ] <= subcellOrd,
                        "ord = " << subcellOrd << " >= '" << m_cell->name <<
                        "'.subcell_count[" << subcellDim << "] = " <<
                        m_cell->subcell_count[ subcellDim ]);
}

void CellTopology::requireNodeMap( const unsigned subcellDim ,
                                   const unsigned subcellOrd ,
                                   const unsigned nodeOrd ) const
{
  const unsigned n =
    m_cell->subcell[subcellDim][subcellOrd].topology->node_count ;
  ThrowInvalidArgMsgIf( n <= nodeOrd,
                        nodeOrd << " >= '" << m_cell->name <<
                        "'.subcell[" << subcellDim << "][" << subcellOrd <<
                        "].topology->node_count = " << n);
}

void CellTopology::requireNodePermutation( const unsigned permutationOrd ,
                                           const unsigned nodeOrd ) const
{
  ThrowInvalidArgMsgIf( m_cell->permutation_count <= permutationOrd,
                        permutationOrd << " >= " << m_cell->permutation_count );
  ThrowInvalidArgMsgIf( m_cell->node_count <= nodeOrd,
                        nodeOrd << " >= " << m_cell->node_count );
}

std::ostream & operator << ( std::ostream & os, const CellTopology & cell)
{
  const CellTopologyData* cell_top_data = cell.getCellTopologyData();
  return shards::operator<<(os, *cell_top_data);
}

} // namespace fem
} // namespace mesh
} // namespace stk
