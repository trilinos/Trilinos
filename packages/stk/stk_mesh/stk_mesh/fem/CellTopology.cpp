//#define HAVE_SHARDS_DEBUG

#include <stdexcept>
#include <sstream>
#include <stk_mesh/fem/CellTopology.hpp>

namespace stk {
namespace mesh {
namespace fem {

typedef CellTopologyData_Subcell Subcell ;

void CellTopology::requireCell() const
{
  if ( m_cell == NULL || m_cell->base == NULL ) {
    std::ostringstream oss ;
    oss << "shards::CellTopology::requireCell() : FAILED " ;
    if ( m_cell == NULL ) {
      oss << "is NULL" ;
    }
    else {
      oss << "'" << m_cell->name << "' has NULL base";
    }
    oss << " ) FAILED";
    throw std::runtime_error( oss.str() );
  }
}

void CellTopology::requireDimension( const unsigned subcellDim ) const
{
  if ( 3 < subcellDim ) {
    std::ostringstream oss ;
    oss << "shards::CellTopology::requireDimension( ERROR: dim = "
        << subcellDim << " > 3 )" ;
    throw std::invalid_argument( oss.str() );
  }
}

void CellTopology::requireSubcell( const unsigned subcellDim ,
                                   const unsigned subcellOrd ) const
{
  if ( m_cell->subcell_count[ subcellDim ] <= subcellOrd ) {
    std::ostringstream oss ;
    oss << "shards::CellTopology::requireSubcell( dim = "
        << subcellDim << " , ERROR: ord = " << subcellOrd
        << " > '" << m_cell->name
        << "'.subcell_count[" << subcellDim
        << "] = " << m_cell->subcell_count[ subcellDim ]
        << " )" ;
    throw std::invalid_argument( oss.str() );
  }
}

void CellTopology::requireNodeMap( const unsigned subcellDim ,
                                   const unsigned subcellOrd ,
                                   const unsigned nodeOrd ) const
{
  const unsigned n =
    m_cell->subcell[subcellDim][subcellOrd].topology->node_count ;

  if ( n <= nodeOrd ) {
    std::ostringstream oss ;
    oss << "shards::CellTopology::requireNodeMap( " 
        << subcellDim << " , "
        << subcellOrd
        << " , ERROR: " << nodeOrd << " >= '"
        << m_cell->name 
        << "'.subcell[" << subcellDim
        << "][" << subcellOrd
        << "].topology->node_count = "
        << n << " )" ;
    throw std::invalid_argument( oss.str() );
  }
}

void CellTopology::requireNodePermutation( const unsigned permutationOrd ,
                                           const unsigned nodeOrd ) const
{
  const bool bad_p = m_cell->permutation_count <= permutationOrd ;
  const bool bad_n = m_cell->node_count        <= nodeOrd ;
  if ( bad_p || bad_n ) {
    std::ostringstream oss ;
    oss << "shards::CellTopology::requireNodePermutation( " ;
    if ( bad_p ) {
      oss << " ERROR: " << permutationOrd << " >= "
          << m_cell->permutation_count ;
    }
    else {
      oss << permutationOrd ;
    }
    oss << " , " ;
    if ( bad_n ) {
      oss << " ERROR: " << nodeOrd << " >= " << m_cell->node_count ;
    }
    else {
      oss << nodeOrd ;
    }
    oss << " )" ;
    throw std::invalid_argument( oss.str() );
  }
}

std::ostream & operator << ( std::ostream & os, const CellTopology & cell)
{
  const CellTopologyData* cell_top_data = cell.getTopologyData();
  return shards::operator<<(os, *cell_top_data);
}

} // namespace fem
} // namespace mesh
} // namespace stk


