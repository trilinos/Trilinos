#include <stk_mesh/base/FieldRestriction.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <sstream>

namespace stk {
namespace mesh {

void FieldRestriction::print(
  std::ostream & os,
  const EntityRank & entity_rank,
  const Selector & selector,
  FieldArrayRank field_rank
  ) const
{
  os << "FieldRestriction[ entity_rank: " << entity_rank << " selector: \"" << selector << "\", dimensions: " ;
  os << m_stride[0] ;
  for ( FieldArrayRank i = 1 ; i < field_rank ; ++i ) {
    if ( ! m_stride[i] ) {
      os << " , 0 " ;
    }
    else if ( m_stride[i] % m_stride[i-1] ) {
      os << " , " << m_stride[i] << " / " << m_stride[i-1] ;
    }
    else {
      os << " , " << m_stride[i] / m_stride[i-1] ;
    }
  }
  os << " ]" ;
}

std::string print_restriction(
  const FieldRestriction & restr,
  const EntityRank & entity_rank,
  const Selector& selector,
  FieldArrayRank field_rank
                              )
{
  std::ostringstream oss;
  restr.print(oss, entity_rank, selector, field_rank);
  return oss.str();
}

} // namespace mesh
} // namespace stk
