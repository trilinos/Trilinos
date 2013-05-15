#include <stk_mesh/base/FieldRestriction.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <sstream>

namespace stk {
namespace mesh {

void FieldRestriction::print(
  std::ostream & os,
  const EntityRank & entity_rank,
  const Part & part,
  FieldArrayRank field_rank
  ) const
{
  MetaData const& meta = MetaData::get(part);
  os << "FieldRestriction[ entity_rank: " << meta.entity_rank_name(entity_rank) << " part: \"" << part.name() << "\", dimensions: " ;
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
  const Part & part,
  FieldArrayRank field_rank
                              )
{
  std::ostringstream oss;
  restr.print(oss, entity_rank, part, field_rank);
  return oss.str();
}

} // namespace mesh
} // namespace stk
