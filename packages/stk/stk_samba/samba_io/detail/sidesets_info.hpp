#ifndef SAMBA_SAMBA_IO_DETAIL_SIDESETS_INFO_HPP
#define SAMBA_SAMBA_IO_DETAIL_SIDESETS_INFO_HPP

#include <samba/mesh.hpp>
#include <samba/field.hpp>

namespace samba {
namespace io {
namespace detail {

struct SideBlock
{
  std::string         m_block_name;
  entity_key_interval m_sides;
  // entity_topology  m_parent_topology;
  // entity_block_key m_entity_block;
};


struct SidesetImpl
{
  typedef boost::shared_ptr<SidesetImpl> Ptr;

  std::string                  m_name;
  entity_block_key             m_entity_block;
  entity_rank                  m_rank;
  std::vector<SideBlock>       m_blocks;
};

} // namespace detail


typedef samba::io::detail::SidesetImpl Sideset;
typedef samba::io::detail::SideBlock   SideBlock;
typedef Sideset::Ptr                   SidesetPtr;


} // namespace io
} // namespace samba

#endif
