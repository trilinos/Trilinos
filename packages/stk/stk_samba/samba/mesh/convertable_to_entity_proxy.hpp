#ifndef SAMBA_SAMBA_MESH_CONVERTABLE_TO_ENTITY_PROXY_HPP
#define SAMBA_SAMBA_MESH_CONVERTABLE_TO_ENTITY_PROXY_HPP

#include <samba/partition_offset.hpp>

namespace samba {

class entity_proxy;

namespace detail {

class partition_impl;

//used to break cyclic dependency
class convertable_to_entity_proxy
{
  public:
    convertable_to_entity_proxy( partition_impl const * arg_partition
                                   ,partition_offset arg_offset
                                  )
      : m_partition(arg_partition)
      , m_offset(arg_offset)
    {}

  private:
    partition_impl const* m_partition;
    partition_offset         m_offset;

    friend class samba::entity_proxy;
};

} } //namespace samba::detail

#endif //SAMBA_SAMBA_MESH_CONVERTABLE_TO_ENTITY_PROXY_HPP


