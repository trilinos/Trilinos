#ifndef SAMBA_SAMBA_MESH_CONVERTABLE_TO_PARTITION_PROXY_HPP
#define SAMBA_SAMBA_MESH_CONVERTABLE_TO_PARTITION_PROXY_HPP

namespace samba {

class partition_proxy;

namespace detail {

class partition_impl;

//used to break cyclic dependency
struct convertable_to_partition_proxy
{
  convertable_to_partition_proxy( partition_impl const* arg_partition )
    : m_partition(arg_partition)
  {}

  partition_impl const* m_partition;
};

}} //namespace samba::detail

#endif //SAMBA_SAMBA_MESH_CONVERTABLE_TO_PARTITION_PROXY_HPP

