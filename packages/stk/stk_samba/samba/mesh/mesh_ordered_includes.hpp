#ifndef SAMBA_SAMBA_MESH_MESH_ORDERED_INCLUDES_HPP
#define SAMBA_SAMBA_MESH_MESH_ORDERED_INCLUDES_HPP

//*****************************************************************************
//include order is chosen to break cycles and
//ensure correct inlining.
//
//DO NOT changing this order
//(unless the cycles changed)
//*****************************************************************************

#include <samba/types.hpp>
#include <samba/mesh/mesh_impl.hpp>

#include <samba/mesh/partition_connectivity.hpp>
#include <samba/mesh/fixed_partition_connectivity.tcc>
#include <samba/mesh/dynamic_partition_connectivity.tcc>

#include <samba/mesh/partition_impl.hpp>
#include <samba/mesh/entity_proxy.hpp>
#include <samba/mesh/partition_proxy.hpp>

#include <samba/set_expression.hpp>

#include <boost/bind.hpp>

#include <samba/mesh/query_impl.tcc>
#include <samba/mesh/conversion_impl.tcc>
#include <samba/mesh/mesh_impl_validate.tcc>
#include <samba/mesh/mesh_impl.tcc>

#include <samba/mesh/partition_impl.tcc>
#include <samba/mesh/entity_proxy.tcc>
#include <samba/mesh/partition_proxy.tcc>
#include <samba/mesh/partition_iterator.tcc>

#endif
