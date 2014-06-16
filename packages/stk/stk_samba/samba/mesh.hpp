#ifndef SAMBA_SAMBA_MESH_HPP
#define SAMBA_SAMBA_MESH_HPP

#include <samba/utility/debug_message.hpp>
#include <samba/mesh/mesh_ordered_includes.hpp>
#include <samba/mesh/query_interface.hpp>
#include <samba/mesh/modification_interface.hpp>

#include <boost/shared_ptr.hpp>

namespace samba {

/**
 * The interface class for this package. Clients will interact with this class
 * to perform mesh operations.
 */
class mesh : public query_interface<mesh>, public modification_interface<mesh>
{
  typedef boost::shared_ptr<detail::mesh_impl> mesh_handle;

public:

  typedef std::vector<entity_key> entity_key_vector;

  //*************************************************************************
  //constructor
  //*************************************************************************

  /**
   * Construct a mesh with the default implementation.
   *
   * arg_connectivity_map configures the connectivity setup for the mesh.
   */
  mesh( samba::connectivity_map const& arg_connectivity_map = samba::connectivity_map::default_map() )
    : m_mesh_impl(new detail::mesh_impl( arg_connectivity_map ) )
    , m_mesh_handle(m_mesh_impl)
  {}


  //// stk_samba, nalus, and their tests all build fine without this.  It seems
  //// to be dead code.
  ////   - Is it for future multithreaded work?
  ////   - There should be a comment about why m_mesh_handle() is called without
  ////       m_mesh_impl as an argument.
  ////
  //called by mesh_impl to convert itself to a mesh
  mesh( detail::mesh_impl & arg_impl)
    : m_mesh_impl(&arg_impl )
    , m_mesh_handle()
  {}

  //*************************************************************************
  //signals
  //*************************************************************************

  mesh_signals & signals()
  { return m_mesh_impl->signals(); }

  mesh_signals const& signals() const
  { return m_mesh_impl->signals(); }

private:

  detail::mesh_impl * m_mesh_impl;
  mesh_handle         m_mesh_handle;

  friend class query_interface<mesh>;
  friend class modification_interface<mesh>;
};

inline
std::ostream& operator<<(std::ostream& os, mesh mesh_arg)
{
  return mesh_arg.detail_get_raw_mesh_impl()->streamit(os);
}

} //namespace samba

#include <samba/mesh/query_interface.tcc>
#include <samba/mesh/modification_interface.tcc>

#endif //SAMBA_SAMBA_MESH_HPP
