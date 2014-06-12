#ifndef SAMBA_SAMBA_DISTRIBUTED_DISTRIBUTED_MESH_HPP
#define SAMBA_SAMBA_DISTRIBUTED_DISTRIBUTED_MESH_HPP

#ifndef SAMBA_ENABLE_PARALLEL
#define SAMBA_ENABLE_PARALLEL true
#endif

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
class distributed_mesh : public query_interface<distributed_mesh>, public modification_interface<distributed_mesh>
{
  typedef boost::shared_ptr<detail::mesh_impl> mesh_handle;

public:

  typedef std::vector<entity_key> entity_key_vector;
  typedef std::vector<std::pair<process_id, std::vector<entity_key> *> > comm_list;

  //*************************************************************************
  //constructor
  //*************************************************************************

  /**
   * Construct a mesh with the default implementation.
   *   @param proc_id processor ID (e.g., MPI rank) for this mesh section
   *   @param num_procs total number of procs (1 greater than max MPI rank)
   *   @param arg_connectivity_map configures the connectivity setup for the mesh.
   */
  distributed_mesh(samba::process_id proc_id, size_t num_procs,
             samba::connectivity_map const& arg_connectivity_map = samba::connectivity_map::default_map() )
    : m_mesh_impl(new detail::mesh_impl( arg_connectivity_map, proc_id ) )
    , m_mesh_handle(m_mesh_impl)
    , m_owned_shared(new std::vector<std::vector<entity_key> >)
  {
    m_owned_shared->resize(num_procs);
  }

  //*************************************************************************
  //signals
  //*************************************************************************

  mesh_signals & signals()
  { return m_mesh_impl->signals(); }

  mesh_signals const& signals() const
  { return m_mesh_impl->signals(); }

  //*************************************************************************
  //Queries
  //*************************************************************************

  void get_comm_list(comm_list &comm_list_out);

  //*************************************************************************
  //modification
  //*************************************************************************

  template <typename EntityKeyIterator, typename EntitySetIterator>
  bool add_unowned_entities( EntityKeyIterator first_key
                             ,EntityKeyIterator last_key
                             ,const entity_state which_shared_type  // either ghosted or shared
                             ,EntitySetIterator first_add_set
                             ,EntitySetIterator last_add_set)
  {
    // WRITE ME!
    return false;
  }

  // Tell the mesh that an entity it owns is being shared or ghosted by another proc.
  bool add_sharer(process_id sharer_or_ghoster_id, entity_key owned);

private:

  detail::mesh_impl * m_mesh_impl;
  mesh_handle         m_mesh_handle;

  boost::shared_ptr<std::vector<std::vector<entity_key> > > m_owned_shared;

  friend class query_interface<distributed_mesh>;
  friend class modification_interface<distributed_mesh>;
};


inline
void distributed_mesh::get_comm_list(comm_list &comm_list_out)
{
  comm_list_out.clear();
  std::vector<std::vector<entity_key> > &comm_vec_vec = *m_owned_shared;
  for (size_t i = 0, end_i = comm_vec_vec.size(); i < end_i; ++i)
  {
    std::vector<entity_key> &ents = comm_vec_vec[i];
    if (ents.empty())
    {
      continue;
    }
    comm_list_out.push_back(std::pair<process_id, std::vector<entity_key> *>(process_id::create(i), &ents));
  }
}

inline
std::ostream& operator<<(std::ostream& os, distributed_mesh mesh_arg)
{
  return mesh_arg.detail_get_raw_mesh_impl()->streamit(os);
}

} //namespace samba

#include <samba/mesh/query_interface.tcc>
#include <samba/mesh/modification_interface.tcc>

namespace samba {

inline
bool distributed_mesh::add_sharer(process_id sharer_or_ghoster_id, entity_key owned)
{
  BOOST_ASSERT_MSG(static_cast<long long>(m_owned_shared->size()) > sharer_or_ghoster_id(),
                   "distributed_mesh::add_sharer --- sharer_or_ghoster_id beyond range of possible sharers");

  if ((process() == sharer_or_ghoster_id) || (owned.process() != process()))
  {
    return false;
  }
  std::vector<entity_key> &shared_with = (*m_owned_shared)[sharer_or_ghoster_id()];
  if (find (shared_with.begin(), shared_with.end(), owned) == shared_with.end())
  {
    shared_with.push_back(owned);
  }

  // THERE MORE TO WRITE!

  return false;
}

} // namespace samba

#endif  // SAMBA_SAMBA_DISTRIBUTED_DISTRIBUTED_MESH_HPP
