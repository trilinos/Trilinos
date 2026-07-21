#include "GhostingCache.hpp"
#include "MetaData.hpp"
#include "BulkData.hpp"
#include "Ghosting.hpp"
#include "stk_util/util/string_case_compare.hpp"
#include "stk_util/parallel/ParallelComm.hpp"
#include "stk_util/parallel/ParallelReduce.hpp"
#include "stk_mesh/baseImpl/PartRepository.hpp"


namespace stk::mesh::impl {

GhostingCache::GhostingCache(BulkData& bulk_data) :
  m_meta_data(bulk_data.mesh_meta_data()),
  m_bulk_data(bulk_data)
{}

GhostingCache::~GhostingCache()
{
  for (auto& ghosting : m_ghosting) {
    delete ghosting;
  }

  for (auto& ghosting : m_deleted_ghostings) {
    delete ghosting;
  }
}

Ghosting* GhostingCache::create_ghosting(const std::string& name)
{
  STK_ThrowRequireMsg( !m_bulk_data.in_synchronized_state(),
                  "NOT in the ok-to-modify state" );

  check_name_debug(name);

  Ghosting* ghosting = get_ghosting_from_cache(name);
  if (ghosting) {
    return ghosting;
  } else {
    return create_new_ghosting(name);
  }

}

Ghosting* GhostingCache::get_ghosting_from_cache(const std::string& name)
{
  for(Ghosting* ghosting : m_ghosting) {
    if (ghosting->name() == name) {
      return ghosting;
    }
  }

  for (auto it = m_deleted_ghostings.begin(); it != m_deleted_ghostings.end(); ++it) {
    if ((*it)->name() == name) {
      Ghosting* ghosting = *it;
      m_ghosting.push_back(ghosting);
      m_deleted_ghostings.erase(it);
      return ghosting;
    }
  }

  return nullptr;
}

Ghosting* GhostingCache::create_new_ghosting(const std::string& name)
{
  Ghosting * const ghosting = new Ghosting( m_bulk_data , name , m_ghost_parts.size() );

  m_ghosting.push_back( ghosting );

  create_new_part(name);

  return ghosting;
}


void GhostingCache::create_new_part(const std::string& ghosting_name)
{
  if (m_ghost_parts.size() == 0) {
    STK_ThrowRequireMsg(equal_case(std::string("shared"), ghosting_name), "Expect shared to be the first ghosting created.");
    m_ghost_parts.push_back(&(m_meta_data.globally_shared_part()));
  }
  else if (m_ghost_parts.size() == 1) {
    STK_ThrowRequireMsg(equal_case(std::string("shared_aura"), ghosting_name), "Expect aura to be the second ghosting created.");
    Part & aura_part = m_meta_data.aura_part();
    aura_part.entity_membership_is_parallel_consistent(false);
    m_ghost_parts.push_back(&aura_part);
  }
  else {
    std::ostringstream oss;
    oss << "custom_ghosting_" << m_ghost_parts.size();
    std::string ghostPartName = stk::mesh::impl::convert_to_internal_name(oss.str());
    Part& ghost_part = m_meta_data.declare_part(ghostPartName);
    ghost_part.entity_membership_is_parallel_consistent(false);
    m_ghost_parts.push_back(&ghost_part);
  }
}

void GhostingCache::destroy_ghosting(Ghosting* ghosting)
{
  if(std::find(m_deleted_ghostings.begin(), m_deleted_ghostings.end(), ghosting) == m_deleted_ghostings.end()) {
    return;
  }

  for (auto it = m_ghosting.begin(); it != m_ghosting.end(); ++it) {
    if ((*it)->name() == ghosting->name())
    {
      m_deleted_ghostings.push_back(*it);
      m_ghosting.erase(it);
      break;
    }
  }
}


void GhostingCache::check_name_debug([[maybe_unused]] const std::string& name)
{
  #ifndef NDEBUG
  // Verify name is the same on all processors,
  // if not then throw an exception on all processors.
  if (m_bulk_data.parallel_size() > 1) {
    CommBroadcast bc( m_bulk_data.parallel() , 0 );

    if ( bc.parallel_rank() == 0 ) {
      bc.send_buffer().pack(name);
    }

    bc.allocate_buffer();

    if ( bc.parallel_rank() == 0 ) {
      bc.send_buffer().pack(name);
    }

    bc.communicate();

    std::string name_recv;
    bc.recv_buffer().unpack(name_recv);
    int error = name_recv != name;

    all_reduce( m_bulk_data.parallel() , ReduceMax<1>( & error ) );

    STK_ThrowRequireMsg(!error, std::string("Ghosting name parallel inconsistent: local name = ") + name + ", rank 0 name = " + name_recv);
  }
#endif
}

const std::vector<Ghosting*>& GhostingCache::get_ghostings() const
{
  return m_ghosting;
}

Part* GhostingCache::get_ghosting_part(const Ghosting* ghosting) const
{
  STK_ThrowRequireMsg(ghosting->ordinal() < m_ghost_parts.size(), "BulkData::ghosting_part ERROR, no part corresponds to ghosting with name="<<ghosting->name()<<" and ordinal="<<ghosting->ordinal());
  return m_ghost_parts[ghosting->ordinal()];
}

bool GhostingCache::is_ghost_part(Part* part) const
{
  return contains(m_ghost_parts, *part);
}


}
