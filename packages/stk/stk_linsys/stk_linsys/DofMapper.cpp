/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <sstream>
#include <stdexcept>

#include <stk_linsys/DofMapper.hpp>

#include <stk_mesh/base/FieldData.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>

#include <stk_linsys/ImplDetails.hpp>

namespace stk {
namespace linsys {

DofMapper::DofMapper(MPI_Comm comm, bool create_reverse_mappings)
 : m_field_id_map(),
   m_fei_vecspace(new fei::VectorSpace(comm)),
   m_reverse_mappings_enabled(create_reverse_mappings),
   m_fei_reversemap(NULL)
{
}

DofMapper::~DofMapper()
{
  delete m_fei_reversemap;
}

namespace {

void throw_fei_err(const std::string& mesg, int err)
{
  std::ostringstream osstr;
  osstr << mesg << err;
  std::string str = osstr.str();
  throw std::runtime_error(str);
}

}//namespace <anonymous>

void
DofMapper::add_dof_mappings(const stk::mesh::BulkData& mesh_bulk,
                            const stk::mesh::Selector& selector,
                            stk::mesh::EntityRank ent_type,
                            const stk::mesh::FieldBase& field)
{
  int idType = static_cast<int>(ent_type);

  m_fei_vecspace->defineIDTypes( 1, &idType );

  int field_id = impl::map_field_to_int(m_field_id_map, field);

  const std::vector<stk::mesh::Bucket*>& all_buckets = mesh_bulk.buckets(ent_type);
  std::vector<stk::mesh::Bucket*> buckets;
  stk::mesh::get_buckets(selector, all_buckets, buckets);

  bool already_declared_fei_field = false;

  for(size_t i=0; i<buckets.size(); ++i) {
    int field_data_size = stk::mesh::field_data_size( field, *buckets[i] );

    if (field_data_size == 0) continue;

    if (!already_declared_fei_field) {
      int field_size = field.max_size(ent_type);
      m_fei_vecspace->defineFields(1, &field_id, &field_size);
      already_declared_fei_field = true;
    }

    std::vector<int> ids(mesh_bulk.bucket_capacity());

    stk::mesh::Bucket::iterator
      iter = buckets[i]->begin(), iter_end = buckets[i]->end();

    int num_ids = 0;

    for(; iter != iter_end; ++iter) {
      stk::mesh::Entity& entity = *iter;
      ids[num_ids++] = impl::entityid_to_int(entity.identifier());
    }

    int err = m_fei_vecspace->addDOFs(field_id, idType, num_ids, &ids[0]);
    if (err != 0) throw_fei_err("stk::linsys::DofMapper::add_dof_mappings ERROR: fei::VectorSpace::addDOFs returned error-code=",err);
  }

  std::vector<int> shared_ids;
  std::vector<int> sharing_procs;

  const std::vector<stk::mesh::Entity*>& entity_comm = mesh_bulk.entity_comm();
  for(size_t i=0; i<entity_comm.size(); ++i) {
    stk::mesh::Entity* ent = entity_comm[i] ;

    //we only care about entities of the right type, and which have data for 'field'.
    if (ent->entity_rank() != ent_type) continue;
    if (!stk::mesh::field_data_valid(field, *ent)) continue;

    const stk::mesh::PairIterEntityComm ec = ent->sharing();

    for ( size_t j = 0 ; j < ec.size() ; ++j ) {
      shared_ids.push_back(impl::entityid_to_int(ent->identifier()));
      sharing_procs.push_back(ec[j].proc);
    }
  }

  if (shared_ids.size() > 0) {
    std::vector<int> num_sharing_procs_per_id(shared_ids.size(), 1);

    int err = m_fei_vecspace->initSharedIDs(shared_ids.size(), idType,
                                           &shared_ids[0], &num_sharing_procs_per_id[0],
                                           &sharing_procs[0]);
    if (err != 0) throw_fei_err("stk::linsys::DofMapper::add_dof_mappings ERROR: fei::VectorSpace::initSharedIDs returned error-code=",err);
  }
}

void
DofMapper::finalize()
{
  int err = m_fei_vecspace->initComplete();
  if (err != 0) throw_fei_err("stk::linsys::DofMapper::finalize ERROR: fei::VectorSpace::initComplete returned error-code=",err);

  if (m_reverse_mappings_enabled) {
    delete m_fei_reversemap;
    m_fei_reversemap = new fei::ReverseMapper(*m_fei_vecspace);
  }
}

int
DofMapper::get_field_id(const stk::mesh::FieldBase& field) const
{
  return impl::query_field_to_int_mapping(m_field_id_map, field);
}

int
DofMapper::get_global_index(stk::mesh::EntityRank ent_type,
                            stk::mesh::EntityId ent_id,
                            stk::mesh::FieldBase& field,
                            int offset_into_field)
{
  int err = 0, index = 0;
  int field_id = get_field_id(field);
  int int_id = impl::entityid_to_int(ent_id);

  try {
    err = m_fei_vecspace->getGlobalIndex(ent_type, int_id, field_id, index);
    if (err != 0) throw_fei_err("fei::VectorSpace::getGlobalIndex error=",err);
  }
  catch (...) {
    std::ostringstream msg;
    msg << "stk::linsys::DofMapper::get_global_index ERROR: "
     << "fei::VectorSpace::getGlobalIndex returned error-code ("<<err
     << ") or threw exception, probably meaning that the entity with type="<<ent_type<<" and id="
     << ent_id<<" was not found.";
    std::string str = msg.str();
    throw std::runtime_error(str);
  }

  return index + offset_into_field;
}

void
DofMapper::get_dof(int global_index,
                   stk::mesh::EntityRank& ent_type,
                   stk::mesh::EntityId& ent_id,
                   const stk::mesh::FieldBase*& field,
                   int& offset_into_field) const
{
  if (!m_reverse_mappings_enabled || m_fei_reversemap == NULL) {
    std::ostringstream msg;
    msg << "stk::linsys::DofMapper::get_dof ERROR: "
     << "either reverse-mappings are disabled or DofMapper::finalize hasn't "
     << "been called yet.";
    std::string str = msg.str();
    throw std::runtime_error(str);
  }

  fei::EqnRecord eqrec = m_fei_reversemap->getEqnRecord(global_index);

  ent_type = eqrec.IDType;
  ent_id   = eqrec.ID;
  offset_into_field = eqrec.offset;
  field = impl::get_field(m_field_id_map, eqrec.fieldID);
}

}//namespace linsys
}//namespace stk

