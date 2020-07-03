// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include <stk_mesh/baseImpl/MeshPrintUtils.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/diag/StringUtil.hpp>
//----------------------------------------------------------------------

namespace stk {
namespace mesh {
namespace impl {

std::ostream &
print_entity_id( std::ostream & os , const MetaData & meta_data ,
                  EntityRank type , EntityId id )
{
  const std::string & name = meta_data.entity_rank_name( type );
  return os << name << "[" << id << "]" ;
}


std::ostream &
print_entity_key( std::ostream & os , const MetaData & meta_data ,
                  const EntityKey & key )
{
  const EntityRank rank   = key.rank();
  const EntityId   id     = key.id();
  return print_entity_id( os , meta_data , rank , id );
}

std::string
print_entity_key( const MetaData & meta_data , const EntityKey & key )
{
  std::ostringstream out;
  print_entity_key(out, meta_data, key);
  return out.str();
}

bool print_comm_data_for_entity_in_ghosting(const BulkData& mesh, const Ghosting& ghosting, EntityKey entityKey, std::ostream& out)
{
    std::vector<int> procs;
    mesh.comm_procs(ghosting, entityKey, procs);
    if (procs.empty()) { return false; }
    out << "        Ghosting " << mesh.ghosting_part(ghosting).name() << " with procs:  " << stk::util::join(procs, ", ") << std::endl;
    return true;
}

void print_comm_data_for_entity(const BulkData& mesh, EntityKey entityKey, std::ostream& out) 
{
    const std::vector<Ghosting*> & ghostLevels = mesh.ghostings();
    bool ghostedAnywhere = false;
    for (size_t ghostI=0 ; ghostI<ghostLevels.size() ; ++ghostI) {
        const Ghosting& ghosting = *ghostLevels[ghostI];
        ghostedAnywhere |= print_comm_data_for_entity_in_ghosting(mesh, ghosting, entityKey, out);
    }
    if (!ghostedAnywhere) {
        out << "        Not Communicated" << std::endl;
    }
}

void print_field_data_for_entity(const BulkData& mesh, const MeshIndex& meshIndex, std::ostream& out)
{
    const Bucket* bucket = meshIndex.bucket;
    size_t b_ord = meshIndex.bucket_ordinal;
    const FieldVector& all_fields = mesh.mesh_meta_data().get_fields();
    for(FieldBase* field : all_fields) { 
        if(static_cast<unsigned>(field->entity_rank()) != bucket->entity_rank()) continue;
        FieldMetaData field_meta_data = field->get_meta_data_for_field()[bucket->bucket_id()];
        unsigned data_size = field_meta_data.m_bytes_per_entity;
        if (data_size > 0) { // entity has this field?
            void* data = field_meta_data.m_data + field_meta_data.m_bytes_per_entity * b_ord;
            out << "        For field: " << *field << ", has data: ";
            field->print_data(out, data, data_size); 
            out << std::endl;
        }
    }
}

void print_field_data_for_entity(const stk::mesh::BulkData& mesh, const stk::mesh::Entity entity, std::ostream& out)
{
  const stk::mesh::MeshIndex meshIndex = mesh.mesh_index(entity);
  print_field_data_for_entity(mesh, meshIndex, out);
}

void print_entity_connectivity(const BulkData& mesh, const MeshIndex& meshIndex, std::ostream& out)
{
    const Bucket* bucket = meshIndex.bucket;
    size_t b_ord = meshIndex.bucket_ordinal;
    const std::vector<std::string> & rank_names = mesh.mesh_meta_data().entity_rank_names();
    EntityRank b_rank = bucket->entity_rank();
    for (EntityRank r = stk::topology::NODE_RANK, re = static_cast<EntityRank>(rank_names.size()); r < re; ++r) {
        out << "        Connectivity to " << rank_names[r] << std::endl;
        Entity const* entities = bucket->begin(b_ord, r);
        ConnectivityOrdinal const* ordinals = bucket->begin_ordinals(b_ord, r);
        const int num_conn         = bucket->num_connectivity(b_ord, r);
        for (int c_itr = 0; c_itr < num_conn; ++c_itr) {
            Entity target_entity = entities[c_itr];
            uint32_t ord = ordinals[c_itr];
            out << "          [" << ord << "]  " << mesh.entity_key(target_entity) << "  ";
            if (r != stk::topology::NODE_RANK) {
                out << mesh.bucket(target_entity).topology();
                if (b_rank != stk::topology::NODE_RANK) {
                    Permutation const *permutations = bucket->begin_permutations(b_ord, r);
                    if (permutations) {
                        out << " permutation index " << permutations[c_itr];
                    }
                }
            }
            out << ", state = " << mesh.state(target_entity);
            out << std::endl;
        }
    }
}

void print_bucket_parts(const BulkData& mesh, const Bucket* bucket, std::ostream& out)
{
  out << "    Found bucket " << bucket->bucket_id() << " with superset parts: { ";
  const PartVector& supersets = bucket->supersets();
  for(const Part* part : supersets) {
    out << part->name() << " ";
  }    
  out << "}" << std::endl;
}

void print_entity_offset_and_state(const BulkData& mesh, const MeshIndex& meshIndex, std::ostream& out)
{
    Entity entity = (*meshIndex.bucket)[meshIndex.bucket_ordinal];
    out << "      " << print_entity_key(mesh.mesh_meta_data(), mesh.entity_key(entity)) << "(offset: " << entity.local_offset() <<
            "), state = " << mesh.state(entity) << std::endl;
}

} // namespace impl
} // namespace mesh
} // namespace stk

//----------------------------------------------------------------------
//----------------------------------------------------------------------

