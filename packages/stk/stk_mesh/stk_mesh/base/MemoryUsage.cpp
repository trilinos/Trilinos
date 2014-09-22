/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_mesh/base/MemoryUsage.hpp>
#include <iostream>                     // for operator<<, basic_ostream, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket
#include "stk_mesh/base/BulkData.hpp"   // for BulkData
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/FieldBase.hpp"  // for FieldBase, etc
#include "stk_mesh/base/FieldRestriction.hpp"  // for FieldRestriction
#include "stk_mesh/base/MetaData.hpp"   // for MetaData
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"   // for Relation
#include "stk_mesh/base/Selector.hpp"   // for Selector
#include "stk_mesh/base/Types.hpp"      // for PartVector, EntityRank
#include "stk_mesh/baseImpl/FieldRepository.hpp"  // for FieldVector
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/environment/ReportHandler.hpp"  // for ThrowErrorMsg


namespace stk {
namespace mesh {

void compute_memory_usage(const BulkData& bulk, MemoryUsage& mem_usage)
{
  mem_usage.entity_rank_names = bulk.mesh_meta_data().entity_rank_names();

  const FieldVector& fields = bulk.mesh_meta_data().get_fields();
  mem_usage.num_fields = fields.size();
  mem_usage.field_bytes = fields.size()*sizeof(FieldBase);
  for(size_t i=0; i<fields.size(); ++i) {
    mem_usage.field_bytes += fields[i]->name().length();
    mem_usage.field_bytes += sizeof(FieldRestriction)*fields[i]->restrictions().size();
  }

  const PartVector& parts = bulk.mesh_meta_data().get_parts();
  mem_usage.num_parts = parts.size();
  mem_usage.part_bytes = parts.size()*sizeof(Part);
  for(size_t i=0; i<parts.size(); ++i) {
    mem_usage.part_bytes += parts[i]->name().length();
    mem_usage.part_bytes += sizeof(Part*)        * parts[i]->supersets().size();
    mem_usage.part_bytes += sizeof(Part*)        * parts[i]->subsets().size();
  }

  size_t total_bytes = mem_usage.field_bytes + mem_usage.part_bytes;

  mem_usage.entity_counts.clear();
  mem_usage.downward_relation_counts.clear();
  mem_usage.upward_relation_counts.clear();
  mem_usage.bucket_counts.clear();
  mem_usage.bucket_bytes.clear();

  Selector all = bulk.mesh_meta_data().universal_part();
  count_entities(all, bulk, mem_usage.entity_counts);

  size_t nranks = mem_usage.entity_counts.size();
  mem_usage.downward_relation_counts.resize(nranks, 0);
  mem_usage.upward_relation_counts.resize(nranks, 0);
  mem_usage.bucket_counts.resize(nranks, 0);
  mem_usage.bucket_bytes.resize(nranks, 0);

  std::vector<Entity> entities;
  for(size_t i=0; i<nranks; ++i) {
    EntityRank rank_i = static_cast<EntityRank>(i);
    total_bytes += mem_usage.entity_counts[rank_i]*sizeof(Entity);

    get_entities(bulk, rank_i, entities);

    for(size_t n=0; n<entities.size(); ++n) {
      Entity entity = entities[n];
      for(EntityRank r=stk::topology::NODE_RANK; r<rank_i; ++r) {
        unsigned num_rels = bulk.num_connectivity(entity, r);
        mem_usage.downward_relation_counts[r] += num_rels;
        ThrowErrorMsg("stk::mesh::compute_memory_usage need to be largely re-written for the new Connectivity scheme but is not needed for this 4.27.7.");
      }
      for(EntityRank r=static_cast<EntityRank>(rank_i+1); r<nranks; ++r) {
        unsigned num_rels = bulk.num_connectivity(entity, r);
        mem_usage.upward_relation_counts[r] += num_rels;
        ThrowErrorMsg("stk::mesh::compute_memory_usage need to be largely re-written for the new Connectivity scheme but is not needed for this 4.27.7.");
      }
    }

    const BucketVector& buckets = bulk.buckets(rank_i);
    mem_usage.bucket_counts[rank_i] = buckets.size();
    for(size_t b=0; b<buckets.size(); ++b) {
      Bucket& bucket = *buckets[b];
      mem_usage.bucket_bytes[rank_i] += bucket.allocation_size();
      total_bytes += bucket.allocation_size();
    }
  }

  mem_usage.total_bytes = total_bytes;
}

void print_memory_usage(const MemoryUsage& mem_usage, std::ostream& os)
{
  os << "----- stk_mesh Memory Usage: ------"<<std::endl;
  os << "Fields:"<<std::endl;
  os << "  "<<mem_usage.num_fields<<" fields, "<<mem_usage.field_bytes<<" bytes"<<std::endl;
  os << "Parts:"<<std::endl;
  os << "  "<<mem_usage.num_parts<<" parts, "<<mem_usage.part_bytes<<" bytes"<<std::endl;
  os << "Entities:"<<std::endl;
  for(size_t i=0; i<mem_usage.entity_counts.size(); ++i) {
    int n = mem_usage.entity_counts[i];
    unsigned bytes = n*sizeof(Entity);
    if (mem_usage.entity_rank_names.size() > i)
      os << "  "<<mem_usage.entity_rank_names[i]<<": ";
    else
      os << "  Rank "<<i<<": ";
    os << n << " entities, "<< bytes<<" bytes"<<std::endl;
  }
  os << "Downward Relations:"<<std::endl;
  for(size_t i=0; i<mem_usage.downward_relation_counts.size(); ++i) {
    int n = mem_usage.downward_relation_counts[i];
    unsigned bytes = n*sizeof(Relation);
    if (mem_usage.entity_rank_names.size() > i)
      os << "  "<<mem_usage.entity_rank_names[i]<<": ";
    else
      os << "  Rank "<<i<<": ";
    os << n << " relations, "<< bytes<<" bytes"<<std::endl;
  }
  os << "Upward Relations:"<<std::endl;
  for(size_t i=0; i<mem_usage.upward_relation_counts.size(); ++i) {
    int n = mem_usage.upward_relation_counts[i];
    unsigned bytes = n*sizeof(Relation);
    if (mem_usage.entity_rank_names.size() > i)
      os << "  "<<mem_usage.entity_rank_names[i]<<": ";
    else
      os << "  Rank "<<i<<": ";
    os << n << " relations, "<< bytes<<" bytes"<<std::endl;
  }
  os << "Buckets:"<<std::endl;
  for(size_t i=0; i<mem_usage.bucket_counts.size(); ++i) {
    int n = mem_usage.bucket_counts[i];
    unsigned bytes = mem_usage.bucket_bytes[i];
    if (mem_usage.entity_rank_names.size() > i)
      os << "  "<<mem_usage.entity_rank_names[i]<<": ";
    else
      os << "  Rank "<<i<<": ";
    os << n << " buckets, "<< bytes<<" bytes"<<std::endl;
  }
  os << "Total bytes: "<<mem_usage.total_bytes<<" ("<<(static_cast<double>(mem_usage.total_bytes))/(1024*1024)<<"MB)"<<std::endl;
}

}//namespace mesh
}//namespace stk


