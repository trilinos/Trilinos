/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stk_mesh/base/MemoryUsage.hpp>
#include <stk_mesh/base/GetEntities.hpp>

#include <iostream>

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
    mem_usage.part_bytes += sizeof(Part*)        * parts[i]->intersection_of().size();
    mem_usage.part_bytes += sizeof(PartRelation) * parts[i]->relations().size();
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

  std::vector<Entity*> entities;
  for(size_t i=0; i<nranks; ++i) {
    EntityRank rank = i;
    total_bytes += mem_usage.entity_counts[rank]*sizeof(Entity);

    get_entities(bulk, rank, entities);

    for(size_t n=0; n<entities.size(); ++n) {
      Entity& entity = *entities[n];
      for(EntityRank r=0; r<i; ++r) {
        unsigned num_rels = entity.relations(r).size();
        mem_usage.downward_relation_counts[r] += num_rels;
        total_bytes += num_rels*sizeof(Relation);
      }
      for(EntityRank r=i+1; r<nranks; ++r) {
        unsigned num_rels = entity.relations(r).size();
        mem_usage.upward_relation_counts[r] += num_rels;
        total_bytes += num_rels*sizeof(Relation);
      }
    }

    const std::vector<Bucket*>& buckets = bulk.buckets(rank);
    mem_usage.bucket_counts[rank] = buckets.size();
    for(size_t b=0; b<buckets.size(); ++b) {
      Bucket& bucket = *buckets[b];
      mem_usage.bucket_bytes[rank] += bucket.allocation_size();
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
  os << "Total bytes: "<<mem_usage.total_bytes<<" ("<<((double)mem_usage.total_bytes)/(1024*1024)<<"MB)"<<std::endl;
}

}//namespace mesh
}//namespace stk


