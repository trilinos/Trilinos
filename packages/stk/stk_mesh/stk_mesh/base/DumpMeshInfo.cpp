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

#include "stk_mesh/base/DumpMeshInfo.hpp"
#include "stk_mesh/base/BulkData.hpp"
#include "stk_mesh/base/MetaData.hpp"
#include "stk_mesh/baseImpl/Partition.hpp"
#include "stk_mesh/baseImpl/PrintEntityState.hpp"
namespace stk::mesh::impl {

std::ostream & print_entity_id(std::ostream & os, const MetaData & meta_data, EntityRank type, EntityId id)
{
  const std::string & name = meta_data.entity_rank_name(type);
  return os << name << "[" << id << "]";
}


std::ostream & print_entity_key(std::ostream & os, const MetaData & meta_data, const EntityKey & key)
{
  const EntityRank rank = key.rank();
  const EntityId id = key.id();
  return print_entity_id(os, meta_data, rank, id);
}

std::string print_entity_key(const MetaData & meta_data, const EntityKey & key)
{
  std::ostringstream out;
  print_entity_key(out, meta_data, key);
  return out.str();
}

bool print_comm_data_for_entity_in_ghosting(const BulkData & mesh, const Ghosting & ghosting,
                                            const EntityKey & entityKey, std::ostream & out)
{
  std::vector<int> procs;
  mesh.comm_procs(ghosting, entityKey, procs);

  if (procs.empty()) {
    return false;
  }

  out << "        Ghosting " << mesh.ghosting_part(ghosting).name() << " procs: "
      << stk::util::join(procs, ", ") << std::endl;

  return true;
}

void print_comm_data_for_entity(const BulkData & mesh, const EntityKey & entityKey, std::ostream & out)
{
  const std::vector<Ghosting*> & ghostLevels = mesh.ghostings();
  bool ghostedAnywhere = false;

  for (size_t ghostI = 0; ghostI < ghostLevels.size(); ++ghostI) {
    const Ghosting & ghosting = *ghostLevels[ghostI];
    ghostedAnywhere |= print_comm_data_for_entity_in_ghosting(mesh, ghosting, entityKey, out);
  }

  if (!ghostedAnywhere) {
    out << "        Not Communicated" << std::endl;
  }
}

void print_field_data_for_entity(const BulkData & mesh, const MeshIndex & meshIndex, std::ostream & out)
{
  const Bucket* bucket = meshIndex.bucket;
  const FieldVector & all_fields = mesh.mesh_meta_data().get_fields();

  for (FieldBase * field : all_fields) {
    if (field->entity_rank() != bucket->entity_rank()) continue;
    if (field_bytes_per_entity(*field, *bucket) > 0) {
      out << "        " << *field << ", ";
      field->print_data(out, meshIndex);
      out << std::endl;
    }
  }
}

void print_field_data_for_entity(const stk::mesh::BulkData & mesh, const stk::mesh::Entity & entity, std::ostream & out)
{
  print_field_data_for_entity(mesh, mesh.mesh_index(entity), out);
}

void print_entity_connectivity(const BulkData & mesh, const MeshIndex & meshIndex, std::ostream & out)
{
  const Bucket * bucket = meshIndex.bucket;
  size_t b_ord = meshIndex.bucket_ordinal;
  const std::vector<std::string> & rank_names = mesh.mesh_meta_data().entity_rank_names();
  EntityRank b_rank = bucket->entity_rank();
  for (EntityRank r = stk::topology::NODE_RANK, re = static_cast<EntityRank>(rank_names.size()); r < re; ++r) {
    out << "        Conn to " << rank_names[r] << std::endl;
    Entity const * entities = bucket->begin(b_ord, r);
    ConnectivityOrdinal const * ordinals = bucket->begin_ordinals(b_ord, r);
    const int num_conn = bucket->num_connectivity(b_ord, r);
    for (int c_itr = 0; c_itr < num_conn; ++c_itr) {
      Entity target_entity = entities[c_itr];
      uint32_t ord = ordinals[c_itr];
      out << "          [" << ord << "] " << mesh.entity_key(target_entity) << " ";
      if (r != stk::topology::NODE_RANK) {
        out << mesh.bucket(target_entity).topology();
        if (b_rank != stk::topology::NODE_RANK) {
          Permutation const *permutations = bucket->begin_permutations(b_ord, r);
          if (permutations) {
            out << " perm " << static_cast<unsigned>(permutations[c_itr]);
          }
        }
      }
      out << ", state = " << mesh.state(target_entity);
      out << std::endl;
    }
  }
}

void print_bucket_parts(const BulkData & /*mesh*/, const Bucket * bucket, std::ostream & out)
{
  out << "    bucket " << bucket->bucket_id() << " parts: { ";
  const PartVector & supersets = bucket->supersets();
  for (const Part * part : supersets) {
    out << part->name() << " ";
  }
  out << "}" << std::endl;
}

void print_entity_offset_and_state(const BulkData & mesh, const MeshIndex & meshIndex, std::ostream & out)
{
  Entity entity = (*meshIndex.bucket)[meshIndex.bucket_ordinal];
  out << "      " << print_entity_key(mesh.mesh_meta_data(), mesh.entity_key(entity)) << "(offset: "
      << entity.local_offset() << ", local_id: " << mesh.local_id(entity) << "), state = "
      << mesh.state(entity) << std::endl;

}

void print_connectivity_of_rank(const BulkData & bulk, const Entity & targetEntity, EntityRank connectedRank,
                                std::ostream & out)
{
  if (bulk.is_valid(targetEntity)) {
    out << connectedRank << "-connectivity(";

    const Entity * beginEntity = bulk.begin(targetEntity, connectedRank);
    const Entity * endEntity = bulk.end(targetEntity, connectedRank);
    std::for_each(beginEntity, endEntity,
                  [&](Entity entity) {
      unsigned numConnected = bulk.num_connectivity(entity, connectedRank);
      for (unsigned i = 0; i < numConnected; ++i) {
        if (bulk.is_valid(entity)) {
          out << "{" << bulk.identifier(entity) << ",topo=" << bulk.bucket(entity).topology() << ",owned="
              << bulk.bucket(entity).owned() << ",shared=" << bulk.bucket(entity).shared() << ",in_aura="
              << bulk.bucket(entity).in_aura() << ",custom-recv-ghost="
              << bulk.in_receive_custom_ghost(bulk.entity_key(entity)) << "}";
        }
        else {
          out << "{invalid entity!}";
        }
      }
      out << "), ";
    });
  }
  else {
    out << "invalid entity!";
  }
}

void dump_all_meta_info(const MetaData & meta, std::ostream & out)
{
  out << "MetaData info...(ptr=" << &meta << ")\n";

  out << "spatial dimension = " << meta.spatial_dimension() << "\n";

  out << "  Entity rank names:\n";
  const auto & entityRankNames = meta.entity_rank_names();
  for (size_t i = 0, e = entityRankNames.size(); i != e; ++i) {
    out << "    " << i << ": " << entityRankNames[i] << "\n";
  }
  out << "  Special Parts:\n";
  out << "    Universal part ord = " << meta.universal_part().mesh_meta_data_ordinal() << "\n";
  out << "    Owns part ord = " << meta.locally_owned_part().mesh_meta_data_ordinal() << "\n";
  out << "    Shared part ord = " << meta.globally_shared_part().mesh_meta_data_ordinal() << "\n";

  out << "  All parts:\n";
  const PartVector& all_parts = meta.get_parts();
  for(const Part* part : all_parts) {
    print(out, "    ", *part);
  }

  out << "  All fields:\n";
  const FieldVector& all_fields = meta.get_fields();
  for(const FieldBase* field : all_fields) {
     print(out, "    ", *field);
  }
}

void dump_mesh_bucket_info(const BulkData & bulk, std::ostream & out, Bucket * bucket)
{
  print_bucket_parts(bulk, bucket, out);

  for (size_t b_ord = 0, b_end = bucket->size(); b_ord < b_end; ++b_ord) {
    MeshIndex meshIndex(bucket, b_ord);
    impl::print_entity_offset_and_state(bulk, meshIndex, out);
    impl::print_entity_connectivity(bulk, meshIndex, out);
    const unsigned numFields = bulk.mesh_meta_data().get_fields().size();
    if (numFields > 0) {
      print_field_data_for_entity(bulk, meshIndex, out);
    }
    print_comm_data_for_entity(bulk, bulk.entity_key((*bucket)[b_ord]), out);
  }
}

void dump_all_mesh_info(const BulkData & bulk, std::ostream & out)
{
  const MetaData & meta = bulk.mesh_meta_data();

  dump_all_meta_info(meta, out);

  out << "BulkData (ptr=" << &bulk << ")\n";

  // Iterate all buckets for all ranks...
  const std::vector<std::string> & rank_names = meta.entity_rank_names();
  for (size_t i = 0, e = rank_names.size(); i < e; ++i) {
    EntityRank rank = static_cast<EntityRank>(i);
    out << "  All " << rank_names[i] << " entities:" << std::endl;

    const BucketVector& buckets = bulk.buckets(rank);
    for (Bucket * bucket : buckets) {
      dump_mesh_bucket_info(bulk, out, bucket);
    }
  }
}

void dump_mesh_per_proc(const BulkData & bulk, const std::string & fileNamePrefix)
{
  std::ostringstream oss;
  oss << fileNamePrefix << "." << bulk.parallel_rank();
  std::ofstream out(oss.str());
  dump_all_mesh_info(bulk, out);
}

class DumpInfoBulkData : public BulkData
{
public:
  const stk::mesh::impl::BucketRepository& get_bucket_repository() const { return m_bucket_repository; }
};

void dump_partition_summary(const BulkData & bulk, std::ostream & out)
{
  const auto & dumpInfoBulk = static_cast<const DumpInfoBulkData&>(bulk);

  const MetaData & meta = bulk.mesh_meta_data();
  const std::vector<std::string> & rankNames = meta.entity_rank_names();
  for (size_t irank = 0; irank < rankNames.size(); ++irank) {
    EntityRank rank = static_cast<EntityRank>(irank);
    out << "All " << rankNames[irank] << " partitions:" << std::endl;

    unsigned partitionIndex = 0;
    for (const stk::mesh::impl::Partition * partition : dumpInfoBulk.get_bucket_repository().get_partitions(rank)) {
      out << "Partition " << partitionIndex++ << ": [ ";
      for (const PartOrdinal partOrdinal : partition->get_legacy_partition_id()) {
        out << meta.get_part(partOrdinal).name() << " ";
      }
      out << "]" << std::endl;

      for (const Bucket * bucket : *partition) {
        out << "  Bucket " << bucket->bucket_id() << ": capacity=" << bucket->capacity()
            << " size=" << bucket->size() << " [ ";
        for (const Entity & entity : *bucket) {
          out << bulk.identifier(entity) << " ";
        }
        out << "]" << std::endl;
      }
    }
    out << std::endl;
  }
}

void dump_partition_summary_per_proc(const BulkData & bulk, const std::string& fileNamePrefix)
{
  std::ostringstream oss;
  oss << fileNamePrefix << "." << bulk.parallel_rank();
  std::ofstream out(oss.str());
  dump_partition_summary(bulk, out);
}

std::vector<int> get_bucket_size_histogram(const BulkData & bulk, unsigned binWidth)
{
  STK_ThrowRequire(binWidth > 0);
  STK_ThrowRequire(binWidth <= bulk.get_maximum_bucket_capacity());

  std::vector<int> sizeHistogram(bulk.get_maximum_bucket_capacity()/binWidth, 0);

  for (int irank = 0; irank < bulk.mesh_meta_data().entity_rank_count(); ++irank) {
    EntityRank rank = static_cast<EntityRank>(irank);
    for (const Bucket * bucket : bulk.buckets(rank)) {
      int bucketSize = std::min<int>(std::max<int>(bucket->size(), 1), bulk.get_maximum_bucket_capacity());
      int binIndex = (bucketSize - 1) / binWidth;
      sizeHistogram[binIndex] += 1;
    }
  }

  return sizeHistogram;
}

void dump_bucket_size_histogram(const BulkData & bulk, std::ostream & out, unsigned binWidth)
{
  const std::vector<int> sizeHistogram = get_bucket_size_histogram(bulk, binWidth);

  out << "size  count" << std::endl;
  int size = binWidth;
  for (int count : sizeHistogram) {
    out << std::left << std::setw(3) << size << "  " << count << std::endl;
    size += binWidth;
  }
}

void dump_bucket_size_histogram_per_proc(const BulkData & bulk, const std::string & fileNamePrefix, unsigned binWidth)
{
  std::ostringstream oss;
  oss << fileNamePrefix << "." << bulk.parallel_rank();
  std::ofstream out(oss.str());
  dump_bucket_size_histogram(bulk, out, binWidth);
}

void dump_global_bucket_size_histogram(const BulkData & bulk, const std::string & fileName, unsigned binWidth)
{
  std::vector<int> sizeHistogram = get_bucket_size_histogram(bulk, binWidth);

  if (bulk.parallel_rank() != 0) {
    MPI_Reduce(sizeHistogram.data(), sizeHistogram.data(), sizeHistogram.size(), MPI_INT, MPI_SUM, 0, bulk.parallel());
  }
  else {
    MPI_Reduce(MPI_IN_PLACE, sizeHistogram.data(), sizeHistogram.size(), MPI_INT, MPI_SUM, 0, bulk.parallel());
    std::ofstream out(fileName);
    out << "size  count" << std::endl;
    int size = binWidth;
    for (int count : sizeHistogram) {
      out << std::left << std::setw(3) << size << "  " << count << std::endl;
      size += binWidth;
    }
  }
}

}
