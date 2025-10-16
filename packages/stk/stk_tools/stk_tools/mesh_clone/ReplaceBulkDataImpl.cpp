// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <stk_tools/mesh_clone/ReplaceBulkData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_io/IossBridge.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace stk {
namespace tools {
namespace impl {

void get_bucket_parts(const stk::mesh::Bucket & bucket, stk::mesh::PartVector & parts)
{
  parts.clear();

  stk::mesh::PartVector const& bucketParts = bucket.supersets();
  for ( stk::mesh::PartVector::const_iterator ip = bucketParts.begin(); ip != bucketParts.end(); ++ip ) {
    stk::mesh::Part & bucketPart = **ip;
    if (bucketPart.primary_entity_rank() != stk::topology::INVALID_RANK &&
        bucketPart.primary_entity_rank() != bucket.entity_rank()) {
      continue;
    }
    if (stk::mesh::is_auto_declared_part(bucketPart) && !stk::mesh::is_topology_root_part(bucketPart)) {
      continue;
    }

    parts.push_back(&bucketPart);
  }
  std::sort( parts.begin(), parts.end(), stk::mesh::PartLess() );
}

stk::mesh::Selector translate_selector(const stk::mesh::Selector & in_selector, const stk::mesh::MetaData & out_meta)
{
  if (in_selector == stk::mesh::Selector())
  {
    return in_selector;
  }
  STK_ThrowRequireMsg(in_selector.is_all_unions(), "Cannot translate selector " << in_selector);
  stk::mesh::PartVector in_parts, out_parts;
  in_selector.get_parts(in_parts);
  translate_parts(in_parts, out_meta, out_parts);
  return stk::mesh::selectUnion(out_parts);
}

void translate_parts(const stk::mesh::PartVector & inParts, const stk::mesh::MetaData & outMeta, stk::mesh::PartVector & outParts)
{
  outParts.clear();
  for ( auto&& inPart : inParts ) {
    stk::mesh::Part* outPart = outMeta.get_part(inPart->name());
    if(nullptr != outPart) {
      stk::mesh::insert(outParts, *outPart);
    }
  }
}

void
clone_meta_data_parts_and_fields(const stk::mesh::MetaData & in_meta, stk::mesh::MetaData & out_meta)
{

  // This is pretty nasty.  We want the part.mesh_meta_data_ordinal to be the same for the in_meta and out_meta.
  // To accomplish this, we must be careful about the order of the part creation.
  // Specifically, we must clone the parts that were created before in_meta.initialize() were called, then
  // call out_meta.initialize(), then clone the rest.

  const stk::mesh::PartVector & in_parts = in_meta.get_parts();

  unsigned ipart = 0;
  bool more_to_do = ipart < in_parts.size();
  while (more_to_do)
  {
    stk::mesh::Part * in_part = in_parts[ipart++];
    if (stk::mesh::is_topology_root_part(*in_part))
    {   
      more_to_do = false;
    }   
    else
    {   
      more_to_do = ipart < in_parts.size();
      STK_ThrowRequire(in_part->primary_entity_rank() == stk::topology::INVALID_RANK);
      stk::mesh::Part & out_part = out_meta.declare_part(in_part->name());
      STK_ThrowRequire(out_part.mesh_meta_data_ordinal() == in_part->mesh_meta_data_ordinal());
    }   
  }

  out_meta.initialize(in_meta.spatial_dimension(), in_meta.entity_rank_names());

  for ( auto&& in_part : in_parts )
  {
    stk::mesh::Part & out_part =
        (in_part->primary_entity_rank() == stk::topology::INVALID_RANK) ?
        out_meta.declare_part(in_part->name()) :
        out_meta.declare_part(in_part->name(), in_part->primary_entity_rank(), in_part->force_no_induce());
    STK_ThrowRequire(out_part.mesh_meta_data_ordinal() == in_part->mesh_meta_data_ordinal());
    if (stk::io::is_part_io_part(*in_part))
    {   
      stk::io::put_io_part_attribute(out_part);
    }         
  }

  for ( auto&& in_part : in_parts )
  {
    stk::mesh::Part & out_part = out_meta.get_part(in_part->mesh_meta_data_ordinal());
    const stk::mesh::PartVector & in_subsets = in_part->subsets();
    for (auto && in_subset : in_subsets)
    {   
      stk::mesh::Part & out_subset = out_meta.get_part(in_subset->mesh_meta_data_ordinal());
      out_meta.declare_part_subset(out_part, out_subset);
    }   
  }

  const stk::mesh::FieldVector & in_fields = in_meta.get_fields();
  for ( auto&& in_field : in_fields )
  {
    if (in_field->state() == stk::mesh::StateNone)
    {   
      stk::mesh::FieldBase * out_field = in_field->clone(out_meta.get_field_repository());

      for ( auto&& in_restriction : in_field->restrictions() )
      {
        const stk::mesh::Selector out_selector = translate_selector(in_restriction.selector(), out_meta);

        out_meta.declare_field_restriction( *out_field, out_selector, in_restriction.num_scalars_per_entity(), in_restriction.dimension() );
      }
    }
  }
}

void copy_field_data(const stk::mesh::BulkData & inMesh, stk::mesh::BulkData & outMesh,
                     const stk::mesh::FieldBase & inField, const stk::mesh::FieldBase & outField)
{
  inField.synchronize<stk::mesh::ReadOnly>();
  outField.synchronize<stk::mesh::ReadWrite>();

  stk::mesh::EntityRank entityRank = outField.entity_rank();
  STK_ThrowRequire(inField.entity_rank() == entityRank);

  const stk::mesh::MetaData & outMeta = outMesh.mesh_meta_data();
  stk::mesh::Selector outFieldSelector = stk::mesh::selectField(outField) & !outMeta.aura_part();
  const stk::mesh::BucketVector & outBuckets = outMesh.get_buckets(entityRank, outFieldSelector);

  auto copy_bytes = [&](auto& inEntityBytes, auto& outEntityBytes) {
    STK_ThrowRequireMsg(inEntityBytes.num_bytes() == outEntityBytes.num_bytes(),
                        "Mismatched field size for field " << inField.name() << " inLength = "
                        << inEntityBytes.num_bytes() << " outLength = " << outEntityBytes.num_bytes() << "\n");

    for (stk::mesh::ByteIdx byte : inEntityBytes.bytes()) {
      outEntityBytes(byte) = inEntityBytes(byte);
    }
  };

  auto inFieldBytes = inField.data_bytes<const std::byte>();
  auto outFieldBytes = outField.data_bytes<std::byte>();

  if (inField.host_data_layout() == stk::mesh::Layout::Right) {
    for (stk::mesh::Bucket* outBucket : outBuckets) {
      for (stk::mesh::Entity outEntity : *outBucket) {
        const stk::mesh::EntityId outEntityId = outMesh.identifier(outEntity);
        const stk::mesh::Entity inEntity = inMesh.get_entity(entityRank, outEntityId);
        STK_ThrowRequireMsg(inMesh.is_valid(inEntity), "Mismatched meshes: " << outMesh.entity_key(outEntity)
                            << " not found in input mesh.");
        auto inEntityBytes = inFieldBytes.entity_bytes<stk::mesh::Layout::Right>(inEntity);
        auto outEntityBytes = outFieldBytes.entity_bytes<stk::mesh::Layout::Right>(outEntity);
        copy_bytes(inEntityBytes, outEntityBytes);
      }
    }
  }
  else if (inField.host_data_layout() == stk::mesh::Layout::Left) {
    for (stk::mesh::Bucket* outBucket : outBuckets) {
      for (stk::mesh::Entity outEntity : *outBucket) {
        const stk::mesh::EntityId outEntityId = outMesh.identifier(outEntity);
        const stk::mesh::Entity inEntity = inMesh.get_entity(entityRank, outEntityId);
        STK_ThrowRequireMsg(inMesh.is_valid(inEntity), "Mismatched meshes: " << outMesh.entity_key(outEntity)
                            << " not found in input mesh.");
        auto inEntityBytes = inFieldBytes.entity_bytes<stk::mesh::Layout::Left>(inEntity);
        auto outEntityBytes = outFieldBytes.entity_bytes<stk::mesh::Layout::Left>(outEntity);
        copy_bytes(inEntityBytes, outEntityBytes);
      }
    }
  }
  else {
    STK_ThrowErrorMsg("Unsupported host Field data layout: " << inField.host_data_layout());
  }
}

void copy_field_data(const stk::mesh::BulkData & inMesh, stk::mesh::BulkData & outMesh)
{
  const stk::mesh::MetaData & inMeta = inMesh.mesh_meta_data();
  stk::mesh::MetaData & outMeta = outMesh.mesh_meta_data();

  const stk::mesh::FieldVector & inFields = inMeta.get_fields();
  const stk::mesh::FieldVector & outFields = outMeta.get_fields();
//  STK_ThrowAssert(inFields.size() == outFields.size());

  unsigned numFields = std::min(inFields.size(), outFields.size());
  for ( unsigned fieldIndex=0; fieldIndex < numFields; ++fieldIndex ) {
    const stk::mesh::FieldBase & inField = *inFields[fieldIndex];
    const stk::mesh::FieldBase & outField = *outFields[fieldIndex];
    copy_field_data(inMesh, outMesh, inField, outField);
  }

  const bool outMeshAuraFromCommunication = outMesh.is_automatic_aura_on() && !inMesh.is_automatic_aura_on();
  if (outMeshAuraFromCommunication) {
    const std::vector<stk::mesh::Ghosting *> ghostings = outMesh.ghostings();
    const std::vector<const stk::mesh::FieldBase *> const_fields(outFields.begin(), outFields.end());
    stk::mesh::communicate_field_data(*ghostings[stk::mesh::BulkData::AURA], const_fields);
  }
}

void copy_field_data(const stk::mesh::BulkData & inMesh, std::vector<stk::mesh::FieldBase*>& inFieldsToCopyFrom, stk::mesh::BulkData & outMesh, std::vector<stk::mesh::FieldBase*>& outFieldsToCopyTo )
{

//  STK_ThrowAssert(inFields.size() == outFields.size());

  unsigned numFields = std::min(inFieldsToCopyFrom.size(), outFieldsToCopyTo.size());
  for ( unsigned fieldIndex=0; fieldIndex < numFields; ++fieldIndex ) {
    const stk::mesh::FieldBase & inField = *inFieldsToCopyFrom[fieldIndex];
    const stk::mesh::FieldBase & outField = *outFieldsToCopyTo[fieldIndex];
    copy_field_data(inMesh, outMesh, inField, outField);
  }

  const bool outMeshAuraFromCommunication = outMesh.is_automatic_aura_on() && !inMesh.is_automatic_aura_on();
  if (outMeshAuraFromCommunication) {
    const std::vector<stk::mesh::Ghosting *> ghostings = outMesh.ghostings();
    const std::vector<const stk::mesh::FieldBase *> const_fields(outFieldsToCopyTo.begin(), outFieldsToCopyTo.end());
    stk::mesh::communicate_field_data(*ghostings[stk::mesh::BulkData::AURA], const_fields);
  }
}

void clone_bulk_data_entities(const stk::mesh::BulkData & inMesh, stk::mesh::BulkData & outMesh)
{
  const stk::mesh::MetaData & inMeta = inMesh.mesh_meta_data();
  const stk::mesh::MetaData & outMeta = outMesh.mesh_meta_data();

  std::vector<int> sharing;

  const stk::mesh::EntityRank highestEntityRank = static_cast<stk::mesh::EntityRank>(inMeta.entity_rank_count()-1);
  for (stk::mesh::EntityRank entityRank = stk::topology::NODE_RANK; entityRank <= highestEntityRank; ++entityRank) {
    stk::mesh::Selector notGhost = inMeta.locally_owned_part() | inMeta.globally_shared_part();
    const stk::mesh::BucketVector & buckets = inMesh.get_buckets(entityRank, notGhost);

    for ( auto&& bucketPtr : buckets ) {
      const stk::mesh::Bucket & b = *bucketPtr;

      stk::mesh::PartVector inBucketParts, outBucketParts;
      get_bucket_parts(b,inBucketParts);

      translate_parts(inBucketParts, outMeta, outBucketParts);

      const bool bucketIsLocallyOwned = b.member(inMeta.locally_owned_part());

      const int length = b.size();
      for (int i = 0; i < length; ++i) {
        stk::mesh::Entity inEntity = b[i];

        stk::mesh::Entity outEntity;

        // Unfortunately, there is a subtle bug that may be present in the input mesh that can cause a face to be present on
        // this processor even when there are no locally owned elements using that face.  We will skip these faces.
        if (!bucketIsLocallyOwned && entityRank == inMeta.side_rank()) {
          bool foundLocallyOwnedElem = false;
          const unsigned numElems = inMesh.num_elements(inEntity);
          const stk::mesh::Entity* elems = inMesh.begin_elements(inEntity);
          for (unsigned ielem=0; ielem<numElems; ++ielem) {
            if (inMesh.bucket(elems[ielem]).member(inMeta.locally_owned_part())) {
              foundLocallyOwnedElem = true;
              break;
            }
          }
          if (!foundLocallyOwnedElem) {
            continue;
          }
        }

        if (!outMesh.is_valid(outEntity)) {
          outEntity = outMesh.declare_entity( entityRank, inMesh.identifier(inEntity), outBucketParts );
        }

        if (stk::topology::NODE_RANK == entityRank && outMesh.parallel_size() > 1) {
          inMesh.comm_shared_procs(inEntity,sharing);
          for ( size_t k = 0 ; k < sharing.size() ; ++k ) {
            outMesh.add_node_sharing(outEntity, sharing[k]);
          }
        }

        for (stk::mesh::EntityRank relativeRank = stk::topology::NODE_RANK; relativeRank < entityRank; ++relativeRank) {
          const unsigned numRelatives = inMesh.num_connectivity(inEntity, relativeRank);
          const stk::mesh::Entity* inRelatives = inMesh.begin(inEntity, relativeRank);
          const stk::mesh::ConnectivityOrdinal * relativeOrdinals = inMesh.begin_ordinals(inEntity, relativeRank);
          const stk::mesh::Permutation * relativePermutations = inMesh.begin_permutations(inEntity, relativeRank);

          for (unsigned relativeIndex=0; relativeIndex<numRelatives; ++relativeIndex) {
            stk::mesh::Entity inRelative = inRelatives[relativeIndex];

            stk::mesh::Entity outRelative = outMesh.get_entity(relativeRank, inMesh.identifier(inRelative));
            if (!outMesh.is_valid(outRelative)) {
              // This relative might live only on another processor
              continue;
            }

            if (inMesh.has_permutation(inEntity, relativeRank)) {
              outMesh.declare_relation( outEntity, outRelative, relativeOrdinals[relativeIndex], relativePermutations[relativeIndex] );
            } else {
              outMesh.declare_relation( outEntity, outRelative, relativeOrdinals[relativeIndex] );
            }
          }
        }
      }
    }
  }
}

}
}
}


