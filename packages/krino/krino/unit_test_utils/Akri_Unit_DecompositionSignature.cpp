// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_Unit_DecompositionSignature.hpp>

#include <algorithm>
#include <iterator>
#include <sstream>
#include <tuple>

#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_util/parallel/CommSparse.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace krino {

bool ElementSignature::operator<(const ElementSignature & other) const
{
  return std::tie(vertexCoordinates, partNames) < std::tie(other.vertexCoordinates, other.partNames);
}

bool ElementSignature::operator==(const ElementSignature & other) const
{
  return vertexCoordinates == other.vertexCoordinates && partNames == other.partNames;
}

static ElementSignature element_signature(const stk::mesh::BulkData & mesh, const stk::mesh::Entity element)
{
  ElementSignature signature;
  for (const stk::mesh::Part * part : mesh.bucket(element).supersets())
  {
    // Leave out stk's own parts (such as the universal and locally owned parts), whose names are in
    // braces, and parts of other ranks.
    if (part->primary_entity_rank() == stk::topology::ELEMENT_RANK && part->name().front() != '{')
      signature.partNames.push_back(part->name());
  }
  std::sort(signature.partNames.begin(), signature.partNames.end());

  const stk::mesh::FieldBase & coordsField = *mesh.mesh_meta_data().coordinate_field();
  const auto coordsData = coordsField.data<double, stk::mesh::ReadOnly>();
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  const unsigned numVertices = mesh.bucket(element).topology().base().num_nodes();
  const stk::mesh::Entity * nodes = mesh.begin_nodes(element);
  for (unsigned n = 0; n < numVertices; ++n)
  {
    const auto x = coordsData.entity_values(nodes[n]);
    std::array<float, 3> vertex{{0.f, 0.f, 0.f}};
    for (unsigned d = 0; d < dim; ++d)
      vertex[d] = static_cast<float>(x(stk::mesh::ComponentIdx(d)));
    signature.vertexCoordinates.push_back(vertex);
  }
  std::sort(signature.vertexCoordinates.begin(), signature.vertexCoordinates.end());
  return signature;
}

static void pack_signature(stk::CommBuffer & buffer, const ElementSignature & signature)
{
  buffer.pack<size_t>(signature.partNames.size());
  for (const std::string & name : signature.partNames)
  {
    buffer.pack<size_t>(name.size());
    for (const char c : name)
      buffer.pack<char>(c);
  }
  buffer.pack<size_t>(signature.vertexCoordinates.size());
  for (const auto & vertex : signature.vertexCoordinates)
    for (const float x : vertex)
      buffer.pack<float>(x);
}

static ElementSignature unpack_signature(stk::CommBuffer & buffer)
{
  ElementSignature signature;
  size_t numParts = 0;
  buffer.unpack<size_t>(numParts);
  signature.partNames.resize(numParts);
  for (std::string & name : signature.partNames)
  {
    size_t length = 0;
    buffer.unpack<size_t>(length);
    name.resize(length);
    for (char & c : name)
      buffer.unpack<char>(c);
  }
  size_t numVertices = 0;
  buffer.unpack<size_t>(numVertices);
  signature.vertexCoordinates.resize(numVertices);
  for (auto & vertex : signature.vertexCoordinates)
    for (float & x : vertex)
      buffer.unpack<float>(x);
  return signature;
}

std::vector<ElementSignature> gather_decomposition_signature(
    const stk::mesh::BulkData & mesh, const stk::mesh::Selector & elementSelector)
{
  std::vector<ElementSignature> localSignatures;
  const stk::mesh::Selector ownedSelected = elementSelector & mesh.mesh_meta_data().locally_owned_part();
  for (const stk::mesh::Bucket * bucket : mesh.get_buckets(stk::topology::ELEMENT_RANK, ownedSelected))
    for (const stk::mesh::Entity element : *bucket)
      localSignatures.push_back(element_signature(mesh, element));

  std::vector<ElementSignature> signatures;
  if (0 == mesh.parallel_rank())
    signatures = localSignatures;

  stk::CommSparse commSparse(mesh.parallel());
  stk::pack_and_communicate(commSparse,
      [&]()
      {
        if (0 == mesh.parallel_rank()) return;
        stk::CommBuffer & buffer = commSparse.send_buffer(0);
        buffer.pack<size_t>(localSignatures.size());
        for (const ElementSignature & signature : localSignatures)
          pack_signature(buffer, signature);
      });
  stk::unpack_communications(commSparse,
      [&](int proc)
      {
        stk::CommBuffer & buffer = commSparse.recv_buffer(proc);
        size_t numSignatures = 0;
        buffer.unpack<size_t>(numSignatures);
        for (size_t i = 0; i < numSignatures; ++i)
          signatures.push_back(unpack_signature(buffer));
      });

  std::sort(signatures.begin(), signatures.end());
  return signatures;
}

size_t global_num_elements_differing_between_signatures(stk::ParallelMachine comm,
    const std::vector<ElementSignature> & signature0,
    const std::vector<ElementSignature> & signature1)
{
  std::vector<ElementSignature> differing;
  std::set_symmetric_difference(signature0.begin(),
      signature0.end(),
      signature1.begin(),
      signature1.end(),
      std::back_inserter(differing));
  const size_t localNumDiffering = differing.size();
  size_t globalNumDiffering = 0;
  stk::all_reduce_sum(comm, &localNumDiffering, &globalNumDiffering, 1);
  return globalNumDiffering;
}

std::string describe_elements_differing_between_signatures(const std::vector<ElementSignature> & signature0,
    const std::vector<ElementSignature> & signature1,
    const size_t maxNum)
{
  std::vector<ElementSignature> only0, only1;
  std::set_difference(signature0.begin(), signature0.end(), signature1.begin(), signature1.end(), std::back_inserter(only0));
  std::set_difference(signature1.begin(), signature1.end(), signature0.begin(), signature0.end(), std::back_inserter(only1));
  std::ostringstream os;
  size_t num = 0;
  for (const auto & [which, elements] : {std::make_pair(0, &only0), std::make_pair(1, &only1)})
  {
    for (const ElementSignature & element : *elements)
    {
      if (num++ == maxNum) return os.str();
      os << "  only in signature " << which << ":";
      for (const std::string & name : element.partNames)
        os << " " << name;
      for (const auto & x : element.vertexCoordinates)
        os << " (" << x[0] << "," << x[1] << "," << x[2] << ")";
      os << "\n";
    }
  }
  return os.str();
}

} // namespace krino
