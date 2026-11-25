// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_UNIT_TEST_UTILS_AKRI_UNIT_MULTIBLOCKMESH_HPP_
#define KRINO_KRINO_UNIT_TEST_UTILS_AKRI_UNIT_MULTIBLOCKMESH_HPP_

#include <algorithm>
#include <memory>
#include <set>
#include <string>
#include <vector>

#include <Akri_BoundingBox.hpp>
#include <Akri_MeshSurface.hpp>

#include <stk_io/IossBridge.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_util/parallel/Parallel.hpp>

namespace krino
{

inline const BoundingBox skinQueryBox(stk::math::Vector3d(-3, -3, -3), stk::math::Vector3d(3, 3, 3));

// A mesh of elements, each in a block of its own, as given by their node ids and the node coordinates,
// with the elements distributed round robin over the first processors, and the sides of the blocks'
// skins: their exposed sides and the sides between them.
struct MultiBlockMesh
{
  MultiBlockMesh(const stk::topology elemTopology,
      const std::vector<stk::math::Vector3d> & nodeLocs,
      const std::vector<stk::mesh::EntityIdVector> & elemNodeIds,
      const stk::ParallelMachine comm = MPI_COMM_WORLD)
  {
    myMeta = stk::mesh::MeshBuilder().set_spatial_dimension(3).create_meta_data();
    for (size_t e = 0; e < elemNodeIds.size(); ++e)
    {
      myBlocks.push_back(&myMeta->declare_part_with_topology("block_" + std::to_string(e + 1), elemTopology));
      stk::io::put_io_part_attribute(*myBlocks.back());
      myMeta->set_part_id(*myBlocks.back(), e + 1); // so stk counts it as an element block, as when read
    }
    myCoords = &myMeta->declare_field<double>(stk::topology::NODE_RANK, "coordinates", 1);
    stk::mesh::put_field_on_mesh(*myCoords, myMeta->universal_part(), 3, nullptr);
    myMeta->set_coordinate_field(myCoords);
    myMeta->commit();

    myMesh = stk::mesh::MeshBuilder(comm).create(myMeta);
    const int me = myMesh->parallel_rank();
    const int numProcs = myMesh->parallel_size();
    auto owner = [&](const size_t e) { return int(e % numProcs); };

    myMesh->modification_begin();
    for (size_t e = 0; e < elemNodeIds.size(); ++e)
      if (owner(e) == me)
        stk::mesh::declare_element(*myMesh, *myBlocks[e], e + 1, elemNodeIds[e]);
    for (size_t e = 0; e < elemNodeIds.size(); ++e)
      for (size_t other = 0; other < elemNodeIds.size(); ++other)
        if (owner(e) == me && owner(other) != me)
          for (auto id : elemNodeIds[e])
            if (std::find(elemNodeIds[other].begin(), elemNodeIds[other].end(), id) != elemNodeIds[other].end())
              myMesh->add_node_sharing(myMesh->get_entity(stk::topology::NODE_RANK, id), owner(other));
    auto coordsData = myCoords->data<stk::mesh::ReadWrite>();
    for (stk::mesh::EntityId id = 1; id <= nodeLocs.size(); ++id)
    {
      const stk::mesh::Entity node = myMesh->get_entity(stk::topology::NODE_RANK, id);
      if (!myMesh->is_valid(node)) continue;
      auto x = coordsData.entity_values(node);
      for (unsigned d = 0; d < 3; ++d)
        x(stk::mesh::ComponentIdx(d)) = nodeLocs[id - 1][d];
    }
    myMesh->modification_end();

    const stk::mesh::Selector allBlocks = stk::mesh::selectUnion(myBlocks);
    stk::mesh::create_exposed_block_boundary_sides(*myMesh, allBlocks, {});
    stk::mesh::create_interior_block_boundary_sides(*myMesh, allBlocks, {});

    for (auto & ids : elemNodeIds)
    {
      stk::math::Vector3d centroid = stk::math::Vector3d::ZERO;
      for (auto id : ids)
        centroid += nodeLocs[id - 1];
      myElemCentroids.push_back(centroid / ids.size());
      myElemNodeIds.emplace_back(ids.begin(), ids.end());
    }
  }

  std::unique_ptr<MeshSurface> build_skin(const std::vector<unsigned> & blockIndices, const int sign) const
  {
    stk::mesh::PartVector blocks;
    for (auto b : blockIndices)
      blocks.push_back(myBlocks[b]);
    auto skin = build_skin_mesh_surface(*myMeta, *myCoords, blocks, sign);
    skin->prepare_to_compute(0., skinQueryBox, 0.);
    return skin;
  }

  std::shared_ptr<stk::mesh::MetaData> myMeta;
  std::unique_ptr<stk::mesh::BulkData> myMesh;
  stk::mesh::PartVector myBlocks;
  stk::mesh::Field<double> * myCoords{nullptr};
  std::vector<stk::math::Vector3d> myElemCentroids;
  std::vector<std::set<stk::mesh::EntityId>> myElemNodeIds;
};

// Two tets, in block_1 and block_2, sharing the face with nodes 2, 3 and 4.
inline MultiBlockMesh two_tet_blocks(const stk::ParallelMachine comm = MPI_COMM_WORLD)
{
  return MultiBlockMesh(stk::topology::TETRAHEDRON_4,
      {{0, 0, 0}, {1, 0, 0}, {0, 1, 0}, {0, 0, 1}, {1, 1, 1}},
      {{1, 2, 3, 4}, {2, 3, 4, 5}},
      comm);
}

// Two unit hexes, in block_1 and block_2, sharing the quadrilateral x = 1 (nodes 2, 3, 7, 6), which is
// not planar if sharedNodeShift moves node 7 along x.
inline MultiBlockMesh two_hex_blocks(const double sharedNodeShift, const stk::ParallelMachine comm = MPI_COMM_WORLD)
{
  std::vector<stk::math::Vector3d> nodeLocs;
  for (double z : {0., 1.})
    for (auto [x, y] : std::array<std::array<double, 2>, 6>{{{0, 0}, {1, 0}, {1, 1}, {0, 1}, {2, 0}, {2, 1}}})
      nodeLocs.push_back(stk::math::Vector3d(x, y, z));
  // ids 1-6 at z=0 and 7-12 at z=1 in the order above; the shared face is x=1: ids 2, 3, 9, 8.
  nodeLocs[8][0] += sharedNodeShift;
  return MultiBlockMesh(stk::topology::HEXAHEDRON_8, nodeLocs, {{1, 2, 3, 4, 7, 8, 9, 10}, {2, 5, 6, 3, 8, 11, 12, 9}}, comm);
}

// Writes the mesh made by buildMesh (given MPI_COMM_SELF) to an Exodus file, on processor 0 alone, so that
// every processor can read it.
template <class BUILD_MESH>
void write_mesh_file_on_processor_0(const std::string & filename, const BUILD_MESH & buildMesh)
{
  if (0 == stk::parallel_machine_rank(MPI_COMM_WORLD))
  {
    const MultiBlockMesh mesh = buildMesh(MPI_COMM_SELF);
    stk::io::StkMeshIoBroker io(MPI_COMM_SELF);
    io.set_bulk_data(*mesh.myMesh);
    const size_t output = io.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    io.write_output_mesh(output);
  }
  stk::parallel_machine_barrier(MPI_COMM_WORLD);
}

} // namespace krino

#endif /* KRINO_KRINO_UNIT_TEST_UTILS_AKRI_UNIT_MULTIBLOCKMESH_HPP_ */
