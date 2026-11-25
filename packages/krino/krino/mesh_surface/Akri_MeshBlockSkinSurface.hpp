// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef KRINO_KRINO_MESH_SURFACE_AKRI_MESHBLOCKSKINSURFACE_HPP_
#define KRINO_KRINO_MESH_SURFACE_AKRI_MESHBLOCKSKINSURFACE_HPP_

#include <Akri_MeshFacetedSurface.hpp>

#include <stk_mesh/base/Types.hpp>
#include <stk_util/parallel/Parallel.hpp>

#include <memory>
#include <optional>
#include <string>
#include <vector>

namespace stk { namespace io { class StkMeshIoBroker; } }
namespace stk { namespace mesh { class BulkData; } }

namespace krino {

// The mesh of an Exodus file whose blocks' skins are surfaces.  It is read when first needed
// (collectively, decomposed by RIB, without aura), with the sides of its blocks' skins: the exposed sides
// and the sides between blocks.  The surfaces of the file share it (see get_block_skin_mesh_file) and let
// it go once they have their facets.
class BlockSkinMeshFile
{
public:
  BlockSkinMeshFile(const stk::ParallelMachine comm, const std::string & filename);
  ~BlockSkinMeshFile();

  const std::string & filename() const { return myFilename; }
  stk::ParallelMachine parallel_comm() const { return myComm; }

  // Collective at the first call, which reads the file.
  const stk::mesh::BulkData & mesh();

  // Collective at the first call.  The element blocks with these names; throws, listing the file's
  // blocks, if any name is not one of them.
  stk::mesh::PartVector get_blocks(const std::vector<std::string> & blockNames);

private:
  void read_mesh_with_skin_sides();

  stk::ParallelMachine myComm;
  std::string myFilename;
  std::unique_ptr<stk::io::StkMeshIoBroker> myIoBroker;
};

// The file of this name, shared with every surface of it that still holds it.
std::shared_ptr<BlockSkinMeshFile> get_block_skin_mesh_file(const stk::ParallelMachine comm, const std::string & filename);

// The skin of blocks of an Exodus file (see build_skin_mesh_surface) as a connected MeshFacetedSurface,
// built from the file's mesh at the first prepare_to_compute, after which the surface lets the mesh go.
// The facets face out of the blocks (sign = -1 reverses them), so the face between two blocks faces out
// of each in its skin, with the same facets.
class MeshBlockSkinSurface : public MeshFacetedSurface
{
public:
  MeshBlockSkinSurface(const unsigned dim,
      std::shared_ptr<BlockSkinMeshFile> file,
      const std::vector<std::string> & blockNames,
      const int sign);

  // Collective.  Accommodates the bounding box of the blocks' nodes (reading the file if it has not been
  // read yet).
  virtual void insert_into(BoundingBox & bbox) const override;

protected:
  virtual void build_local_facets(const std::vector<BoundingBox> & procQueryBboxes,
      VertexIdsAndCoords & vertexIdsAndCoords,
      FacetVertexIds & facetVertexIds) override;

private:
  const BoundingBox & blocks_bounding_box() const;

  unsigned myDim;
  std::shared_ptr<BlockSkinMeshFile> myFile;
  std::string myFilename;
  std::vector<std::string> myBlockNames;
  int mySign;
  mutable std::optional<BoundingBox> myBlocksBoundingBox;
};

} // namespace krino

#endif /* KRINO_KRINO_MESH_SURFACE_AKRI_MESHBLOCKSKINSURFACE_HPP_ */
