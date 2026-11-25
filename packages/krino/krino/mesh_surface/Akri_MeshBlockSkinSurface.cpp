// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include <Akri_MeshBlockSkinSurface.hpp>

#include <Akri_DiagWriter.hpp>
#include <Akri_MeshSurface.hpp>

#include <Ioss_Property.h>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/ExodusTranslator.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_util/util/ReportHandler.hpp>

#include <map>
#include <sstream>

namespace krino {

BlockSkinMeshFile::BlockSkinMeshFile(const stk::ParallelMachine comm, const std::string & filename)
    : myComm(comm), myFilename(filename)
{
}

BlockSkinMeshFile::~BlockSkinMeshFile() = default;

void BlockSkinMeshFile::read_mesh_with_skin_sides()
{
  krinolog << "Reading the mesh of " << myFilename << " for the skins of its blocks." << stk::diag::dendl;
  myIoBroker = std::make_unique<stk::io::StkMeshIoBroker>(myComm);
  myIoBroker->property_add(Ioss::Property("DECOMPOSITION_METHOD", "RIB"));
  myIoBroker->add_mesh_database(myFilename, stk::io::READ_MESH);
  myIoBroker->create_input_mesh();
  std::shared_ptr<stk::mesh::BulkData> mesh = stk::mesh::MeshBuilder(myComm)
                                                  .set_aura_option(stk::mesh::BulkData::NO_AUTO_AURA)
                                                  .create(myIoBroker->meta_data_ptr());
  myIoBroker->set_bulk_data(mesh);
  myIoBroker->populate_bulk_data();

  const stk::mesh::Part & universal = mesh->mesh_meta_data().universal_part();
  stk::mesh::create_exposed_block_boundary_sides(*mesh, universal, {});
  stk::mesh::create_interior_block_boundary_sides(*mesh, universal, {});
}

const stk::mesh::BulkData & BlockSkinMeshFile::mesh()
{
  if (!myIoBroker)
    read_mesh_with_skin_sides();
  return myIoBroker->bulk_data();
}

static std::string names_of_element_blocks(const stk::mesh::MetaData & meta)
{
  std::ostringstream names;
  for (const auto * part : meta.get_parts())
    if (stk::mesh::is_element_block(*part))
      names << " " << part->name();
  return names.str();
}

stk::mesh::PartVector BlockSkinMeshFile::get_blocks(const std::vector<std::string> & blockNames)
{
  const stk::mesh::MetaData & meta = mesh().mesh_meta_data();
  stk::mesh::PartVector blocks;
  for (const auto & name : blockNames)
  {
    stk::mesh::Part * block = meta.get_part(name);
    STK_ThrowRequireMsg(nullptr != block && stk::mesh::is_element_block(*block),
        "The mesh file " << myFilename << " has no block named " << name << "; its blocks are:"
                         << names_of_element_blocks(meta) << ".");
    blocks.push_back(block);
  }
  return blocks;
}

std::shared_ptr<BlockSkinMeshFile> get_block_skin_mesh_file(const stk::ParallelMachine comm, const std::string & filename)
{
  static std::map<std::string, std::weak_ptr<BlockSkinMeshFile>> filesByName;
  std::weak_ptr<BlockSkinMeshFile> & sharedFile = filesByName[filename];
  if (auto file = sharedFile.lock(); file && file->parallel_comm() == comm)
    return file;
  auto file = std::make_shared<BlockSkinMeshFile>(comm, filename);
  sharedFile = file;
  return file;
}

static const stk::mesh::Field<double> & coordinates_field_of(const stk::mesh::BulkData & mesh)
{
  const stk::mesh::FieldBase * coords = mesh.mesh_meta_data().coordinate_field();
  STK_ThrowRequireMsg(nullptr != coords, "Mesh has no coordinates field.");
  return *static_cast<const stk::mesh::Field<double> *>(coords);
}

static BoundingBox compute_bounding_box_of_block_nodes(const stk::mesh::BulkData & mesh, const stk::mesh::PartVector & blocks)
{
  const unsigned dim = mesh.mesh_meta_data().spatial_dimension();
  const auto coordsData = coordinates_field_of(mesh).data();
  BoundingBox bbox;
  for (auto * bucket : mesh.get_buckets(stk::topology::NODE_RANK, stk::mesh::selectUnion(blocks)))
    for (auto node : *bucket)
      bbox.accommodate(stk::math::Vector3d(coordsData.entity_values(node), dim));
  bbox.global_reduce();
  return bbox;
}

MeshBlockSkinSurface::MeshBlockSkinSurface(const unsigned dim,
    std::shared_ptr<BlockSkinMeshFile> file,
    const std::vector<std::string> & blockNames,
    const int sign)
    : MeshFacetedSurface(file->parallel_comm(), dim, true),
      myDim(dim),
      myFile(file),
      myFilename(file->filename()),
      myBlockNames(blockNames),
      mySign(sign)
{
}

const BoundingBox & MeshBlockSkinSurface::blocks_bounding_box() const
{
  if (!myBlocksBoundingBox)
  {
    const stk::mesh::BulkData & mesh = myFile->mesh();
    myBlocksBoundingBox = compute_bounding_box_of_block_nodes(mesh, myFile->get_blocks(myBlockNames));
  }
  return *myBlocksBoundingBox;
}

void MeshBlockSkinSurface::insert_into(BoundingBox & bbox) const
{
  bbox.accommodate(blocks_bounding_box());
}

void MeshBlockSkinSurface::build_local_facets(const std::vector<BoundingBox> & /*procQueryBboxes*/,
    VertexIdsAndCoords & vertexIdsAndCoords,
    FacetVertexIds & facetVertexIds)
{
  STK_ThrowRequireMsg(myFile, "The skin of the blocks of " << myFilename << " was already built.");
  blocks_bounding_box(); // before the file is let go, for later insert_into
  const stk::mesh::BulkData & mesh = myFile->mesh();
  STK_ThrowRequireMsg(mesh.mesh_meta_data().spatial_dimension() == myDim,
      "The mesh file " << myFilename << " is " << mesh.mesh_meta_data().spatial_dimension()
                       << "D, but the skins of its blocks bound a " << myDim << "D mesh.");
  const stk::mesh::PartVector blocks = myFile->get_blocks(myBlockNames);
  krinolog << "Building facets for the skin of " << stk::mesh::selectUnion(blocks) << " of " << myFilename << stk::diag::dendl;
  append_skin_facets(mesh, coordinates_field_of(mesh), blocks, mesh.mesh_meta_data().locally_owned_part(), mySign, vertexIdsAndCoords, facetVertexIds);
  myFile.reset();
}

} // namespace krino
