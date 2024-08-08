/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef MOCK_MESH_UTILS_HPP
#define MOCK_MESH_UTILS_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include "StkMesh.hpp"
#include "SparcMesh.hpp"
#include <memory>
#include <string>

namespace mock_utils {

inline
void read_mesh(MPI_Comm comm,
               const std::string& fileName,
               const std::string& partName,
               const std::vector<std::string>& fieldNames,
               std::shared_ptr<mock::StkMesh>& mesh)
{
  std::shared_ptr<stk::mesh::BulkData> bulk = stk::mesh::MeshBuilder(comm).create();
  stk::mesh::MetaData& meta = bulk->mesh_meta_data();

  stk::io::StkMeshIoBroker ioBroker(comm);
  ioBroker.set_bulk_data(bulk);
  ioBroker.add_mesh_database(fileName, stk::io::READ_MESH);
  ioBroker.create_input_mesh();

  stk::mesh::Part* surfacePart = meta.get_part(partName);
  STK_ThrowRequireMsg(surfacePart != nullptr, std::string("Error, didn't find part: ") + partName);
  constexpr double initialValue = 0.0;
  for(const std::string& fieldName : fieldNames) {
    stk::mesh::Field<double>& field = meta.declare_field<double>(stk::topology::FACE_RANK, fieldName);
    stk::mesh::put_field_on_mesh(field, *surfacePart, &initialValue);
  }

  ioBroker.populate_bulk_data();

  mesh.reset(new mock::StkMesh(bulk, surfacePart->name()));
}

inline
void read_mesh(MPI_Comm comm,
               const std::string& fileName,
               const std::string& partName,
               const std::vector<std::string>& fieldNames,
               std::shared_ptr<mock::SparcMesh>& mesh)
{
  std::shared_ptr<mock::StkMesh> stkMesh;
  read_mesh(comm, fileName, partName, fieldNames, stkMesh);

  std::shared_ptr<stk::mesh::BulkData> stkBulk = stkMesh->get_stk_mesh();
  stk::mesh::Part* surfacePart = stkBulk->mesh_meta_data().get_part(partName);
  STK_ThrowRequireMsg(surfacePart != nullptr, std::string("Error, didn't find assumed part: ") + partName);

  stk::mesh::Selector ownedSurface = stkBulk->mesh_meta_data().locally_owned_part() & *surfacePart;
  stk::mesh::EntityVector sides;
  stk::mesh::get_entities(*stkBulk, stk::topology::FACE_RANK, ownedSurface, sides);

  const stk::mesh::FieldBase* coordField = stkBulk->mesh_meta_data().coordinate_field();
  std::vector<mock::SparcSide> sparcSides(sides.size());

  unsigned i = 0;
  for(stk::mesh::Entity side : sides) {
    const stk::mesh::Entity* nodes = stkBulk->begin_nodes(side);
    const unsigned numNodes = stkBulk->num_nodes(side);
    mock::SparcSide& sparcSide = sparcSides[i++];
    sparcSide.key = stkBulk->identifier(side);
    sparcSide.numNodes = numNodes;
    for(unsigned n=0; n<numNodes; ++n) {
      const double* coordData = reinterpret_cast<const double*>(stk::mesh::field_data(*coordField, nodes[n]));
      for(unsigned d=0; d<3; ++d) {
        sparcSide.nodeCoords[n][d] = coordData[d];
      }
    }
  }

  mesh.reset(new mock::SparcMesh(comm, sparcSides));
}

}

#endif // MOCK_MESH_UTILS_HPP

