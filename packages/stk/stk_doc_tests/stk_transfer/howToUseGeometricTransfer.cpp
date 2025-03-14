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

#include <gtest/gtest.h>
#include <stk_util/parallel/Parallel.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_transfer/GeometricTransfer.hpp>

namespace {

//BEGIN_types
struct FieldConfigData {
  std::string name;
  stk::mesh::EntityRank rank;
  std::vector<double> initialValues;
};

using FieldConfig = std::vector<FieldConfigData>;
using BulkDataPtr = std::shared_ptr<stk::mesh::BulkData>;
//END_types

//BEGIN_send_adapter
class StkSendAdapter
{
public:
  using EntityKey = stk::mesh::EntityKey;
  using EntityProc = stk::search::IdentProc<EntityKey, int>;
  using EntityProcVec = std::vector<EntityProc>;
  using BoundingBox = std::pair<stk::search::Box<double>, EntityProc>;

  using Coords = std::array<double, 3>;

  StkSendAdapter(MPI_Comm globalComm, BulkDataPtr & bulk,
                 const std::string & partName, const FieldConfig & fieldConfig)
    : m_globalComm(globalComm),
      m_bulk(bulk),
      m_meta(bulk->mesh_meta_data()),
      m_part(m_meta.get_part(partName)),
      m_ghosting(nullptr)
  {
    for (const FieldConfigData & fieldConf : fieldConfig) {
      m_fields.push_back(m_meta.get_field<double>(fieldConf.rank, fieldConf.name));
    }
  }

  MPI_Comm comm() const { return m_globalComm; }

  void bounding_boxes(std::vector<BoundingBox> & searchDomain) const
  {
    stk::mesh::Selector ownedSelector = m_meta.locally_owned_part() & *m_part;
    const auto elements = stk::mesh::get_entities(*m_bulk, stk::topology::ELEM_RANK,
                                                  ownedSelector);
    searchDomain.clear();
    const int procInSearchComm = stk::parallel_machine_rank(m_globalComm);
    for (stk::mesh::Entity element : elements) {
      EntityProc entityProc(m_bulk->entity_key(element), procInSearchComm);
      searchDomain.emplace_back(get_box(element), entityProc);
    }
  }

  void copy_entities(const EntityProcVec & entitiesToSend, const std::string & /*name*/)
  {
    m_ghostedEntities.clear();

    for (auto keyProc : entitiesToSend) {
      const stk::mesh::EntityKey key = keyProc.id();
      const unsigned proc = keyProc.proc();
      m_ghostedEntities.emplace_back(m_bulk->get_entity(key), proc);
    }

    unsigned hasEntitiesToGhost = not m_ghostedEntities.empty();

    stk::all_reduce(m_globalComm, stk::ReduceSum<1>(&hasEntitiesToGhost));

    if (hasEntitiesToGhost) {
      stk::util::sort_and_unique(m_ghostedEntities);

      m_bulk->modification_begin();
      if (m_ghosting == nullptr) {
        m_ghosting = &m_bulk->create_ghosting("transfer_ghosting");
      }
      m_bulk->change_ghosting(*m_ghosting, m_ghostedEntities);
      m_bulk->modification_end();
    }
  }

  void update_values()
  {
    std::vector<const stk::mesh::FieldBase*> commFields;
    for (stk::mesh::Field<double> * field : m_fields) {
      commFields.push_back(static_cast<stk::mesh::FieldBase*>(field));
    }
    stk::mesh::communicate_field_data(*m_bulk, commFields);
  }

  Coords parametric_coords(EntityKey /*entityKey*/, const double * /*spatialCoordinates*/,
                           double & distance) const
  {
    distance = 0.0;
    return Coords{0.0, 0.0, 0.0};
  }

  void interpolate_fields(const Coords & /*parametricCoords*/, EntityKey entityKey,
                          unsigned numFields, const std::vector<unsigned> & fieldSizes,
                          const std::vector<double *> & recvFieldPtrs) const
  {
    // This is where the actual application-specific shape function interpolation
    // operation would go.  For simplicity, this example uses zeroth-order
    // interpolation from only the first node's value.
    const stk::mesh::Entity targetElement = m_bulk->get_entity(entityKey);
    const stk::mesh::Entity firstNode = m_bulk->begin_nodes(targetElement)[0];
    for (unsigned n = 0; n < numFields; ++n) {
      const double * fieldData = stk::mesh::field_data(*m_fields[n], firstNode);
      for (unsigned idx = 0; idx < fieldSizes[n]; ++idx) {
        recvFieldPtrs[n][idx] = fieldData[idx];
      }
    }
  }

private:
  stk::search::Box<double> get_box(stk::mesh::Entity element) const
  {
    constexpr double minDouble = std::numeric_limits<double>::lowest();
    constexpr double maxDouble = std::numeric_limits<double>::max();
    double minXYZ[3] = {maxDouble, maxDouble, maxDouble};
    double maxXYZ[3] = {minDouble, minDouble, minDouble};
    const auto * coordField =
        static_cast<const stk::mesh::Field<double>*>(m_meta.coordinate_field());

    const stk::mesh::Entity * nodes = m_bulk->begin_nodes(element);
    const unsigned numNodes = m_bulk->num_nodes(element);
    for (unsigned i = 0; i < numNodes; ++i) {
      const double * coords = stk::mesh::field_data(*coordField, nodes[i]);
      minXYZ[0] = std::min(minXYZ[0], coords[0]);
      minXYZ[1] = std::min(minXYZ[1], coords[1]);
      minXYZ[2] = std::min(minXYZ[2], coords[2]);
      maxXYZ[0] = std::max(maxXYZ[0], coords[0]);
      maxXYZ[1] = std::max(maxXYZ[1], coords[1]);
      maxXYZ[2] = std::max(maxXYZ[2], coords[2]);
    }

    constexpr double tol = 1.e-5;
    return stk::search::Box<double>(minXYZ[0]-tol, minXYZ[1]-tol, minXYZ[2]-tol,
                                    maxXYZ[0]+tol, maxXYZ[1]+tol, maxXYZ[2]+tol);
  }

  MPI_Comm m_globalComm;
  BulkDataPtr m_bulk;
  stk::mesh::MetaData & m_meta;
  stk::mesh::Part* m_part;
  std::vector<stk::mesh::Field<double>*> m_fields;
  stk::mesh::Ghosting * m_ghosting;
  stk::mesh::EntityProcVec m_ghostedEntities;

};
//END_send_adapter

//BEGIN_recv_adapter
class StkRecvAdapter
{
public:
  using EntityKey = stk::mesh::EntityKey;
  using EntityProc = stk::search::IdentProc<EntityKey>;
  using EntityProcVec = std::vector<EntityProc>;
  using BoundingBox = std::pair<stk::search::Sphere<double>, EntityProc>;

  using Point = stk::search::Point<double>;
  using Coords = std::array<double, 3>;

  StkRecvAdapter(MPI_Comm globalComm, BulkDataPtr & bulk,
                 const std::string & partName, const FieldConfig & fieldConfig)
    : m_globalComm(globalComm),
      m_bulk(bulk),
      m_meta(m_bulk->mesh_meta_data()),
      m_part(m_meta.get_part(partName))
  {
    for (const FieldConfigData & fieldConf : fieldConfig) {
      m_fields.push_back(m_meta.get_field<double>(fieldConf.rank, fieldConf.name));
    }
  }

  MPI_Comm comm() const { return m_globalComm; }

  void bounding_boxes(std::vector<BoundingBox> & searchRange) const
  {
    stk::mesh::Selector ownedSelector = m_meta.locally_owned_part() & *m_part;
    const auto nodes = stk::mesh::get_entities(*m_bulk, stk::topology::NODE_RANK,
                                               ownedSelector);
    constexpr double radius = 1.e-6;
    searchRange.clear();
    const int procInSearchComm = stk::parallel_machine_rank(m_globalComm);
    for (const stk::mesh::Entity & node : nodes) {
      EntityProc entityProc(m_bulk->entity_key(node), procInSearchComm);
      searchRange.emplace_back(stk::search::Sphere<double>(get_location(node), radius),
                               entityProc);
    }
  }

  void update_values()
  {
    std::vector<const stk::mesh::FieldBase*> commFields;
    for (stk::mesh::Field<double> * field : m_fields) {
      commFields.push_back(static_cast<stk::mesh::FieldBase*>(field));
    }
    stk::mesh::communicate_field_data(*m_bulk, commFields);
  }

  const double * node_coords(EntityKey entityKey) const
  {
    const stk::mesh::Entity node = m_bulk->get_entity(entityKey);
    const auto & coordField =
        *static_cast<const stk::mesh::Field<double>*>(m_meta.coordinate_field());
    return stk::mesh::field_data(coordField, node);
  }

  void save_parametric_coords(EntityKey entityKey, const Coords & parametricCoords)
  {
    m_sendParametricCoords[entityKey] = parametricCoords;
  }

  unsigned num_fields() { return m_fields.size(); }

  const Coords & get_parametric_coords(EntityKey entityKey)
  {
    return m_sendParametricCoords.at(entityKey);
  }

  double * field_values(EntityKey entityKey, unsigned fieldIndex)
  {
    const stk::mesh::Entity node = m_bulk->get_entity(entityKey);
    return stk::mesh::field_data(*m_fields[fieldIndex], node);
  }

  unsigned field_size(EntityKey entityKey, unsigned fieldIndex)
  {
    const stk::mesh::Entity node = m_bulk->get_entity(entityKey);
    return stk::mesh::field_scalars_per_entity(*m_fields[fieldIndex], node);
  }

private:
  Point get_location(stk::mesh::Entity node) const
  {
    const auto & coordField =
        *static_cast<const stk::mesh::Field<double>*>(m_meta.coordinate_field());
    const double * coords = stk::mesh::field_data(coordField, node);

    return Point(coords[0], coords[1], coords[2]);
  }

  MPI_Comm m_globalComm;
  BulkDataPtr m_bulk;
  stk::mesh::MetaData & m_meta;
  stk::mesh::Part * m_part;
  std::vector<stk::mesh::Field<double>*> m_fields;
  std::map<EntityKey, Coords> m_sendParametricCoords;
};
//END_recv_adapter

//BEGIN_interpolate
template<typename SendAdapter, typename RecvAdapter>
class Interpolate
{
public:
  using MeshA = SendAdapter;
  using MeshB = RecvAdapter;
  using EntityKeyA = typename MeshA::EntityKey;
  using EntityKeyB = typename MeshB::EntityKey;
  using EntityProcA = typename MeshA::EntityProc;
  using EntityProcB = typename MeshB::EntityProc;
  using EntityKeyMap = std::multimap<EntityKeyB, EntityKeyA>;
  using EntityProcRelation = std::pair<EntityProcB, EntityProcA>;
  using EntityProcRelationVec = std::vector<EntityProcRelation>;

  using Coords = typename MeshA::Coords;

  static void filter_to_nearest(EntityKeyMap & localRangeToDomain,
                                MeshA & sendMesh, MeshB & recvMesh)
  {
    using iterator = typename EntityKeyMap::iterator;
    using const_iterator = typename EntityKeyMap::const_iterator;

    for (const_iterator key = localRangeToDomain.begin(); key != localRangeToDomain.end();)
    {
      const EntityKeyB recvEntityKey = key->first;
      double closestDistance = std::numeric_limits<double>::max();

      const double * recvCoords = recvMesh.node_coords(recvEntityKey);

      std::pair<iterator, iterator> sendEntities = localRangeToDomain.equal_range(recvEntityKey);
      iterator nearest = sendEntities.second;

      for (iterator ii = sendEntities.first; ii != sendEntities.second; ++ii) {
        const EntityKeyA sendEntity = ii->second;

        double distance = 0;
        const Coords parametricCoords = sendMesh.parametric_coords(sendEntity, recvCoords,
                                                                   distance);

        if (distance < closestDistance) {
          closestDistance = distance;
          recvMesh.save_parametric_coords(recvEntityKey, parametricCoords);
          nearest = ii;
        }
      }

      key = sendEntities.second;
      if (nearest != sendEntities.first) {
        localRangeToDomain.erase(sendEntities.first, nearest);
      }
      if (nearest != sendEntities.second) {
        localRangeToDomain.erase(++nearest, sendEntities.second);
      }
    }
  }

  static void apply(MeshB & recvMesh, MeshA & sendMesh, EntityKeyMap & localRangeToDomain)
  {
    const unsigned numFields = recvMesh.num_fields();

    std::vector<double *> fieldPtrs(numFields);
    std::vector<unsigned> fieldSizes(numFields);

    typename EntityKeyMap::const_iterator ii;
    for (ii = localRangeToDomain.begin(); ii != localRangeToDomain.end(); ++ii) {
      const EntityKeyB recvNode = ii->first;
      const EntityKeyA sendElem = ii->second;

      const Coords & sendParametricCoords = recvMesh.get_parametric_coords(recvNode);

      for (unsigned n = 0; n < numFields; ++n) {
        fieldPtrs[n] = recvMesh.field_values(recvNode, n);
        fieldSizes[n] = recvMesh.field_size(recvNode, n);
      }

      sendMesh.interpolate_fields(sendParametricCoords, sendElem, numFields,
                                  fieldSizes, fieldPtrs);
    }
  }

};
//END_interpolate

//BEGIN_supporting_functions
BulkDataPtr read_mesh(MPI_Comm comm,
                      const std::string & fileName,
                      const FieldConfig & fieldConfig)
{
  BulkDataPtr bulk = stk::mesh::MeshBuilder(comm).create();

  stk::io::StkMeshIoBroker ioBroker(comm);
  ioBroker.set_bulk_data(bulk);
  ioBroker.add_mesh_database(fileName, stk::io::READ_MESH);
  ioBroker.create_input_mesh();

  stk::mesh::MetaData& meta = bulk->mesh_meta_data();
  for (const FieldConfigData & fieldConf : fieldConfig) {
    auto & field = meta.declare_field<double>(fieldConf.rank, fieldConf.name);
    stk::mesh::put_field_on_mesh(field, meta.universal_part(), fieldConf.initialValues.size(),
                                 fieldConf.initialValues.data());
  }

  ioBroker.populate_bulk_data();

  return bulk;
}

bool all_field_values_equal(BulkDataPtr & bulk, const FieldConfig & fieldConfig)
{
  stk::mesh::MetaData & meta = bulk->mesh_meta_data();

  for (const FieldConfigData & fieldConf : fieldConfig) {
    const auto & field = *meta.get_field<double>(fieldConf.rank, fieldConf.name);
    stk::mesh::Selector fieldSelector(*meta.get_part("block_1"));
    const auto nodes = stk::mesh::get_entities(*bulk, fieldConf.rank, fieldSelector);
    for (stk::mesh::Entity node : nodes) {
      const double* fieldData = stk::mesh::field_data(field, node);
      for (unsigned i = 0; i < fieldConf.initialValues.size(); ++i) {
        if (std::abs(fieldData[i] - fieldConf.initialValues[i]) > 1.e-6) {
          return false;
        }
      }
    }
  }

  return true;
}
//END_supporting_functions

//BEGIN_main_application
template <typename INTERPOLATE>
using GeomTransfer = stk::transfer::GeometricTransfer<INTERPOLATE>;

using TransferType = GeomTransfer<Interpolate<StkSendAdapter, StkRecvAdapter>>;

std::shared_ptr<TransferType> setup_transfer(MPI_Comm globalComm,
                                             BulkDataPtr & sendBulk, BulkDataPtr & recvBulk,
                                             const FieldConfig & sendFieldConfig,
                                             const FieldConfig & recvFieldConfig)
{
  auto sendAdapter = std::make_shared<StkSendAdapter>(globalComm, sendBulk, "block_1",
                                                      sendFieldConfig);
  auto recvAdapter = std::make_shared<StkRecvAdapter>(globalComm, recvBulk, "block_1",
                                                      recvFieldConfig);

  auto transfer = std::make_shared<TransferType>(sendAdapter, recvAdapter,
                                                 "demoTransfer", globalComm);

  transfer->initialize();

  return transfer;
}

TEST(StkTransferHowTo, useGeometricTransfer)
{
  MPI_Comm commWorld = MPI_COMM_WORLD;

  FieldConfig sendFieldConfig {{"temperature", stk::topology::NODE_RANK, {300.0}},
                               {"velocity", stk::topology::NODE_RANK, {1.0, 2.0, 3.0}}};
  FieldConfig recvFieldConfig {{"temperature", stk::topology::NODE_RANK, {0.0}},
                               {"velocity", stk::topology::NODE_RANK, {0.0, 0.0, 0.0}}};

  BulkDataPtr sendBulk = read_mesh(commWorld, "generated:1x1x4", sendFieldConfig);
  BulkDataPtr recvBulk = read_mesh(commWorld, "generated:1x1x4", recvFieldConfig);

  auto transfer = setup_transfer(commWorld, sendBulk, recvBulk,
                                 sendFieldConfig, recvFieldConfig);

  transfer->apply();
  EXPECT_TRUE(all_field_values_equal(recvBulk, sendFieldConfig));
}
//END_main_application

}
