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
#include <stk_coupling/SplitComms.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_transfer/ReducedDependencyGeometricTransfer.hpp>

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

  StkSendAdapter(MPI_Comm globalComm, BulkDataPtr & bulk,
                 const std::string & partName, const FieldConfig & fieldConfig)
    : m_globalComm(globalComm),
      m_bulk(bulk),
      m_meta(bulk->mesh_meta_data()),
      m_part(m_meta.get_part(partName))
  {
    unsigned totalFieldSize = 0;
    for (const FieldConfigData & fieldConf : fieldConfig) {
      m_fields.push_back(m_meta.get_field<double>(fieldConf.rank, fieldConf.name));
      totalFieldSize += fieldConf.initialValues.size();
    }
    m_totalFieldSize = totalFieldSize;
  }

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

  void update_values()
  {
    std::vector<const stk::mesh::FieldBase*> commFields;
    for (stk::mesh::Field<double> * field : m_fields) {
      commFields.push_back(static_cast<stk::mesh::FieldBase*>(field));
    }
    stk::mesh::communicate_field_data(*m_bulk, commFields);
  }

  void interpolate_fields(const std::array<double, 3> & parametricCoords,
                          EntityKey entityKey, double * interpValues) const
  {
    // This is where the actual application-specific shape function interpolation
    // operation would go.  For simplicity, this example uses zeroth-order
    // interpolation from only the first node's value.
    const stk::mesh::Entity targetElement = m_bulk->get_entity(entityKey);
    const stk::mesh::Entity firstNode = m_bulk->begin_nodes(targetElement)[0];
    unsigned offset = 0;
    for (const stk::mesh::Field<double> * field : m_fields) {
      const double * fieldData = stk::mesh::field_data(*field, firstNode);
      for (unsigned idx = 0; idx < field->max_size(); ++idx) {
        interpValues[offset++] = fieldData[idx];
      }
    }
  }

  unsigned total_field_size() const { return m_totalFieldSize; }

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
  unsigned m_totalFieldSize;
};
//END_send_adapter

//BEGIN_remote_send_adapter
class RemoteSendAdapter
{
public:
  using EntityKey = uint64_t;
  using EntityProc = stk::search::IdentProc<EntityKey>;
  using EntityProcVec = std::vector<EntityProc>;
  using BoundingBox = std::pair<stk::search::Box<double>, EntityProc>;

  void bounding_boxes(std::vector<BoundingBox> & ) const {}
  void update_values() {}
};
//END_remote_send_adapter

//BEGIN_recv_adapter
class StkRecvAdapter
{
public:
  using EntityKey = stk::mesh::EntityKey;
  using EntityProc = stk::search::IdentProc<EntityKey>;
  using EntityProcVec = std::vector<EntityProc>;

  using BoundingBox = std::pair<stk::search::Sphere<double>, EntityProc>;
  using Point = stk::search::Point<double>;
  using ToPointsContainer = std::vector<Point>;
  using ToPointsDistanceContainer = std::vector<double>;

  StkRecvAdapter(MPI_Comm globalComm, BulkDataPtr & bulk,
                 const std::string & partName, const FieldConfig & fieldConfig)
    : m_globalComm(globalComm),
      m_bulk(bulk),
      m_meta(m_bulk->mesh_meta_data()),
      m_part(m_meta.get_part(partName))
  {
    unsigned totalFieldSize = 0;
    for (const FieldConfigData & fieldConf : fieldConfig) {
      m_fields.push_back(m_meta.get_field<double>(fieldConf.rank, fieldConf.name));
      totalFieldSize += fieldConf.initialValues.size();
    }
    m_totalFieldSize = totalFieldSize;
  }

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

  void get_to_points_coordinates(const EntityProcVec & toEntityKeys,
                                 ToPointsContainer & toPoints)
  {
    toPoints.clear();
    for (EntityProc entityProc : toEntityKeys) {
      toPoints.push_back(get_location(m_bulk->get_entity(entityProc.id())));
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

  void set_field_values(const EntityKey & entityKey, const double * recvInterpValues)
  {
    stk::mesh::Entity node = m_bulk->get_entity(entityKey);
    unsigned offset = 0;
    for (const stk::mesh::Field<double> * field : m_fields) {
      double * fieldData = stk::mesh::field_data(*field, node);
      for (unsigned idx = 0; idx < field->max_size(); ++idx) {
        fieldData[idx] = recvInterpValues[offset++];
      }
    }
  }

  unsigned total_field_size() const { return m_totalFieldSize; }

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
  unsigned m_totalFieldSize;
};
//END_recv_adapter

//BEGIN_remote_recv_adapter
class RemoteRecvAdapter
{
public:
  using EntityKey = uint64_t;
  using EntityProc = stk::search::IdentProc<EntityKey>;
  using EntityProcVec = std::vector<EntityProc>;
  using BoundingBox = std::pair<stk::search::Sphere<double>, EntityProc>;

  using Point = stk::search::Point<double>;
  using ToPointsContainer = std::vector<Point>;
  using ToPointsDistance = double;
  using ToPointsDistanceContainer = std::vector<ToPointsDistance>;

  void bounding_boxes(std::vector<BoundingBox> & ) const {}
  void get_to_points_coordinates(const EntityProcVec & , ToPointsContainer & ) {}
  void update_values() {}
};
//END_remote_recv_adapter

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
  using EntityProcRelation = std::pair<EntityProcB, EntityProcA>;
  using EntityProcRelationVec = std::vector<EntityProcRelation>;

  void obtain_parametric_coords(
      const typename MeshA::EntityProcVec & elemsToInterpolateFrom,
      const MeshA & sendAdapter,
      const typename MeshB::ToPointsContainer & pointsToInterpolateTo,
      typename MeshB::ToPointsDistanceContainer & distanceToInterpolationPoints)
  {
    for (unsigned i = 0; i < elemsToInterpolateFrom.size(); ++i) {
      m_parametricCoords.push_back({0, 0, 0});
      distanceToInterpolationPoints.push_back(0.0);
    }
  }

  void mask_parametric_coords(const std::vector<int> & filterMaskFrom, int fromCount)
  {
    for (unsigned i = 0; i < filterMaskFrom.size(); ++i) {
      if (filterMaskFrom[i]) {
        m_maskedParametricCoords.push_back(m_parametricCoords[i]);
      }
    }
  }

  void apply(MeshB * recvAdapter, MeshA * sendAdapter,
             const typename MeshB::EntityProcVec & toEntityKeysMasked,
             const typename MeshA::EntityProcVec & fromEntityKeysMasked,
             const stk::transfer::ReducedDependencyCommData & commData)
  {
    const unsigned totalFieldSize = sendAdapter->total_field_size();
    std::vector<double> sendInterpValues(fromEntityKeysMasked.size() * totalFieldSize);
    std::vector<double> recvInterpValues(toEntityKeysMasked.size() * totalFieldSize);

    interpolate_from_send_mesh(fromEntityKeysMasked, *sendAdapter, sendInterpValues);
    stk::transfer::do_communication(commData, sendInterpValues, recvInterpValues,
                                    totalFieldSize);
    write_to_recv_mesh(recvInterpValues, toEntityKeysMasked, *recvAdapter);
  }

private:
  void interpolate_from_send_mesh(const typename MeshA::EntityProcVec & fromEntityKeysMasked,
                                  const MeshA & sendAdapter,
                                  std::vector<double> & sendInterpValues)
  {
    unsigned offset = 0;
    for (unsigned i = 0; i < fromEntityKeysMasked.size(); ++i) {
      typename MeshA::EntityKey key = fromEntityKeysMasked[i].id();
      sendAdapter.interpolate_fields(m_maskedParametricCoords[i], key, &sendInterpValues[offset]);
      offset += sendAdapter.total_field_size();
    }
  }

  void write_to_recv_mesh(const std::vector<double> & recvInterpValues,
                          const typename MeshB::EntityProcVec & toEntityKeysMasked,
                          MeshB & recvAdapter)
  {
    unsigned offset = 0;
    for (unsigned i = 0; i < toEntityKeysMasked.size(); ++i) {
      typename MeshB::EntityKey key = toEntityKeysMasked[i].id();
      recvAdapter.set_field_values(key, &recvInterpValues[offset]);
      offset += recvAdapter.total_field_size();
    }
  }
  std::vector<std::array<double, 3>> m_parametricCoords;
  std::vector<std::array<double, 3>> m_maskedParametricCoords;
};
//END_interpolate

//BEGIN_send_interpolate
template<typename SendAdapter, typename RecvAdapter>
class SendInterpolate
{
public:
  using MeshA = SendAdapter;
  using MeshB = RecvAdapter;
  using EntityKeyA = typename MeshA::EntityKey;
  using EntityKeyB = typename MeshB::EntityKey;
  using EntityProcA = typename MeshA::EntityProc;
  using EntityProcB = typename MeshB::EntityProc;
  using EntityProcRelation = std::pair<EntityProcB, EntityProcA>;
  using EntityProcRelationVec = std::vector<EntityProcRelation>;

  void obtain_parametric_coords(
      const typename MeshA::EntityProcVec & elemsToInterpolateFrom,
      const MeshA & sendAdapter,
      const typename MeshB::ToPointsContainer & pointsToInterpolateTo,
      typename MeshB::ToPointsDistanceContainer & distanceToInterpolationPoints)
  {
    for (unsigned i = 0; i < elemsToInterpolateFrom.size(); ++i) {
      m_parametricCoords.push_back({0, 0, 0});
      distanceToInterpolationPoints.push_back(0.0);
    }
  }

  void mask_parametric_coords(const std::vector<int> & filterMaskFrom, int fromCount)
  {
    for (unsigned i = 0; i < filterMaskFrom.size(); ++i) {
      if (filterMaskFrom[i]) {
        m_maskedParametricCoords.push_back(m_parametricCoords[i]);
      }
    }
  }

  void apply(MeshB * /*recvAdapter*/, MeshA * sendAdapter,
             const typename MeshB::EntityProcVec & /*toEntityKeysMasked*/,
             const typename MeshA::EntityProcVec & fromEntityKeysMasked,
             const stk::transfer::ReducedDependencyCommData & commData)
  {
    const unsigned totalFieldSize = sendAdapter->total_field_size();
    std::vector<double> sendInterpValues(fromEntityKeysMasked.size() * totalFieldSize);
    interpolate_from_send_mesh(fromEntityKeysMasked, *sendAdapter, sendInterpValues);

    std::vector<double> recvInterpValues;  // Unused
    stk::transfer::do_communication(commData, sendInterpValues, recvInterpValues,
                                    totalFieldSize);
  }

private:
  void interpolate_from_send_mesh(const typename MeshA::EntityProcVec & fromEntityKeysMasked,
                                  const MeshA & sendAdapter,
                                  std::vector<double> & sendInterpValues)
  {
    unsigned offset = 0;
    for (unsigned i = 0; i < fromEntityKeysMasked.size(); ++i) {
      typename MeshA::EntityKey key = fromEntityKeysMasked[i].id();
      sendAdapter.interpolate_fields(m_maskedParametricCoords[i], key, &sendInterpValues[offset]);
      offset += sendAdapter.total_field_size();
    }
  }

  std::vector<std::array<double, 3>> m_parametricCoords;
  std::vector<std::array<double, 3>> m_maskedParametricCoords;
};
//END_send_interpolate

//BEGIN_recv_interpolate
template<typename SendAdapter, typename RecvAdapter>
class RecvInterpolate
{
public:
  using MeshA = SendAdapter;
  using MeshB = RecvAdapter;
  using EntityKeyA = typename MeshA::EntityKey;
  using EntityKeyB = typename MeshB::EntityKey;
  using EntityProcA = typename MeshA::EntityProc;
  using EntityProcB = typename MeshB::EntityProc;
  using EntityProcRelation = std::pair<EntityProcB, EntityProcA>;
  using EntityProcRelationVec = std::vector<EntityProcRelation>;

  void obtain_parametric_coords(const typename MeshA::EntityProcVec , MeshA & ,
                                const typename MeshB::ToPointsContainer & ,
                                const typename MeshB::ToPointsDistanceContainer & ) {}

  void mask_parametric_coords(const std::vector<int> & , int ) {}

  void apply(MeshB * recvAdapter, MeshA * /*sendAdapter*/,
             const typename MeshB::EntityProcVec & toEntityKeysMasked,
             const typename MeshA::EntityProcVec & /*fromEntityKeysMasked*/,
             const stk::transfer::ReducedDependencyCommData & comm_data)
  {
    const unsigned totalFieldSize = recvAdapter->total_field_size();
    std::vector<double> sendInterpValues;  // Unused
    std::vector<double> recvInterpValues(toEntityKeysMasked.size() * totalFieldSize);
    stk::transfer::do_communication(comm_data, sendInterpValues, recvInterpValues,
                                    totalFieldSize);

    write_to_recv_mesh(recvInterpValues, toEntityKeysMasked, *recvAdapter);
  }

private:
  void write_to_recv_mesh(const std::vector<double> & recvInterpValues,
                          const typename MeshB::EntityProcVec & toEntityKeysMasked,
                          MeshB & recvAdapter)
  {
    unsigned offset = 0;
    for (unsigned i = 0; i < toEntityKeysMasked.size(); ++i) {
      typename MeshB::EntityKey key = toEntityKeysMasked[i].id();
      recvAdapter.set_field_values(key, &recvInterpValues[offset]);
      offset += recvAdapter.total_field_size();
    }
  }

};
//END_recv_interpolate

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

namespace spmd {

//BEGIN_main_application_spmd
template <typename INTERPOLATE>
using RDGeomTransfer = stk::transfer::ReducedDependencyGeometricTransfer<INTERPOLATE>;

using TransferType = RDGeomTransfer<Interpolate<StkSendAdapter, StkRecvAdapter>>;

std::shared_ptr<TransferType> setup_transfer(MPI_Comm globalComm,
                                             BulkDataPtr & sendBulk,
                                             BulkDataPtr & recvBulk,
                                             const FieldConfig & sendFieldConfig,
                                             const FieldConfig & recvFieldConfig)
{
  auto sendAdapter = std::make_shared<StkSendAdapter>(globalComm, sendBulk, "block_1", sendFieldConfig);
  auto recvAdapter = std::make_shared<StkRecvAdapter>(globalComm, recvBulk, "block_1", recvFieldConfig);

  auto transfer = std::make_shared<TransferType>(sendAdapter, recvAdapter, "demoTransfer", globalComm);

  transfer->initialize();

  return transfer;
}

TEST(StkTransferHowTo, useReducedDependencyGeometricTransferSPMD)
{
  MPI_Comm commWorld = MPI_COMM_WORLD;

  FieldConfig sendFieldConfig {{"temperature", stk::topology::NODE_RANK, {300.0}},
                               {"velocity", stk::topology::NODE_RANK, {1.0, 2.0, 3.0}}};
  FieldConfig recvFieldConfig {{"temperature", stk::topology::NODE_RANK, {0.0}},
                               {"velocity", stk::topology::NODE_RANK, {0.0, 0.0, 0.0}}};

  BulkDataPtr sendBulk = read_mesh(commWorld, "generated:1x1x4", sendFieldConfig);
  BulkDataPtr recvBulk = read_mesh(commWorld, "generated:1x1x4", recvFieldConfig);

  auto transfer = setup_transfer(commWorld, sendBulk, recvBulk, sendFieldConfig, recvFieldConfig);

  transfer->apply();
  EXPECT_TRUE(all_field_values_equal(recvBulk, sendFieldConfig));
}
//END_main_application_spmd

}

namespace mpmd {

//BEGIN_main_application_mpmd
template <typename INTERPOLATE>
using RDGeomTransfer = stk::transfer::ReducedDependencyGeometricTransfer<INTERPOLATE>;

using SendTransferType = RDGeomTransfer<SendInterpolate<StkSendAdapter, RemoteRecvAdapter>>;
using RecvTransferType = RDGeomTransfer<RecvInterpolate<RemoteSendAdapter, StkRecvAdapter>>;

std::shared_ptr<SendTransferType> setup_send_transfer(MPI_Comm globalComm, BulkDataPtr & bulk,
                                                      const FieldConfig & fieldConfig)
{
  auto sendAdapter = std::make_shared<StkSendAdapter>(globalComm, bulk, "block_1",
                                                      fieldConfig);
  std::shared_ptr<RemoteRecvAdapter> nullRecvAdapter;

  auto sendTransfer = std::make_shared<SendTransferType>(sendAdapter, nullRecvAdapter,
                                                         "SendTransfer", globalComm);

  sendTransfer->initialize();

  return sendTransfer;
}

std::shared_ptr<RecvTransferType> setup_recv_transfer(MPI_Comm globalComm,
                                                      BulkDataPtr & mesh,
                                                      const FieldConfig & fieldConfig)
{
  std::shared_ptr<RemoteSendAdapter> nullSendAdapter;
  auto recvAdapter = std::make_shared<StkRecvAdapter>(globalComm, mesh, "block_1",
                                                      fieldConfig);

  auto recvTransfer = std::make_shared<RecvTransferType>(nullSendAdapter, recvAdapter,
                                                         "RecvTransfer", globalComm);
  recvTransfer->initialize();

  return recvTransfer;
}

TEST(StkTransferHowTo, useReducedDependencyGeometricTransferMPMD)
{
  MPI_Comm commWorld = MPI_COMM_WORLD;
  int numProcs = stk::parallel_machine_size(commWorld);
  int myRank = stk::parallel_machine_rank(commWorld);
  if (numProcs < 2) return;

  int color = myRank < numProcs/2 ? 0 : 1;

  stk::coupling::SplitComms splitComms(commWorld, color);
  splitComms.set_free_comms_in_destructor(true);
  MPI_Comm myComm = splitComms.get_split_comm();

  FieldConfig sendFieldConfig {{"temperature", stk::topology::NODE_RANK, {300.0}},
                               {"velocity", stk::topology::NODE_RANK, {1.0, 2.0, 3.0}}};
  FieldConfig recvFieldConfig {{"temperature", stk::topology::NODE_RANK, {0.0}},
                               {"velocity", stk::topology::NODE_RANK, {0.0, 0.0, 0.0}}};

  if (color == 0) {
    BulkDataPtr bulk = read_mesh(myComm, "generated:1x1x4", sendFieldConfig);
    auto sendTransfer = setup_send_transfer(commWorld, bulk, sendFieldConfig);
    sendTransfer->apply();
  }
  else if (color == 1) {
    BulkDataPtr bulk = read_mesh(myComm, "generated:1x1x4", recvFieldConfig);
    auto recvTransfer = setup_recv_transfer(commWorld, bulk, recvFieldConfig);
    recvTransfer->apply();
    EXPECT_TRUE(all_field_values_equal(bulk, sendFieldConfig));
  }
}
//END_main_application_mpmd

}
}
