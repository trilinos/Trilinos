#include <gtest/gtest.h>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/MeshBuilder.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldBase.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Types.hpp>

#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_util/util/PairIter.hpp>   // for PairIter

#include <stk_unit_test_utils/MeshFixture.hpp>

#include "stk_io/StkMeshIoBroker.hpp"
#include "stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp"

TEST(StkMesh, well_defined_local_ids_no_Aura)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> meshPtr = stk::mesh::MeshBuilder(comm)
                                                  .set_maintain_local_ids(true).create();

  std::string meshSpec = "generated:1x1x2";
  stk::io::fill_mesh(meshSpec, *meshPtr);

  const int myRank = stk::parallel_machine_rank(comm);
  using VecPairs = std::vector<std::pair<stk::mesh::EntityId,unsigned>>;
  VecPairs expected = (myRank==0) ?
  VecPairs{{1,0}, {2,1}, {3,2}, {4,3}, {5,4}, {6,5}, {7,6}, {8,7}}
  :
  VecPairs{{5,4}, {6,5}, {7,6}, {8,7}, {9,0}, {10,1}, {11,2}, {12,3}};

  for(const auto& nodeIdLocalId : expected) {
    stk::mesh::EntityId nodeId = nodeIdLocalId.first;
    unsigned expectedLocalId = nodeIdLocalId.second;
    stk::mesh::Entity node = meshPtr->get_entity(stk::topology::NODE_RANK, nodeId);
    ASSERT_TRUE(meshPtr->is_valid(node));
    EXPECT_EQ(expectedLocalId, meshPtr->local_id(node))<<"P"<<myRank<<" nodeId"<<nodeId;
  }
}

TEST(StkMesh, well_defined_local_ids_with_Aura)
{
  MPI_Comm comm = stk::parallel_machine_world();
  if (stk::parallel_machine_size(comm) != 2) { GTEST_SKIP(); }

  std::shared_ptr<stk::mesh::BulkData> meshPtr = stk::mesh::MeshBuilder(comm)
                                                  .set_maintain_local_ids(true).create();

  std::string meshSpec = "generated:1x1x2";
  stk::io::fill_mesh(meshSpec, *meshPtr);

  const int myRank = stk::parallel_machine_rank(comm);
  using VecPairs = std::vector<std::pair<stk::mesh::EntityId,unsigned>>;
  VecPairs expected = (myRank==0) ?
  VecPairs{{1,0}, {2,1}, {3,2}, {4,3}, {5,4}, {6,5}, {7,6}, {8,7}, {9,8}, {10,9}, {11,10}, {12,11}}
  :
  VecPairs{{1,8}, {2,9}, {3,10}, {4,11}, {5,4}, {6,5}, {7,6}, {8,7}, {9,0}, {10,1}, {11,2}, {12,3}};

  for(const auto& nodeIdLocalId : expected) {
    stk::mesh::EntityId nodeId = nodeIdLocalId.first;
    unsigned expectedLocalId = nodeIdLocalId.second;
    stk::mesh::Entity node = meshPtr->get_entity(stk::topology::NODE_RANK, nodeId);
    ASSERT_TRUE(meshPtr->is_valid(node));
    EXPECT_EQ(expectedLocalId, meshPtr->local_id(node))<<"P"<<myRank<<" nodeId"<<nodeId;
  }
}

class LocalIds : public stk::unit_test_util::MeshFixture
{
protected:
  LocalIds()
    : stk::unit_test_util::MeshFixture(3)
  {}
  virtual ~LocalIds() {}
};


class LocalIdBulkData
{
public:
  LocalIdBulkData(stk::mesh::BulkData& bulkData) : m_bulkData(bulkData), bulkDataLocalIdMapper()
  {
    bulkDataLocalIdMapper.set_size(bulkData);
    m_localIdToElement.resize(num_elements());
    fillIds(stk::topology::ELEM_RANK, m_localIdToElement);
    m_localIdToNode.resize(num_nodes());
    fillIds(stk::topology::NODE_RANK, m_localIdToNode);
    m_localIdToNode.resize(num_nodes());
    m_coords = m_bulkData.mesh_meta_data().get_field<double>(stk::topology::NODE_RANK, "coordinates");
  }

  ~LocalIdBulkData() {}

  unsigned num_elements() const
  {
    return num(stk::topology::ELEM_RANK);
  }

  unsigned num_nodes() const
  {
    return num(stk::topology::NODE_RANK);
  }

  std::vector<unsigned> get_nodes(unsigned elem_index) const
  {
    std::vector<unsigned> node_conn(m_bulkData.num_nodes(m_localIdToElement[elem_index]));
    const stk::mesh::Entity* nodes = m_bulkData.begin_nodes(m_localIdToElement[elem_index]);
    for(size_t i=0;i<node_conn.size();++i)
      node_conn[i] = get_local_id(nodes[i]);

    return node_conn;
  }

  stk::mesh::EntityValues<const double> get_coordinates(unsigned node_index) const
  {
    stk::mesh::Entity node = m_localIdToNode[node_index];
    return m_coords->data().entity_values(node);
  }

  void reset_local_ids()
  {
    bulkDataLocalIdMapper.set_size(m_bulkData);
    fillIds(stk::topology::ELEM_RANK, m_localIdToElement);
    fillIds(stk::topology::NODE_RANK, m_localIdToNode);
  }

private:

  unsigned get_local_id(stk::mesh::Entity entity) const
  {
    return bulkDataLocalIdMapper.entity_to_local(entity);
  }

  void fillIds(stk::mesh::EntityRank rank, std::vector<stk::mesh::Entity>& localIdToEntity)
  {
    size_t counter = 0;
    const stk::mesh::BucketVector& elemBuckets =m_bulkData.buckets(rank);
    for(size_t i=0;i<elemBuckets.size();++i)
    {
      const stk::mesh::Bucket& bucket = *elemBuckets[i];
      for(size_t j=0;j<bucket.size();++j)
      {
        localIdToEntity[counter] = bucket[j];
        bulkDataLocalIdMapper.add_new_entity_with_local_id(bucket[j], counter);
        counter++;
      }
    }
  }

  unsigned num(stk::mesh::EntityRank rank) const
  {
    return stk::mesh::count_selected_entities(m_bulkData.mesh_meta_data().universal_part(), m_bulkData.buckets(rank));
  }

  typedef stk::mesh::Field<double> CoordFieldType;

  stk::mesh::BulkData& m_bulkData;
  std::vector<stk::mesh::Entity> m_localIdToElement; // size num local elements
  std::vector<stk::mesh::Entity> m_localIdToNode;    // size num local nodes
  CoordFieldType *m_coords = nullptr;
  stk::mesh::impl::LocalIdMapper bulkDataLocalIdMapper;
};

std::vector<std::vector<double> > gold_x_coordinates =
{
  {0, 1, 1, 0, 0, 1, 1, 0 },
  {1, 2, 2, 1, 1, 2, 2, 1 },
  {0, 1, 1, 0, 0, 1, 1, 0 },
  {1, 2, 2, 1, 1, 2, 2, 1 },
  {0, 1, 1, 0, 0, 1, 1, 0 },
  {1, 2, 2, 1, 1, 2, 2, 1 },
  {0, 1, 1, 0, 0, 1, 1, 0 },
  {1, 2, 2, 1, 1, 2, 2, 1 }
};

TEST_F(LocalIds, using_local_ids)
{
  if (get_parallel_size() == 1)
  {
    setup_mesh("generated:2x2x2", stk::mesh::BulkData::AUTO_AURA);

    LocalIdBulkData localIdBulkData(get_bulk());

    for(size_t i=0;i<localIdBulkData.num_elements();++i)
    {
      std::vector<unsigned> connectivity = localIdBulkData.get_nodes(i);
      for(size_t j=0;j<connectivity.size();++j)
      {
        auto node_data = localIdBulkData.get_coordinates(connectivity[j]);
        EXPECT_NEAR(gold_x_coordinates[i][j], node_data(0_comp), 1.e-6);
      }
    }
  }
}

typedef stk::PairIter<const stk::mesh::Entity*> Entities;

class BulkDataHelper
{
public:
  BulkDataHelper(const stk::mesh::BulkData& bulkData) : m_bulkData(bulkData)
  {}

  Entities get_nodes(stk::mesh::Entity element) const
  {
    return Entities(m_bulkData.begin_nodes(element), m_bulkData.end_nodes(element));
  }
private:
  const stk::mesh::BulkData& m_bulkData;
};

TEST_F(LocalIds, using_entities)
{
  if (get_parallel_size() == 1)
  {
    setup_mesh("generated:2x2x2", stk::mesh::BulkData::AUTO_AURA);
    BulkDataHelper bulkDataHelper(get_bulk());

    typedef stk::mesh::Field<double> CoordFieldType;
    CoordFieldType *coords = get_meta().get_field<double>(stk::topology::NODE_RANK, "coordinates");
    auto coordData = coords->data();

    unsigned elemIndex = 0;
    const stk::mesh::BucketVector& elemBuckets = get_bulk().buckets(stk::topology::ELEM_RANK);
    for(size_t i=0;i<elemBuckets.size();++i)
    {
      const stk::mesh::Bucket& bucket = *elemBuckets[i];
      for(size_t j=0;j<bucket.size();++j)
      {
        Entities nodes = bulkDataHelper.get_nodes(bucket[j]);
        for(unsigned k=0;k<nodes.size();++k)
        {
          auto nodeData = coordData.entity_values(nodes[k]);
          EXPECT_NEAR(gold_x_coordinates[elemIndex][k], nodeData(0_comp), 1.e-6);
        }
        ++elemIndex;
      }
    }
  }
}
