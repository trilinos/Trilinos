#include <gtest/gtest.h>

#include <stk_mesh/base/BulkData.hpp>
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

  double* get_coordinates(unsigned node_index) const
  {
    stk::mesh::Entity node = m_localIdToNode[node_index];
    return stk::mesh::field_data(*m_coords, node);
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
        double* node_data = localIdBulkData.get_coordinates(connectivity[j]);
        EXPECT_NEAR(gold_x_coordinates[i][j], node_data[0], 1.e-6);
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
          double *node_data = stk::mesh::field_data(*coords, nodes[k]);
          EXPECT_NEAR(gold_x_coordinates[elemIndex][k], node_data[0], 1.e-6);
        }
        ++elemIndex;
      }
    }
  }
}
