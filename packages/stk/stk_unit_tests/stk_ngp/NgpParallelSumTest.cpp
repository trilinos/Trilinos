#include <stk_ngp_test/ngp_test.hpp>
#include <stk_ngp/Ngp.hpp>
#include <stk_ngp/NgpFieldManager.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/Entity.hpp>
#include "stk_mesh/base/FieldParallel.hpp"
#include "stk_ngp/NgpFieldParallel.hpp"

namespace  {

class NgpParallelSum : public stk::unit_test_util::MeshFixture
{
protected:
  NgpParallelSum() : stk::unit_test_util::MeshFixture(3)
  {
  }

  void initialize_shared_values(stk::mesh::FieldBase & field)
  {
    const stk::mesh::BucketVector & buckets = get_bulk().get_buckets(stk::topology::NODE_RANK, get_meta().globally_shared_part());
    for (stk::mesh::Bucket * bucket : buckets) {
      for (const stk::mesh::Entity & node : *bucket) {
        std::vector<int> sharingProcs;
        stk::mesh::EntityKey nodeKey = get_bulk().entity_key(node);
        double * block2Values = static_cast<double*>(stk::mesh::field_data(field, node));
        get_bulk().comm_procs(nodeKey, sharingProcs);
        const size_t numSharers = sharingProcs.size() + 1;
        *block2Values /= numSharers;
      }
    }
  }

};

template <typename T>
void check_field_on_device(stk::mesh::BulkData & bulk,
                           ngp::Mesh &mesh,
                           ngp::FieldManager &fieldManager,
                           unsigned fieldOrdinal,
                           T expectedFieldValue)
{
    ngp::Field<T> & field = fieldManager.get_field<T>(fieldOrdinal);

    ngp::for_each_entity_run(mesh, stk::topology::NODE_RANK, bulk.mesh_meta_data().universal_part(), KOKKOS_LAMBDA(ngp::Mesh::MeshIndex entity)
    {
        NGP_EXPECT_NEAR(field.get(entity, 0), expectedFieldValue, 1.e-12);
    });
}

NGP_TEST_F(NgpParallelSum, simpleVersion)
{
  const double initValue = 1.0;
  const int numStates = 1;
  stk::mesh::Field<double> & field1 = get_meta().declare_field<stk::mesh::Field<double>>(stk::topology::NODE_RANK, "field1", numStates);
  stk::mesh::put_field_on_mesh(field1, get_meta().universal_part(), &initValue);

  setup_mesh("generated:1x1x4", stk::mesh::BulkData::NO_AUTO_AURA);

  initialize_shared_values(field1);

  ngp::Mesh ngpMesh(get_bulk());
  ngp::FieldManager fieldManager(get_bulk());
  ngp::Field<double> & deviceField1 = fieldManager.get_field<double>(field1.mesh_meta_data_ordinal());

  ngp::parallel_sum<double>(get_bulk(), std::vector<ngp::Field<double>*>{&deviceField1});

  check_field_on_device(get_bulk(), ngpMesh, fieldManager, field1.mesh_meta_data_ordinal(), 1.0);
}

}
