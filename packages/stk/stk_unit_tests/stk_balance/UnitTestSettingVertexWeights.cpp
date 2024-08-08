#include <stk_unit_test_utils/unittestMeshUtils.hpp>
#include <stk_mesh/base/Field.hpp>
#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_mesh/baseImpl/elementGraph/BulkDataIdMapper.hpp>

namespace
{

////////////////////////////////////////////////////////////////////////////////////////////

class VertexWeightSettings : public stk::unit_test_util::MeshFixture
{
public:

  VertexWeightSettings() : MeshFixture(), vertexWeightField(nullptr)
  {
    setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
    vertexWeightField = &get_meta().declare_field<double>(stk::topology::ELEMENT_RANK, "vertex_weights", 1);
    stk::mesh::put_field_on_mesh(*vertexWeightField, get_meta().universal_part(), nullptr);
  }

  virtual ~VertexWeightSettings() {}

  void setup_mesh_with_global_element_count_set(const std::string &meshSpecification, stk::mesh::BulkData::AutomaticAuraOption auraOption)
  {
    stk::io::fill_mesh(meshSpecification, get_bulk());
    set_global_element_count();
    fill_vertex_weight_field_data();
  }

  void test_setting_of_vertex_weights_via_field()
  {
    stk::mesh::impl::LocalIdMapper localIds(get_bulk(), stk::topology::ELEM_RANK);

    Zoltan2ParallelGraph zoltan2Graph;
    stk::balance::FieldVertexWeightSettings graphSettings(get_bulk(), *vertexWeightField, 0.0);
    stk::mesh::Selector sel = get_meta().universal_part();
    stk::balance::internal::createZoltanParallelGraph(get_bulk(), sel, get_comm(), graphSettings, zoltan2Graph);

    zoltan2Graph.adjust_vertex_weights(graphSettings, get_bulk(), sel, localIds);
    check_graph_vertex_weights(zoltan2Graph.get_vertex_ids(), zoltan2Graph.get_vertex_weights());
  }

private:

  void set_global_element_count()
  {
    std::vector<size_t> counts;
    stk::mesh::comm_mesh_counts(get_bulk(), counts);
    num_global_elements = counts[stk::topology::ELEM_RANK];
  }

  //////////////////////////////////////////////////////////////////////////////////////

  void fill_vertex_weight_field_data()
  {
    std::vector<double> allWeights = get_all_element_weights();
    putVertexWeightsOnFieldForLocalElements(allWeights);
  }

  std::vector<double> get_all_element_weights() const
  {
    std::vector<double> weights(get_num_global_elements(), 1);
    for(size_t i=0;i<weights.size();++i)
      weights[i] = i+1;
    return weights;
  }

  void putVertexWeightsOnFieldForLocalElements(const std::vector<double>& vertex_weights_for_all_elements)
  {
    const stk::mesh::BucketVector &buckets = get_bulk().get_buckets(stk::topology::ELEMENT_RANK, get_meta().locally_owned_part());
    for(size_t i = 0; i < buckets.size(); i++)
      set_vertex_weights_for_bucket(*buckets[i], vertex_weights_for_all_elements);
  }

  void set_vertex_weights_for_bucket(stk::mesh::Bucket& bucket, const std::vector<double>& vector_of_data)
  {
    for(size_t j = 0; j < bucket.size(); j++)
      set_vertex_weight_for_local_element(vector_of_data, bucket[j]);
  }

  void set_vertex_weight_for_local_element(const std::vector<double> &vector_of_data, stk::mesh::Entity element)
  {
    double*vertex_weight = stk::mesh::field_data(*vertexWeightField, element);
    unsigned id = get_bulk().identifier(element);
    *vertex_weight = vector_of_data[id-1];
  }

  //////////////////////////////////////////////////////////////////////////////////////

  void check_graph_vertex_weights(const std::vector<BalanceGlobalNumber>& vertexIds, const std::vector<double>& vertexWeights)
  {
    for(size_t i=0; i<vertexIds.size(); ++i) {
      stk::mesh::Entity elem = get_bulk().get_entity(stk::topology::ELEM_RANK, vertexIds[i]);
      double*vertex_weight = stk::mesh::field_data(*vertexWeightField, elem);
      EXPECT_EQ(*vertex_weight, vertexWeights[i]);
    }
  }

  size_t get_num_global_elements() const
  {
    STK_ThrowRequireWithSierraHelpMsg(num_global_elements>0);
    return num_global_elements;
  }

  stk::balance::DoubleFieldType * vertexWeightField;
  size_t num_global_elements = 0;
};


TEST_F(VertexWeightSettings, checkVertexWeights)
{
  if(stk::parallel_machine_size(get_comm()) == 2)
  {
    setup_mesh_with_global_element_count_set("generated:1x1x6", stk::mesh::BulkData::AUTO_AURA);
    test_setting_of_vertex_weights_via_field();
  }
}

}
