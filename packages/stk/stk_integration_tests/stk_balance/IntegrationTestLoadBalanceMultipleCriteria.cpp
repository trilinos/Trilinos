#include <stk_mesh/base/GetEntities.hpp>
#include <stk_unit_test_utils/unittestMeshUtils.hpp>

#include <stk_balance/balance.hpp>
#include <stk_balance/balanceUtils.hpp>
#include <stk_balance/internal/privateDeclarations.hpp>

#include <stk_unit_test_utils/MeshFixture.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_io/StkMeshIoBroker.hpp>

namespace
{
void set_criteria_on_mesh(std::vector<const stk::mesh::Field<double>*> critFields, const stk::mesh::BulkData& bulkData)
{
  const stk::mesh::BucketVector& buckets = bulkData.buckets(stk::topology::ELEM_RANK);
  for(size_t i=0;i<buckets.size();++i)
  {
    const stk::mesh::Bucket& bucket = *buckets[i];
    double* data = stk::mesh::field_data(*(critFields[bulkData.parallel_rank()]), bucket);
    for(size_t j=0;j<bucket.size();++j)
    {
      data[j] = 1.0 + bulkData.parallel_rank();
    }
  }
}

void write_mesh_and_results(stk::mesh::BulkData& bulkData, const std::string& filename, std::vector<stk::mesh::Field<double>*> out_fields)
{
  stk::io::StkMeshIoBroker broker(bulkData.parallel());
  broker.set_bulk_data(bulkData);
  size_t fh = broker.create_output_mesh(filename, stk::io::WRITE_RESULTS);
  for(size_t i=0;i<out_fields.size();++i)
    broker.add_field(fh, *out_fields[i]);
  broker.begin_output_step(fh, 0);
  broker.write_defined_output_fields(fh);
  broker.end_output_step(fh);
}

class MultipleCriteria : public stk::unit_test_util::MeshFixture
{
public:

  void put_decomposition_on_mesh(stk::mesh::Field<double>& procId)
  {
    stk::mesh::EntityVector elements;
    stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part(), elements);
    for(stk::mesh::Entity element : elements)
    {
      double* data = stk::mesh::field_data(procId, element);
      *data = get_bulk().parallel_rank();
    }
  }

  void test_multiple_criteria_rebalance(const std::string& decompMethod)
  {
    if(get_parallel_size() == 2)
    {
      setup_empty_mesh(stk::mesh::BulkData::AUTO_AURA);
      stk::mesh::Field<double> &criteria1 = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "criteria1");
      stk::mesh::Field<double> &criteria2 = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "criteria2");
      stk::mesh::Field<double> &procId = get_meta().declare_field<double>(stk::topology::ELEM_RANK, "proc_id");

      stk::mesh::put_field_on_mesh(criteria1, get_meta().universal_part(), nullptr);
      stk::mesh::put_field_on_mesh(criteria2, get_meta().universal_part(), nullptr);
      stk::mesh::put_field_on_mesh(procId, get_meta().universal_part(),    nullptr);

      stk::io::fill_mesh("generated:4x4x4", get_bulk());

      std::vector<const stk::mesh::Field<double>*> critFields = {&criteria1, &criteria2};

      set_criteria_on_mesh(critFields, get_bulk());

      put_decomposition_on_mesh(procId);

      write_mesh_and_results(get_bulk(), "before.e", std::vector<stk::mesh::Field<double>*>{&criteria1, &criteria2, &procId});

      std::ostringstream os;
      stk::mesh::EntityVector elements;
      stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part(), elements);
      size_t weight1 = 0;
      size_t weight2 = 0;

      for(stk::mesh::Entity element : elements)
      {
        os << "Before - Proc " << get_bulk().parallel_rank() << " has element " << get_bulk().identifier(element) << std::endl;
        weight1 += *stk::mesh::field_data(criteria1, element);
        weight2 += *stk::mesh::field_data(criteria2, element);
      }

      os << "Before - Proc " << get_bulk().parallel_rank() << " Sum weight1: " << weight1 << " and weight2: " << weight2 << std::endl;

      stk::balance::MultipleCriteriaSettings graphSettings(critFields);
      graphSettings.setDecompMethod(decompMethod);
      os << "decomp method: " << graphSettings.getDecompMethod() << std::endl;
      stk::balance::balanceStkMesh(graphSettings, get_bulk());

      stk::mesh::get_entities(get_bulk(), stk::topology::ELEM_RANK, get_meta().locally_owned_part(), elements);

      weight1 = 0;
      weight2 = 0;
      for(stk::mesh::Entity element : elements)
      {
        os << "After - Proc " << get_bulk().parallel_rank() << " has element " << get_bulk().identifier(element) << std::endl;
        weight1 += *stk::mesh::field_data(criteria1, element);
        weight2 += *stk::mesh::field_data(criteria2, element);
      }

      os << "After - Proc " << get_bulk().parallel_rank() << " Sum weight1: " << weight1 << " and weight2: " << weight2 << std::endl;
      std::cerr << os.str();

      std::vector<size_t> counts;
      stk::mesh::Selector sel = get_meta().locally_owned_part();
      stk::mesh::comm_mesh_counts(get_bulk(), counts, &sel);

      size_t gold_weight1 = counts[stk::topology::ELEM_RANK]/2;
      size_t gold_weight2 = 2*gold_weight1;

      gold_weight1 /= get_bulk().parallel_size();
      gold_weight2 /= get_bulk().parallel_size();

      EXPECT_EQ(gold_weight1, weight1);
      EXPECT_EQ(gold_weight2, weight2);

      put_decomposition_on_mesh(procId);

      write_mesh_and_results(get_bulk(), "after.e", std::vector<stk::mesh::Field<double>*>{&criteria1, &criteria2, &procId});
    }
  }
};


TEST_F(MultipleCriteria, using_fields_rcb)
{
  test_multiple_criteria_rebalance("rcb");
}

TEST_F(MultipleCriteria, using_fields_parmetis)
{
  test_multiple_criteria_rebalance("parmetis");
}


TEST(ThreeDimArray, getIndexOfFlatArray)
{
  const size_t num_selectors = 5;
  const size_t num_entities = 8;
  const size_t num_criteria = 3;
  size_t counter = 0;

  // num_weights vertex 0, num_weights vertex 1, ...

  for(size_t i=0;i<num_entities;++i)
  {
    for(size_t j=0;j<num_selectors;++j)
    {
      for(size_t k=0;k<num_criteria;++k)
      {
        EXPECT_EQ(counter, stk::balance::internal::get_index(num_selectors, num_criteria, i, j, k));
        counter++;
      }
    }
  }
}

}
