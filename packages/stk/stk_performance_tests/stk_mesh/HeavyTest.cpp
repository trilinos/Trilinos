
#include <stk_io/util/IO_Fixture.hpp>
#include <stk_mesh/fem/SkinMesh.hpp>
//#include <stk_rebalance/Rebalance.hpp>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <string>

#if 0
STKUNIT_UNIT_TEST( HeavyTest, heavytest) {

  stk::io::util::IO_Fixture fixture;

  // time meta_data initialize
  std::string input_base_exodus_filename;

  fixture.parallel_initialize_meta_data( input_base_exodus_filename );
  stk::mesh::fem::FEMMetaData & meta_data = fixture.meta_data();

  meta_data.commit();

  // time bulk_data initialize
  fixture.parallel_initialize_bulk_data( input_base_exodus_filename );
  stk::mesh::BulkData & bulk_data = fixture.bulk_data();

  // time rebalance bulk_data

  // time skin bulk_data
  stk::mesh::skin_mesh( bulk_data, meta_data.spatial_dimension());


  // time exodus file creation
  std::string output_base_exodus_filename;
  fixture.create_output_mesh( output_base_exodus_filename );

  double time_step = 0;
  fixture.process_output_mesh( time_step );

}
#endif
