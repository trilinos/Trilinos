
#include <unit_tests/stk_utest_macros.hpp>

#include <stk_mesh/fem/EntityTypes.hpp>
#include <stk_mesh/fem/FieldDeclarations.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <stk_linsys/FieldIdMap.hpp>
#include <stk_linsys/ImplDetails.hpp>

namespace stk_linsys_unit_tests {

//------------- here is the unit-test of stk::linsys::impl functions... -----------------------

void testImpl( MPI_Comm comm )
{
  //First create and initialize a MetaData to use in our testing:

  stk::mesh::MetaData meta_data( stk::mesh::fem_entity_type_names() );

  const unsigned number_of_states = 1;

  stk::mesh::ScalarField& temperature_field =
    meta_data.declare_field<stk::mesh::ScalarField>( "temperature", number_of_states );
  stk::mesh::ScalarField& pressure_field =
    meta_data.declare_field<stk::mesh::ScalarField>( "pressure", number_of_states );
  stk::mesh::VectorField& velocity_field =
    meta_data.declare_field<stk::mesh::VectorField>( "velocity",    number_of_states );

  meta_data.commit();


  //Now test the linsys::impl::map_field_to_int function:

  stk::linsys::FieldIdMap field_id_map;

  int fieldid = stk::linsys::impl::map_field_to_int(field_id_map, temperature_field);
  STKUNIT_ASSERT_EQUAL( fieldid, (int)0 );

  fieldid = stk::linsys::impl::query_field_to_int_mapping(field_id_map, temperature_field);
  STKUNIT_ASSERT_EQUAL( fieldid, (int)0 );

  //test query_field_to_int_mapping for the case where the specified field is not in the map.
  //(an exception should be thrown in this case)
  std::cout << "\nTesting query_field_to_int_mapping with un-mapped field..." << std::endl;
  bool threw_exc0 = false;
  try {
    fieldid = stk::linsys::impl::query_field_to_int_mapping(field_id_map, pressure_field);
  }
  catch(std::exception& exc) {
    threw_exc0 = true;
    std::cout << "Caught expected exception: " << exc.what() << std::endl;
  }
  STKUNIT_ASSERT_EQUAL( threw_exc0, (bool)true );

  //now add the pressure field...
  fieldid = stk::linsys::impl::map_field_to_int(field_id_map, pressure_field);
  STKUNIT_ASSERT_EQUAL( fieldid, (int)1 );

  fieldid = stk::linsys::impl::map_field_to_int(field_id_map, velocity_field);
  STKUNIT_ASSERT_EQUAL( fieldid, (int)2 );

  //here the pressure field is already in the map.
  fieldid = stk::linsys::impl::map_field_to_int(field_id_map, pressure_field);
  STKUNIT_ASSERT_EQUAL( fieldid, (int)1 );

  //confirm that the map is the expected size:
  STKUNIT_ASSERT_EQUAL( field_id_map.size(), (size_t)3 );

  //test the linsys::impl::get_field function:
  const stk::mesh::FieldBase* field = stk::linsys::impl::get_field(field_id_map, 2);

  STKUNIT_ASSERT_EQUAL( field->name(), velocity_field.name() );

  //call linsys::impl::get_field with an unknown field-id to confirm that an
  //exception is thrown.
  bool threw_exc1 = false;
  try {
    field = stk::linsys::impl::get_field(field_id_map, 5);
  }
  catch(std::exception& exc) {
    threw_exc1 = true;
    std::cout << "Caught expected exception: " << exc.what() << std::endl;
  }
  STKUNIT_ASSERT_EQUAL( threw_exc1, (bool)true );

  //Now test the linsys::impl::entityid_to_int function:

  stk::mesh::EntityId id = 42;
  int int_id1 = stk::linsys::impl::entityid_to_int(id);

  STKUNIT_ASSERT_EQUAL( int_id1, (int)42 );

  //3 billion is too large to represent as int...
  uint64_t big_id = 3000000000u;
  STKUNIT_ASSERT_THROW( stk::linsys::impl::verify_convertible_to_int(big_id, "linsys unit-test"), std::runtime_error);
}

} // namespace stk_linsys_unit_tests

STKUNIT_UNIT_TEST(UnitTestingOfLinsysImplementation, testUnit)
{
  MPI_Barrier( MPI_COMM_WORLD );
  stk_linsys_unit_tests::testImpl ( MPI_COMM_WORLD );
}

