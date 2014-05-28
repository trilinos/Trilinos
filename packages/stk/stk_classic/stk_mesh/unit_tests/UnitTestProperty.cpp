/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <sstream>
#include <stdexcept>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Property.hpp>
#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/fem/FEMHelpers.hpp>

using stk::mesh::MetaData;
using stk::mesh::Part;
using stk::mesh::Property;

namespace {

STKUNIT_UNIT_TEST(UnitTestProperty, testProperty)
{
  const int spatial_dimension = 3;
  MetaData meta_data(stk::mesh::fem::entity_rank_names(spatial_dimension));
  MetaData meta_data2(stk::mesh::fem::entity_rank_names(spatial_dimension));
  unsigned size = 0;
  unsigned EntityRank = 0;

  Property<int> & pi = meta_data.declare_property<int>("my_i");
  Property<double> & pf = meta_data.declare_property<double>("my_d");
  Property<double> & px = meta_data.declare_property<double>("my_x");
  Property<double> & pa = meta_data.declare_property<double>("my_array",5);

  const Property<int> & pi2 = meta_data.declare_property<int>("my_i5");
  Property<double> & pb = meta_data.declare_property<double>("my_y",0);
  Property<int> & pi3 = meta_data.declare_property<int>("my_i", pi.size());
  Property<unsigned> & pi_unsigned = meta_data.declare_property<unsigned>("string_unsigned", size);
  const Property<unsigned> & pi_unsigned_const = meta_data.declare_property<unsigned>("string_unsigned", size);
  Property<unsigned> & pi_unsigned_2 = meta_data.declare_property<unsigned>("string_unsigned_2");
  const Property<unsigned> & pi_unsigned_const_2 = meta_data.declare_property<unsigned>("string_unsigned_2");
  Property<double> & pProp = meta_data.declare_property<double>("string_correct_new_double", pb.size());
  Part & part_not_equal_to_pi2 =  meta_data2.declare_part( "part_not_equal_to_pi2" );
  Part & part_not_equal_to_pi = meta_data2.declare_part( "part_not_equal_to_pi");
  Part & part_not_equal_to_pi2_2 = meta_data2.declare_part( "part_not_equal_to_pi2_2", EntityRank );
  Part & part = meta_data.declare_part( "part", EntityRank );

  meta_data.commit();
  meta_data2.commit();

  STKUNIT_ASSERT( pi.type_is<int>() );
  STKUNIT_ASSERT( pf.type_is<double>() );
  STKUNIT_ASSERT( px.type_is<double>() );
  STKUNIT_ASSERT( pa.type_is<double>() );

  STKUNIT_ASSERT( ! pi.type_is<double>() );
  STKUNIT_ASSERT( ! pa.type_is<int>() );

  STKUNIT_ASSERT_EQUAL( pi.size() , 1u );
  STKUNIT_ASSERT_EQUAL( pf.size() , 1u );
  STKUNIT_ASSERT_EQUAL( px.size() , 1u );
  STKUNIT_ASSERT_EQUAL( pa.size() , 5u );

  //Covers add_property in Property.hpp and put_property in MetaData.hpp
  meta_data.put_property( pi , meta_data.locally_owned_part() );
  //covers property_data in Property.hpp
  STKUNIT_ASSERT( stk::mesh::property_data( pi , meta_data.locally_owned_part() ) != NULL);
  STKUNIT_ASSERT( !stk::mesh::property_data( px , meta_data.locally_owned_part() ) );

  //Coverage of virtual const data_type * data( unsigned key ) const in Property.hpp
  STKUNIT_ASSERT_FALSE( stk::mesh::property_data( pi2 , meta_data.locally_owned_part() ) != NULL);

  //Cover unsigned data type in Property.hpp
  STKUNIT_ASSERT_FALSE( stk::mesh::property_data( pi_unsigned , meta_data.locally_owned_part() ) != NULL);
  STKUNIT_ASSERT_FALSE( stk::mesh::property_data( pi_unsigned_2 , meta_data.locally_owned_part() ) != NULL);
  STKUNIT_ASSERT( ! pi_unsigned.type_is<int>() );

  //Cover unsigned const data type in Property.hpp
  STKUNIT_ASSERT_FALSE( stk::mesh::property_data( pi_unsigned_const , meta_data.locally_owned_part() ) != NULL);
  STKUNIT_ASSERT_FALSE( stk::mesh::property_data( pi_unsigned_const_2 , meta_data.locally_owned_part() ) != NULL);

  //Coverage of Property.hpp to test two unequal parts for const property_type and cover of property_data_throw
  STKUNIT_ASSERT_THROW(property_data( pi2 , part_not_equal_to_pi2 ),
		       std::logic_error);
  STKUNIT_ASSERT_THROW(property_data( pi2 , part_not_equal_to_pi2_2 ),
		       std::logic_error);

  //Test get_property_base with an correct type in Property.cpp
  const std::string& string_correct_double = "my_d";
  meta_data.get_property<double>( string_correct_double );

  //Test get_property_base with an incorrect type and size in Property.cpp
  const std::string& string_incorrect_double = "my_i";
  STKUNIT_ASSERT_THROW(
      meta_data.get_property<double>( string_incorrect_double ),
      std::runtime_error
      );

  //More coverage of Property.hpp to test two parts with different meta data
  STKUNIT_ASSERT_THROW(stk::mesh::property_data( pi , part_not_equal_to_pi ),
                       std::logic_error);

  //Final coverage of MetaData.hpp - declare_property
  const std::string& string_correct_new_double = "my_y";
  meta_data.get_property<double>( string_correct_new_double );
  STKUNIT_ASSERT( (pb).type_is<double>() );

  //Coverage of add_property in Property.hpp
  meta_data.put_property( pProp, part);

  //Coverage of declare_property in MetaData.hpp
  STKUNIT_ASSERT( pb.type_is<double>() );

  //Further coverage of declare_property in MetaData.hpp ( pv != NULL)
  STKUNIT_ASSERT( pi3.type_is<int>() );

}
}
