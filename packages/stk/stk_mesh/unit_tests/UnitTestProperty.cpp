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
#include <stk_mesh/baseImpl/PartRepository.hpp>
#include <stk_mesh/baseImpl/EntityRepository.hpp>

namespace stk {
namespace mesh {

class UnitTestProperty {
public:
  UnitTestProperty() {}
  void testProperty();
};

}//namespace mesh
}//namespace stk

namespace {

STKUNIT_UNIT_TEST(UnitTestProperty, testUnit)
{
  stk::mesh::UnitTestProperty umeta;
  umeta.testProperty();
}

}//namespace <anonymous>

namespace stk {
namespace mesh {

void UnitTestProperty::testProperty()
{
  std::vector< std::string > dummy_names(1);
  dummy_names[0].assign("dummy");

  stk::mesh::MetaData meta_data( dummy_names );
  stk::mesh::MetaData meta_data2( dummy_names );

  stk::mesh::Property<int> & pi = meta_data.declare_property<int>("my_i");
  stk::mesh::Property<double> & pf = meta_data.declare_property<double>("my_d");
  stk::mesh::Property<double> & px = meta_data.declare_property<double>("my_x");
  stk::mesh::Property<double> & pa = meta_data.declare_property<double>("my_array",5);

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

  meta_data.put_property( pi , meta_data.locally_owned_part() );

  STKUNIT_ASSERT( stk::mesh::property_data( pi , meta_data.locally_owned_part() ) != NULL);
  STKUNIT_ASSERT( ! stk::mesh::property_data( px , meta_data.locally_owned_part() ) );

  stk::mesh::Property<int> * property();

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
  stk::mesh::Part & part_not_equal_to_pi = meta_data2.declare_part( "part_not_equal_to_pi", 0 );
  STKUNIT_ASSERT_THROW(stk::mesh::property_data( pi , part_not_equal_to_pi ),
                       std::logic_error);

  //Final coverage of MetaData.hpp - declare_property
  stk::mesh::Property<double> & pb = meta_data.declare_property<double>("my_y",0);
  const std::string& string_correct_new_double = "my_y";
  meta_data.get_property<double>( string_correct_new_double ); 
  STKUNIT_ASSERT( (pb).type_is<double>() );

  //Coverage of add_property in Property.hpp
  //  const std::string    & name = "";
  //const std::type_info & type = pi.type_is<int>();
  //unsigned         size = 5;

  /* const std::vector< PropertyBase* > properties = const std::vector< PropertyBase* > &get_ properties();

  properties.push_back( pi );
  properties.push_back( pf );
  properties.push_back( px );
  properties.push_back( pd );
  properties.push_back( pb );

  std::vector<PropertyBase * >::iterator j = properties.begin();

  for ( ; j != m_properties.end() ; ++j ) { 
   stk::mesh::Part & part = meta_data.declare_part( "part", 0 );
   meta_data.put_property( *j, part);
  }*/

}

}
}
