/*------------------------------------------------------------------------*/
/*                 Copyright (c) 2013, Sandia Corporation.
/*                 Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
/*                 the U.S. Governement retains certain rights in this software.
/*                 
/*                 Redistribution and use in source and binary forms, with or without
/*                 modification, are permitted provided that the following conditions are
/*                 met:
/*                 
/*                     * Redistributions of source code must retain the above copyright
/*                       notice, this list of conditions and the following disclaimer.
/*                 
/*                     * Redistributions in binary form must reproduce the above
/*                       copyright notice, this list of conditions and the following
/*                       disclaimer in the documentation and/or other materials provided
/*                       with the distribution.
/*                 
/*                     * Neither the name of Sandia Corporation nor the names of its
/*                       contributors may be used to endorse or promote products derived
/*                       from this software without specific prior written permission.
/*                 
/*                 THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
/*                 "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
/*                 LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
/*                 A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
/*                 OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
/*                 SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
/*                 LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
/*                 DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
/*                 THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
/*                 (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
/*                 OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
/*                 
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/


#include <stddef.h>                     // for NULL
#include <stdexcept>                    // for logic_error, runtime_error
#include <stk_mesh/base/MetaData.hpp>   // for MetaData
#include <stk_mesh/base/Property.hpp>   // for property_data
#include <gtest/gtest.h>
#include <string>                       // for string
#include "stk_mesh/base/PropertyBase.hpp"  // for Property, etc
#include "stk_topology/topology.hpp"    // for topology, etc
namespace stk { namespace mesh { class Part; } }




using stk::mesh::MetaData;
using stk::mesh::Part;
using stk::mesh::Property;

namespace {

TEST(UnitTestProperty, testProperty)
{
  const int spatial_dimension = 3;
  MetaData meta_data(spatial_dimension);
  MetaData meta_data2(spatial_dimension);
  unsigned size = 0;

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
  Part & part_not_equal_to_pi2_2 = meta_data2.declare_part( "part_not_equal_to_pi2_2", stk::topology::NODE_RANK );
  Part & part = meta_data.declare_part( "part", stk::topology::NODE_RANK );

  meta_data.commit();
  meta_data2.commit();

  ASSERT_TRUE( pi.type_is<int>() );
  ASSERT_TRUE( pf.type_is<double>() );
  ASSERT_TRUE( px.type_is<double>() );
  ASSERT_TRUE( pa.type_is<double>() );

  ASSERT_TRUE( ! pi.type_is<double>() );
  ASSERT_TRUE( ! pa.type_is<int>() );

  ASSERT_EQ( pi.size() , 1u );
  ASSERT_EQ( pf.size() , 1u );
  ASSERT_EQ( px.size() , 1u );
  ASSERT_EQ( pa.size() , 5u );

  //Covers add_property in Property.hpp and put_property in MetaData.hpp
  meta_data.put_property( pi , meta_data.locally_owned_part() );
  //covers property_data in Property.hpp
  ASSERT_TRUE( stk::mesh::property_data( pi , meta_data.locally_owned_part() ) != NULL);
  ASSERT_TRUE( !stk::mesh::property_data( px , meta_data.locally_owned_part() ) );

  //Coverage of virtual const data_type * data( unsigned key ) const in Property.hpp
  ASSERT_FALSE( stk::mesh::property_data( pi2 , meta_data.locally_owned_part() ) != NULL);

  //Cover unsigned data type in Property.hpp
  ASSERT_FALSE( stk::mesh::property_data( pi_unsigned , meta_data.locally_owned_part() ) != NULL);
  ASSERT_FALSE( stk::mesh::property_data( pi_unsigned_2 , meta_data.locally_owned_part() ) != NULL);
  ASSERT_TRUE( ! pi_unsigned.type_is<int>() );

  //Cover unsigned const data type in Property.hpp
  ASSERT_FALSE( stk::mesh::property_data( pi_unsigned_const , meta_data.locally_owned_part() ) != NULL);
  ASSERT_FALSE( stk::mesh::property_data( pi_unsigned_const_2 , meta_data.locally_owned_part() ) != NULL);

  //Coverage of Property.hpp to test two unequal parts for const property_type and cover of property_data_throw
  ASSERT_THROW(property_data( pi2 , part_not_equal_to_pi2 ),
		       std::logic_error);
  ASSERT_THROW(property_data( pi2 , part_not_equal_to_pi2_2 ),
		       std::logic_error);

  //Test get_property_base with an correct type in Property.cpp
  const std::string& string_correct_double = "my_d";
  meta_data.get_property<double>( string_correct_double );

  //Test get_property_base with an incorrect type and size in Property.cpp
  const std::string& string_incorrect_double = "my_i";
  ASSERT_THROW(
      meta_data.get_property<double>( string_incorrect_double ),
      std::runtime_error
      );

  //More coverage of Property.hpp to test two parts with different meta data
  ASSERT_THROW(stk::mesh::property_data( pi , part_not_equal_to_pi ),
                       std::logic_error);

  //Final coverage of MetaData.hpp - declare_property
  const std::string& string_correct_new_double = "my_y";
  meta_data.get_property<double>( string_correct_new_double );
  ASSERT_TRUE( (pb).type_is<double>() );

  //Coverage of add_property in Property.hpp
  meta_data.put_property( pProp, part);

  //Coverage of declare_property in MetaData.hpp
  ASSERT_TRUE( pb.type_is<double>() );

  //Further coverage of declare_property in MetaData.hpp ( pv != NULL)
  ASSERT_TRUE( pi3.type_is<int>() );

}
}
