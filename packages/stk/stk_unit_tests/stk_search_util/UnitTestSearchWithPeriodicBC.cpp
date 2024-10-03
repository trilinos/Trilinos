// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
// 
//     * Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
// 
//     * Redistributions in binary form must reproduce the above
//       copyright notice, this list of conditions and the following
//       disclaimer in the documentation and/or other materials provided
//       with the distribution.
// 
//     * Neither the name of NTESS nor the names of its contributors
//       may be used to endorse or promote products derived from this
//       software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
// OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
// DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
// THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
// (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
// OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
// 

#include "stk_unit_test_utils/stk_mesh_fixtures/HexFixture.hpp"
#include <cmath>
#include <gtest/gtest.h>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <stk_mesh/base/Field.hpp>
#include <stk_search_util/PeriodicBoundarySearch.hpp>

typedef stk::mesh::fixtures::HexFixture::CoordFieldType CoordFieldType;
typedef stk::mesh::GetCoordinates<CoordFieldType> CoordinateFunctor;
typedef stk::mesh::PeriodicBoundarySearch<CoordinateFunctor> PeriodicSearch;

namespace {
const double PI     = M_PI;
const double TWO_PI = 2 * PI;

void expect_eq_for_shared_or_owned_node(const stk::mesh::BulkData & bulk_data, stk::mesh::Entity node,
                                        const stk::mesh::Field<double> & theField, double expected_value)
{
  if (bulk_data.is_valid(node) && (bulk_data.bucket(node).owned() || bulk_data.bucket(node).shared() ) )
  {
    const double * const vol = stk::mesh::field_data(theField, node);
    EXPECT_EQ(*vol, expected_value);
  }
}

void do_periodic_assembly(stk::mesh::BulkData & bulk_data, PeriodicSearch & pbc_search,
                          stk::mesh::Field<double> & volField)
{
  //gather to domain node from possibly multiple ranges
  for (size_t i = 0; i < pbc_search.size(); ++i)
  {
    std::pair<stk::mesh::Entity, stk::mesh::Entity> entityPair = pbc_search.get_node_pair(i);
    STK_ThrowRequire(bulk_data.is_valid(entityPair.first));
    STK_ThrowRequire(bulk_data.is_valid(entityPair.second));
    double * domainField = stk::mesh::field_data(volField, entityPair.first);
    double * rangeField = stk::mesh::field_data(volField, entityPair.second);
    *domainField += *rangeField;
  }

  std::vector< stk::mesh::FieldBase const * > ghosted_fields;
  ghosted_fields.push_back(&volField);
  stk::mesh::communicate_field_data( pbc_search.get_ghosting(), ghosted_fields);

  //set ranges equal to domain value
  for (size_t i = 0; i < pbc_search.size(); ++i)
  {
    std::pair<stk::mesh::Entity, stk::mesh::Entity> entityPair = pbc_search.get_node_pair(i);
    STK_ThrowRequire(bulk_data.is_valid(entityPair.first));
    STK_ThrowRequire(bulk_data.is_valid(entityPair.second));
    double * domainField = stk::mesh::field_data(volField, entityPair.first);
    double * rangeField = stk::mesh::field_data(volField, entityPair.second);
    *rangeField = *domainField;
  }
}

void do_volume_assembly(stk::mesh::BulkData & bulk_data, stk::mesh::Field<double> & volField)
{
  const stk::mesh::BucketVector & elemBuckets = bulk_data.buckets(stk::topology::ELEMENT_RANK);
  for (size_t elemBucket = 0; elemBucket < elemBuckets.size(); ++elemBucket)
  {
    const stk::mesh::Bucket & bucket = *elemBuckets[elemBucket];
    for (size_t i = 0; i < bucket.size(); ++i)
    {
      stk::mesh::Entity elem = bucket[i];

      const stk::mesh::Entity * aNode = bulk_data.begin_nodes(elem);
      const size_t numNodes = bulk_data.num_nodes(elem);
      for (size_t in = 0; in < numNodes; ++in)
      {
        double * vol = stk::mesh::field_data(volField, aNode[in]);
        *vol += 0.125;
      }
    }
  }

}

template<typename SearchPairVector>
void check_gold( const SearchPairVector & search_results )
  // check search result
{
  typedef std::vector<std::pair<stk::mesh::EntityId,stk::mesh::EntityId> > GoldVector;
  GoldVector gold;
  gold.push_back(std::make_pair(1,4));
  gold.push_back(std::make_pair(5,8));
  gold.push_back(std::make_pair(9,12));
  gold.push_back(std::make_pair(13,16));

  gold.push_back(std::make_pair(17,20));
  gold.push_back(std::make_pair(21,24));
  gold.push_back(std::make_pair(25,28));
  gold.push_back(std::make_pair(29,32));

  gold.push_back(std::make_pair(33,36));
  gold.push_back(std::make_pair(37,40));
  gold.push_back(std::make_pair(41,44));
  gold.push_back(std::make_pair(45,48));

  gold.push_back(std::make_pair(49,52));
  gold.push_back(std::make_pair(53,56));
  gold.push_back(std::make_pair(57,60));
  gold.push_back(std::make_pair(61,64));

  //make sure search result shows up in gold
  for (size_t i=0, size=search_results.size(); i<size; ++i) {
    stk::mesh::EntityId domain_node = search_results[i].first.id().id();
    stk::mesh::EntityId range_node = search_results[i].second.id().id();

    //entry in search is found in gold

    bool found = std::find(gold.begin(),gold.end(),std::make_pair(domain_node,range_node)) != gold.end();
    EXPECT_TRUE(found) << "We can't find domain/range: " << search_results[i].first << "/" << search_results[i].second  << std::endl;
  }
}

template<typename SearchPairVector>
void check_gold_two_way_multiperiodic( const SearchPairVector & search_results )
  // check search result
{
  typedef std::vector<std::pair<stk::mesh::EntityId,stk::mesh::EntityId> > GoldVector;
  GoldVector gold;
  gold.push_back(std::make_pair(1,4));
  gold.push_back(std::make_pair(5,8));
  gold.push_back(std::make_pair(9,12));

  gold.push_back(std::make_pair(17,20));
  gold.push_back(std::make_pair(21,24));
  gold.push_back(std::make_pair(25,28));

  gold.push_back(std::make_pair(33,36));
  gold.push_back(std::make_pair(37,40));
  gold.push_back(std::make_pair(41,44));

  gold.push_back(std::make_pair(49,52));
  gold.push_back(std::make_pair(53,56));
  gold.push_back(std::make_pair(57,60));

  //top and bottom
  gold.push_back(std::make_pair(1,13));
  gold.push_back(std::make_pair(2,14));
  gold.push_back(std::make_pair(3,15));

  gold.push_back(std::make_pair(17,29));
  gold.push_back(std::make_pair(18,30));
  gold.push_back(std::make_pair(19,31));

  gold.push_back(std::make_pair(33,45));
  gold.push_back(std::make_pair(34,46));
  gold.push_back(std::make_pair(35,47));

  gold.push_back(std::make_pair(49,61));
  gold.push_back(std::make_pair(50,62));
  gold.push_back(std::make_pair(51,63));

  //edge cases
  gold.push_back(std::make_pair(1,16));
  gold.push_back(std::make_pair(17,32));
  gold.push_back(std::make_pair(33,48));
  gold.push_back(std::make_pair(49,64));

  for (size_t i=0, size=search_results.size(); i<size; ++i) {
    stk::mesh::EntityId domain_node = search_results[i].first.id().id();
    stk::mesh::EntityId range_node = search_results[i].second.id().id();

    EXPECT_TRUE((std::find(gold.begin(), gold.end(), std::make_pair(domain_node,range_node) ) ) != gold.end());
  }
}

template<typename SearchPairVector>
void check_gold_three_way_multiperiodic( const SearchPairVector & search_results )
  // check search result
{
  typedef std::vector<std::pair<stk::mesh::EntityId,stk::mesh::EntityId> > GoldVector;
  GoldVector gold;
  gold.push_back(std::make_pair(1,  4));
  gold.push_back(std::make_pair(5,  8));
  gold.push_back(std::make_pair(9, 12));
  gold.push_back(std::make_pair(17, 20));
  gold.push_back(std::make_pair(21, 24));
  gold.push_back(std::make_pair(25, 28));
  gold.push_back(std::make_pair(33, 36));
  gold.push_back(std::make_pair(37, 40));
  gold.push_back(std::make_pair(41, 44));
  gold.push_back(std::make_pair(1, 13));
  gold.push_back(std::make_pair(2, 14));
  gold.push_back(std::make_pair(3, 15));
  gold.push_back(std::make_pair(17, 29));
  gold.push_back(std::make_pair(18, 30));
  gold.push_back(std::make_pair(19, 31));
  gold.push_back(std::make_pair(33, 45));
  gold.push_back(std::make_pair(34, 46));
  gold.push_back(std::make_pair(35, 47));
  gold.push_back(std::make_pair(1, 49));
  gold.push_back(std::make_pair(2, 50));
  gold.push_back(std::make_pair(3, 51));
  gold.push_back(std::make_pair(5, 53));
  gold.push_back(std::make_pair(6, 54));
  gold.push_back(std::make_pair(7, 55));
  gold.push_back(std::make_pair(9, 57));
  gold.push_back(std::make_pair(10, 58));
  gold.push_back(std::make_pair(11, 59));
  gold.push_back(std::make_pair(1, 16));
  gold.push_back(std::make_pair(17, 32));
  gold.push_back(std::make_pair(33, 48));
  gold.push_back(std::make_pair(1, 61));
  gold.push_back(std::make_pair(2, 62));
  gold.push_back(std::make_pair(3, 63));
  gold.push_back(std::make_pair(1, 52));
  gold.push_back(std::make_pair(5, 56));
  gold.push_back(std::make_pair(9, 60));
  gold.push_back(std::make_pair(1, 64));

  for (size_t i=0, size=search_results.size(); i<size; ++i) {
    stk::mesh::EntityId domain_node = search_results[i].first.id().id();
    stk::mesh::EntityId range_node = search_results[i].second.id().id();

    EXPECT_TRUE((std::find(gold.begin(), gold.end(), std::make_pair(domain_node,range_node) ) ) != gold.end());
  }
}

void check_single_periodic_assembly(const stk::mesh::BulkData & bulk_data,
                                    const stk::mesh::fixtures::HexFixture & fixture,
                                    const stk::mesh::Field<double> & volField,
                                    unsigned x,
                                    unsigned y,
                                    unsigned z)
{
  //interior of domain should be 1.0
  for (unsigned i=0; i<x+1u; ++i) {
    for (unsigned j=1; j<y; ++j) {
      for (unsigned k=1; k<z; ++k) {
        stk::mesh::Entity node = fixture.node(i,j,k);
        expect_eq_for_shared_or_owned_node(bulk_data, node, volField, 1.0);
      }
    }
  }

  //faces (not edges) should be 0.5
  //edges should be 0.25
  //there are no "corners" since it is periodic in the x direction
  for (unsigned i=0; i<x+1u; ++i) {
    //top and bottom
    for (unsigned k=1; k<z; ++k) {
      const unsigned jTop = 3;
      const unsigned jBot = 0;
      stk::mesh::Entity node = fixture.node(i,jTop,k);
      expect_eq_for_shared_or_owned_node(bulk_data, node, volField, 0.5);
      node = fixture.node(i,jBot,k);
      expect_eq_for_shared_or_owned_node(bulk_data, node, volField, 0.5);
    }
    //front and back
    for (unsigned j=1; j<y; ++j) {
      const unsigned kFront = 0;
      const unsigned kBack = 3;
      stk::mesh::Entity node = fixture.node(i,j,kFront);
      expect_eq_for_shared_or_owned_node(bulk_data, node, volField, 0.5);
      node = fixture.node(i,j,kBack);
      expect_eq_for_shared_or_owned_node(bulk_data, node, volField, 0.5);
    }
    //edges
    stk::mesh::Entity node = fixture.node(i, 0, 0);
    expect_eq_for_shared_or_owned_node(bulk_data, node, volField, 0.25);
    node = fixture.node(i, 0, 3);
    expect_eq_for_shared_or_owned_node(bulk_data, node, volField, 0.25);
    node = fixture.node(i, 3, 0);
    expect_eq_for_shared_or_owned_node(bulk_data, node, volField, 0.25);
    node = fixture.node(i, 3, 3);
    expect_eq_for_shared_or_owned_node(bulk_data, node, volField, 0.25);
  }

}

void check_rotation_matrix(const PeriodicSearch & pbc_search, double rotationAngle)
{
  //stored column major
  std::vector<double> matrix(9);
  pbc_search.get_search_row_major_rotation(matrix);
  const double tol = 1.0e-12;
  EXPECT_NEAR(matrix[0], std::cos(rotationAngle), tol);
  EXPECT_NEAR(matrix[1], -std::sin(rotationAngle), tol);
  EXPECT_NEAR(matrix[2], 0.0, tol);

  EXPECT_NEAR(matrix[3], std::sin(rotationAngle), tol);
  EXPECT_NEAR(matrix[4], std::cos(rotationAngle), tol);
  EXPECT_NEAR(matrix[5], 0.0, tol);

  EXPECT_NEAR(matrix[6], 0.0, tol);
  EXPECT_NEAR(matrix[7], 0.0, tol);
  EXPECT_NEAR(matrix[8], 1.0, tol);
}

class CylindricalCoordinateMappingWithOffset : public stk::mesh::fixtures::CoordinateMapping
{
public:
  typedef double Scalar;
  CylindricalCoordinateMappingWithOffset(Scalar radius, Scalar theta,
      unsigned numTheta, Scalar offsetX, Scalar offsetY, Scalar offsetZ)
      : CoordinateMapping(),
        m_radius(radius),
        m_theta(theta),
        m_numTheta(numTheta),
        m_offsetX(offsetX),
        m_offsetY(offsetY),
        m_offsetZ(offsetZ)
  { }
  virtual void getNodeCoordinates(Scalar * field, const size_t nx, const size_t ny, const size_t nz) const
  {
    Scalar fracTheta = nx/(m_numTheta - 1);

    //we want the angle to go from pi/2 to pi/2 - theta so we do not
    //invert any elements
    Scalar angle = PI/2.0 + m_theta*fracTheta;
    field[0] = (m_radius + ny)*std::cos(angle) + m_offsetX;
    field[1] = (m_radius + ny)*std::sin(angle) + m_offsetY;
    field[2] = nz + m_offsetZ;
  }
private:
  Scalar m_radius;
  Scalar m_theta;
  unsigned m_numTheta;
  Scalar m_offsetX;
  Scalar m_offsetY;
  Scalar m_offsetZ;


};

}// namespace

TEST(CoarseSearch, PeriodicBC)
{
  const unsigned x = 3, y = 3, z = 3;

  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, x, y, z);

  stk::mesh::BulkData & bulk_data = fixture.m_bulk_data;
  stk::mesh::MetaData & meta_data = fixture.m_meta;
  CoordFieldType & coords_field = *fixture.m_coord_field;

  stk::mesh::Part & side_0 = meta_data.declare_part("side_0", stk::topology::NODE_RANK);
  stk::mesh::Part & side_3 = meta_data.declare_part("side_3", stk::topology::NODE_RANK);

  stk::mesh::Field<double> & volField = meta_data.declare_field<double>(stk::topology::NODE_RANK, "volume");
  stk::mesh::put_field_on_mesh(volField, meta_data.universal_part(), nullptr);

  meta_data.commit();

  fixture.generate_mesh();

  stk::mesh::PartVector independent_parts(1,&side_0);
  stk::mesh::PartVector dependent_parts(1,&side_3);

  bulk_data.modification_begin();
  for (unsigned i=0; i<y+1u; ++i) {
  for (unsigned j=0; j<z+1u; ++j) {
    stk::mesh::Entity node_0 = fixture.node(0,i,j);
    if (bulk_data.is_valid(node_0)  && bulk_data.bucket(node_0).owned()) {
      bulk_data.change_entity_parts( fixture.node(0,i,j), independent_parts);
    }
    stk::mesh::Entity node_3 = fixture.node(3,i,j);
    if (bulk_data.is_valid(node_3)  && bulk_data.bucket(node_3).owned()) {
      bulk_data.change_entity_parts( fixture.node(3,i,j), dependent_parts);
    }
  }}
  bulk_data.modification_end();

  //do periodic search
  typedef stk::mesh::GetCoordinates<CoordFieldType> CoordinateFunctor;
  typedef stk::mesh::PeriodicBoundarySearch<CoordinateFunctor> PeriodicSearch;
  PeriodicSearch pbc_search(bulk_data, CoordinateFunctor(bulk_data, coords_field));

  const stk::mesh::Selector side_0_selector = side_0 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());
  const stk::mesh::Selector side_3_selector = side_3 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());

  pbc_search.add_linear_periodic_pair(side_0_selector, side_3_selector);
  pbc_search.find_periodic_nodes(bulk_data.parallel());


  check_gold(pbc_search.get_pairs() );

  bulk_data.modification_begin();
  pbc_search.create_ghosting("periodic_ghosts");
  bulk_data.modification_end();

  do_volume_assembly(bulk_data, volField);

  std::vector< stk::mesh::FieldBase const * > ghosted_fields;
  ghosted_fields.push_back(&volField);
  stk::mesh::communicate_field_data( pbc_search.get_ghosting(), ghosted_fields);

  do_periodic_assembly(bulk_data, pbc_search, volField);

  check_single_periodic_assembly(bulk_data, fixture, volField, x, y, z);
}


void assign_to_parts_for_two_way(const unsigned x, const unsigned y, const unsigned z,
                                 stk::mesh::fixtures::HexFixture &fixture,
                                 stk::mesh::BulkData &bulk_data,
                                 stk::mesh::PartVector &side_0_parts,
                                 stk::mesh::PartVector &side_1_parts,
                                 stk::mesh::PartVector &side_2_parts,
                                 stk::mesh::PartVector &side_3_parts)
{

  //add periodic pair for side 0 and side 2
  bulk_data.modification_begin();
  for (unsigned i=0; i<y+1u; ++i) {
  for (unsigned j=0; j<z+1u; ++j) {
    stk::mesh::Entity node_0 = fixture.node(0,i,j);
    if (bulk_data.is_valid(node_0)  && bulk_data.bucket(node_0).owned()) {
      bulk_data.change_entity_parts( fixture.node(0,i,j), side_0_parts);
    }
    stk::mesh::Entity node_2 = fixture.node(3,i,j);
    if (bulk_data.is_valid(node_2)  && bulk_data.bucket(node_2).owned()) {
      bulk_data.change_entity_parts( fixture.node(3,i,j), side_2_parts);
    }
  }}
  bulk_data.modification_end();

  //add periodic pair for side 1 and side 3
  bulk_data.modification_begin();
  for (unsigned i=0; i<y+1u; ++i) {
  for (unsigned j=0; j<z+1u; ++j) {
    stk::mesh::Entity node_1 = fixture.node(i,0,j);
    if (bulk_data.is_valid(node_1)  && bulk_data.bucket(node_1).owned()) {
      bulk_data.change_entity_parts( fixture.node(i,0,j), side_1_parts);
    }
    stk::mesh::Entity node_3 = fixture.node(i,3,j);
    if (bulk_data.is_valid(node_3)  && bulk_data.bucket(node_3).owned()) {
      bulk_data.change_entity_parts( fixture.node(i,3,j), side_3_parts);
    }
  }}
  bulk_data.modification_end();
}

TEST(CoarseSearch, TwoWayMultiPeriodicBC)
{
  const unsigned x = 3, y = 3, z = 3;


  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, x, y, z);

  stk::mesh::BulkData & bulk_data = fixture.m_bulk_data;
  stk::mesh::MetaData & meta_data = fixture.m_meta;
  CoordFieldType & coords_field = *fixture.m_coord_field;

  stk::mesh::Part & side_0 = meta_data.declare_part("side_0", stk::topology::NODE_RANK);
  stk::mesh::Part & side_2 = meta_data.declare_part("side_2", stk::topology::NODE_RANK);

  stk::mesh::Part & side_1 = meta_data.declare_part("side_1", stk::topology::NODE_RANK);
  stk::mesh::Part & side_3 = meta_data.declare_part("side_3", stk::topology::NODE_RANK);

  stk::mesh::Field<double> & volField = meta_data.declare_field<double>(stk::topology::NODE_RANK, "volume");
  stk::mesh::put_field_on_mesh(volField, meta_data.universal_part(), nullptr);

  meta_data.commit();

  fixture.generate_mesh();

  stk::mesh::PartVector side_0_parts(1,&side_0);
  stk::mesh::PartVector side_2_parts(1,&side_2);

  stk::mesh::PartVector side_1_parts(1,&side_1);
  stk::mesh::PartVector side_3_parts(1,&side_3);

  assign_to_parts_for_two_way(x, y, z, fixture, bulk_data,
                                side_0_parts, side_1_parts, side_2_parts, side_3_parts);

  //do periodic search
  typedef stk::mesh::GetCoordinates<CoordFieldType> CoordinateFunctor;
  typedef stk::mesh::PeriodicBoundarySearch<CoordinateFunctor> PeriodicSearch;
  PeriodicSearch pbc_search(bulk_data, CoordinateFunctor(bulk_data, coords_field));

  const stk::mesh::Selector side_0_selector = side_0 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());
  const stk::mesh::Selector side_1_selector = side_1 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());
  const stk::mesh::Selector side_2_selector = side_2 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());
  const stk::mesh::Selector side_3_selector = side_3 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());

  pbc_search.add_linear_periodic_pair(side_0_selector, side_2_selector );
  pbc_search.add_linear_periodic_pair(side_1_selector, side_3_selector );
  pbc_search.find_periodic_nodes(bulk_data.parallel());

  //check to see if re-entrant
  pbc_search.find_periodic_nodes(bulk_data.parallel());

  check_gold_two_way_multiperiodic(pbc_search.get_pairs());

  //now we ghost everything to do a local search
  bulk_data.modification_begin();
  pbc_search.create_ghosting("periodic_ghosts");
  bulk_data.modification_end();

  do_volume_assembly(bulk_data, volField);

  std::vector< stk::mesh::FieldBase const * > ghosted_fields;
  ghosted_fields.push_back(&volField);
  stk::mesh::communicate_field_data( pbc_search.get_ghosting(), ghosted_fields);

  do_periodic_assembly(bulk_data, pbc_search, volField);

  //interior of domain should be 1.0
  for (unsigned i=0; i<x+1u; ++i) {
    for (unsigned j=0; j<y+1u; ++j) {
      for (unsigned k=1; k<z; ++k) {
        stk::mesh::Entity node = fixture.node(i,j,k);
        expect_eq_for_shared_or_owned_node(bulk_data, node, volField, 1.0);
      }
    }
  }
  //faces (not edges) should be 0.5
  //there are no corners or edges since this is two way periodic
  for (unsigned i=0; i<x+1u; ++i) {
    //front and back
    for (unsigned j=0; j<y+1u; ++j) {
      const unsigned kFront = 0;
      const unsigned kBack = 3;
      stk::mesh::Entity node = fixture.node(i,j,kFront);
      expect_eq_for_shared_or_owned_node(bulk_data, node, volField, 0.5);
      node = fixture.node(i,j,kBack);
      expect_eq_for_shared_or_owned_node(bulk_data, node, volField, 0.5);
    }
  }
}

void assign_to_parts_for_three_way(const unsigned x, const unsigned y, const unsigned z,
                                   stk::mesh::fixtures::HexFixture &fixture,
                                   stk::mesh::BulkData &bulk_data,
                                   stk::mesh::PartVector &side_0_parts,
                                   stk::mesh::PartVector &side_1_parts,
                                   stk::mesh::PartVector &side_2_parts,
                                   stk::mesh::PartVector &side_3_parts,
                                   stk::mesh::PartVector &side_4_parts,
                                   stk::mesh::PartVector &side_5_parts)
{
  assign_to_parts_for_two_way(x, y, z, fixture, bulk_data,
                              side_0_parts, side_1_parts, side_2_parts, side_3_parts);

  //add periodic pair for side 4 and side 5
  bulk_data.modification_begin();
  for (unsigned i=0; i<y+1u; ++i) {
  for (unsigned j=0; j<z+1u; ++j) {
    stk::mesh::Entity node_4 = fixture.node(i,j,0);
    if (bulk_data.is_valid(node_4)  && bulk_data.bucket(node_4).owned()) {
      bulk_data.change_entity_parts( fixture.node(i,j,0), side_4_parts);
    }
    stk::mesh::Entity node_5 = fixture.node(i,j,3);
    if (bulk_data.is_valid(node_5)  && bulk_data.bucket(node_5).owned()) {
      bulk_data.change_entity_parts( fixture.node(i,j,3), side_5_parts);
    }
  }}
  bulk_data.modification_end();
}

TEST(CoarseSearch, ThreeWayMultiPeriodicBC)
{
  const unsigned x = 3, y = 3, z = 3;

  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, x, y, z);

  stk::mesh::BulkData & bulk_data = fixture.m_bulk_data;
  stk::mesh::MetaData & meta_data = fixture.m_meta;
  CoordFieldType & coords_field = *fixture.m_coord_field;

  stk::mesh::Part & side_0 = meta_data.declare_part("side_0", stk::topology::NODE_RANK);
  stk::mesh::Part & side_2 = meta_data.declare_part("side_2", stk::topology::NODE_RANK);

  stk::mesh::Part & side_1 = meta_data.declare_part("side_1", stk::topology::NODE_RANK);
  stk::mesh::Part & side_3 = meta_data.declare_part("side_3", stk::topology::NODE_RANK);

  stk::mesh::Part & side_4 = meta_data.declare_part("side_4", stk::topology::NODE_RANK);
  stk::mesh::Part & side_5 = meta_data.declare_part("side_5", stk::topology::NODE_RANK);

  stk::mesh::Field<double> & volField = meta_data.declare_field<double>(stk::topology::NODE_RANK, "volume");
  stk::mesh::put_field_on_mesh(volField, meta_data.universal_part(), nullptr);

  meta_data.commit();

  fixture.generate_mesh();

  stk::mesh::PartVector side_0_parts(1, &side_0);
  stk::mesh::PartVector side_2_parts(1, &side_2);

  stk::mesh::PartVector side_1_parts(1, &side_1);
  stk::mesh::PartVector side_3_parts(1, &side_3);

  stk::mesh::PartVector side_4_parts(1, &side_4);
  stk::mesh::PartVector side_5_parts(1, &side_5);

  assign_to_parts_for_three_way(x, y, z, fixture, bulk_data,
                                  side_0_parts, side_1_parts,
                                  side_2_parts, side_3_parts,
                                  side_4_parts, side_5_parts);

  //do periodic search
  typedef stk::mesh::GetCoordinates<CoordFieldType> CoordinateFunctor;
  typedef stk::mesh::PeriodicBoundarySearch<CoordinateFunctor> PeriodicSearch;
  PeriodicSearch pbc_search(bulk_data, CoordinateFunctor(bulk_data, coords_field));

  const stk::mesh::Selector side_0_selector = side_0 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());
  const stk::mesh::Selector side_1_selector = side_1 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());
  const stk::mesh::Selector side_2_selector = side_2 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());
  const stk::mesh::Selector side_3_selector = side_3 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());
  const stk::mesh::Selector side_4_selector = side_4 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());
  const stk::mesh::Selector side_5_selector = side_5 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());

  pbc_search.add_linear_periodic_pair(side_0_selector, side_2_selector);
  pbc_search.add_linear_periodic_pair(side_1_selector, side_3_selector);
  pbc_search.add_linear_periodic_pair(side_4_selector, side_5_selector);
  pbc_search.find_periodic_nodes(bulk_data.parallel());

  check_gold_three_way_multiperiodic(pbc_search.get_pairs());

  //now we ghost everything to do a local search
  bulk_data.modification_begin();
  pbc_search.create_ghosting("periodic_ghosts");
  bulk_data.modification_end();

  do_volume_assembly(bulk_data, volField);

  std::vector< stk::mesh::FieldBase const * > ghosted_fields;
  ghosted_fields.push_back(&volField);
  stk::mesh::communicate_field_data( pbc_search.get_ghosting(), ghosted_fields);

  do_periodic_assembly(bulk_data, pbc_search, volField);

  //interior of domain should be 1.0
  for (unsigned i=0; i<x+1u; ++i) {
    for (unsigned j=0; j<y+1u; ++j) {
      for (unsigned k=0; k<z+1u; ++k) {
        stk::mesh::Entity node = fixture.node(i,j,k);
        expect_eq_for_shared_or_owned_node(bulk_data, node, volField, 1.0);
      }
    }
  }
}



TEST(CoarseSearch, MultiPeriodicBCDisallowRotational)
{
  const unsigned x = 3, y = 3, z = 3;

  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, x, y, z);

  stk::mesh::BulkData & bulk_data = fixture.m_bulk_data;
  stk::mesh::MetaData & meta_data = fixture.m_meta;
  CoordFieldType & coords_field = *fixture.m_coord_field;

  stk::mesh::Part & side_0 = meta_data.declare_part("side_0", stk::topology::NODE_RANK);
  stk::mesh::Part & side_2 = meta_data.declare_part("side_2", stk::topology::NODE_RANK);

  stk::mesh::Part & side_1 = meta_data.declare_part("side_1", stk::topology::NODE_RANK);
  stk::mesh::Part & side_3 = meta_data.declare_part("side_3", stk::topology::NODE_RANK);

  stk::mesh::Part & side_4 = meta_data.declare_part("side_4", stk::topology::NODE_RANK);
  stk::mesh::Part & side_5 = meta_data.declare_part("side_5", stk::topology::NODE_RANK);

  meta_data.commit();

  fixture.generate_mesh();

  stk::mesh::PartVector side_0_parts(1, &side_0);
  stk::mesh::PartVector side_2_parts(1, &side_2);

  stk::mesh::PartVector side_1_parts(1, &side_1);
  stk::mesh::PartVector side_3_parts(1, &side_3);

  stk::mesh::PartVector side_4_parts(1, &side_4);
  stk::mesh::PartVector side_5_parts(1, &side_5);

  // Any assignment is okay, since we don't care about the search results.
  assign_to_parts_for_three_way(x, y, z, fixture, bulk_data,
                                  side_0_parts, side_1_parts,
                                  side_2_parts, side_3_parts,
                                  side_4_parts, side_5_parts);

  const double rotationAngle = TWO_PI/4.0;
  const double rotationAxis[3] = {0.0, 0.0, 1.0};
  const double axisLocation[3] = {0.0, 0.0, 0.0};

  typedef stk::mesh::GetCoordinates<CoordFieldType> CoordinateFunctor;
  typedef stk::mesh::PeriodicBoundarySearch<CoordinateFunctor> PeriodicSearch;

  PeriodicSearch pbc_search_caseA(bulk_data, CoordinateFunctor(bulk_data, coords_field));
  pbc_search_caseA.add_rotational_periodic_pair(side_0 & meta_data.locally_owned_part(),
                                                side_2 & meta_data.locally_owned_part(),
                                                rotationAngle, rotationAxis, axisLocation);

  EXPECT_THROW( pbc_search_caseA.add_rotational_periodic_pair(side_1 & meta_data.locally_owned_part(),
                                                side_3 & meta_data.locally_owned_part(),
                                                rotationAngle, rotationAxis, axisLocation),
                                                std::exception);
}


TEST(CoarseSearch, RotationalPeriodicBC)
{
  const unsigned x = 3, y = 3, z = 3;

  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, x, y, z);

  stk::mesh::BulkData & bulk_data = fixture.m_bulk_data;
  stk::mesh::MetaData & meta_data = fixture.m_meta;
  CoordFieldType & coords_field = *fixture.m_coord_field;

  stk::mesh::Part & side_0 = meta_data.declare_part("side_0", stk::topology::NODE_RANK);
  stk::mesh::Part & side_3 = meta_data.declare_part("side_3", stk::topology::NODE_RANK);

  stk::mesh::Field<double> & volField = meta_data.declare_field<double>(stk::topology::NODE_RANK, "volume");
  stk::mesh::put_field_on_mesh(volField, meta_data.universal_part(), nullptr);

  meta_data.commit();

  const double rotationAngle = -TWO_PI/4.0;
  const double rotationAxis[3] = {0.0, 0.0, 1.0};
  const double axisLocation[3] = {0.0, 0.0, 0.0};
  stk::mesh::fixtures::CylindricalCoordinateMapping coordMap(1.0, rotationAngle, 4);
  fixture.generate_mesh(coordMap);

  stk::mesh::PartVector independent_parts(1,&side_0);
  stk::mesh::PartVector dependent_parts(1,&side_3);

  bulk_data.modification_begin();
  for (unsigned i=0; i<y+1u; ++i) {
  for (unsigned j=0; j<z+1u; ++j) {
    stk::mesh::Entity node_0 = fixture.node(0,i,j);
    if (bulk_data.is_valid(node_0)  && bulk_data.bucket(node_0).owned()) {
      bulk_data.change_entity_parts( fixture.node(0,i,j), independent_parts);
    }
    stk::mesh::Entity node_3 = fixture.node(3,i,j);
    if (bulk_data.is_valid(node_3)  && bulk_data.bucket(node_3).owned()) {
      bulk_data.change_entity_parts( fixture.node(3,i,j), dependent_parts);
    }
  }}
  bulk_data.modification_end();

  //do periodic search
  typedef stk::mesh::GetCoordinates<CoordFieldType> CoordinateFunctor;
  typedef stk::mesh::PeriodicBoundarySearch<CoordinateFunctor> PeriodicSearch;
  PeriodicSearch pbc_search(bulk_data, CoordinateFunctor(bulk_data, coords_field));

  const stk::mesh::Selector side_0_selector = side_0 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());
  const stk::mesh::Selector side_3_selector = side_3 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());


  pbc_search.add_rotational_periodic_pair(side_0_selector,
                                          side_3_selector,
                                          rotationAngle,
                                          rotationAxis,
                                          axisLocation);

  pbc_search.find_periodic_nodes(bulk_data.parallel());


  check_gold(pbc_search.get_pairs() );

  bulk_data.modification_begin();
  pbc_search.create_ghosting("periodic_ghosts");
  bulk_data.modification_end();

  do_volume_assembly(bulk_data, volField);

  std::vector< stk::mesh::FieldBase const * > ghosted_fields;
  ghosted_fields.push_back(&volField);
  stk::mesh::communicate_field_data( pbc_search.get_ghosting(), ghosted_fields);

  do_periodic_assembly(bulk_data, pbc_search, volField);

  check_single_periodic_assembly(bulk_data, fixture, volField, x, y, z);

  check_rotation_matrix(pbc_search, rotationAngle);
}

TEST(CoarseSearch, OffsetRotationalPeriodicBC)
{
  const unsigned x = 3, y = 3, z = 3;

  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, x, y, z);

  stk::mesh::BulkData & bulk_data = fixture.m_bulk_data;
  stk::mesh::MetaData & meta_data = fixture.m_meta;
  CoordFieldType & coords_field = *fixture.m_coord_field;

  stk::mesh::Part & side_0 = meta_data.declare_part("side_0", stk::topology::NODE_RANK);
  stk::mesh::Part & side_3 = meta_data.declare_part("side_3", stk::topology::NODE_RANK);

  stk::mesh::Field<double> & volField = meta_data.declare_field<double>(stk::topology::NODE_RANK, "volume");
  stk::mesh::put_field_on_mesh(volField, meta_data.universal_part(), nullptr);

  meta_data.commit();

  const double rotationAngle = -TWO_PI/4.0;
  const double rotationAxis[3] = {0.0, 0.0, 1.0};
  const double axisLocation[3] = {1.0, 0.0, 0.0};
  CylindricalCoordinateMappingWithOffset coordMap(1.0, rotationAngle, 4, 1.0, 0.0, 0.0);
  fixture.generate_mesh(coordMap);

  stk::mesh::PartVector independent_parts(1,&side_0);
  stk::mesh::PartVector dependent_parts(1,&side_3);

  bulk_data.modification_begin();
  for (unsigned i=0; i<y+1u; ++i) {
  for (unsigned j=0; j<z+1u; ++j) {
    stk::mesh::Entity node_0 = fixture.node(0,i,j);
    if (bulk_data.is_valid(node_0)  && bulk_data.bucket(node_0).owned()) {
      bulk_data.change_entity_parts( fixture.node(0,i,j), independent_parts);
    }
    stk::mesh::Entity node_3 = fixture.node(3,i,j);
    if (bulk_data.is_valid(node_3)  && bulk_data.bucket(node_3).owned()) {
      bulk_data.change_entity_parts( fixture.node(3,i,j), dependent_parts);
    }
  }}
  bulk_data.modification_end();

  //do periodic search
  typedef stk::mesh::GetCoordinates<CoordFieldType> CoordinateFunctor;
  typedef stk::mesh::PeriodicBoundarySearch<CoordinateFunctor> PeriodicSearch;
  PeriodicSearch pbc_search(bulk_data, CoordinateFunctor(bulk_data, coords_field));

  const stk::mesh::Selector side_0_selector = side_0 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());
  const stk::mesh::Selector side_3_selector = side_3 & (meta_data.locally_owned_part() | meta_data.globally_shared_part());


  pbc_search.add_rotational_periodic_pair(side_0_selector,
                                          side_3_selector,
                                          rotationAngle,


                                          rotationAxis,
                                          axisLocation);

  pbc_search.find_periodic_nodes(bulk_data.parallel());


  check_gold(pbc_search.get_pairs() );

  bulk_data.modification_begin();
  pbc_search.create_ghosting("periodic_ghosts");
  bulk_data.modification_end();

  do_volume_assembly(bulk_data, volField);

  std::vector< stk::mesh::FieldBase const * > ghosted_fields;
  ghosted_fields.push_back(&volField);
  stk::mesh::communicate_field_data( pbc_search.get_ghosting(), ghosted_fields);

  do_periodic_assembly(bulk_data, pbc_search, volField);

  check_single_periodic_assembly(bulk_data, fixture, volField, x, y, z);

  check_rotation_matrix(pbc_search, rotationAngle);
}

TEST(PeriodicBoundarySearch, testRotationMatrixForRotationAboutZ)
{
    double angleInRadians = PI/2;
    double axis[3] = { 0, 0, 1 };

    stk::mesh::matrix3x3 matrix1;
    stk::mesh::fillRotationMatrix(angleInRadians, axis[0], axis[1], axis[2], matrix1);

    stk::mesh::matrix3x3 goldMatrix;

    goldMatrix.setData(0, 0);
    goldMatrix.setData(1, -1);
    goldMatrix.setData(2, 0);

    goldMatrix.setData(3, 1);
    goldMatrix.setData(4, 0);
    goldMatrix.setData(5, 0);

    goldMatrix.setData(6, 0);
    goldMatrix.setData(7, 0);
    goldMatrix.setData(8, 1);

    EXPECT_EQ(goldMatrix.getData(0), 0);
    EXPECT_EQ(goldMatrix.getData(1), -1);
    EXPECT_EQ(goldMatrix.getData(2), 0);
    EXPECT_EQ(goldMatrix.getData(3), 1);
    EXPECT_EQ(goldMatrix.getData(4), 0);
    EXPECT_EQ(goldMatrix.getData(5), 0);
    EXPECT_EQ(goldMatrix.getData(6), 0);
    EXPECT_EQ(goldMatrix.getData(7), 0);
    EXPECT_EQ(goldMatrix.getData(8), 1);

    EXPECT_EQ(goldMatrix.getData(0,0), 0);
    EXPECT_EQ(goldMatrix.getData(0,1), -1);
    EXPECT_EQ(goldMatrix.getData(0,2), 0);
    EXPECT_EQ(goldMatrix.getData(1,0), 1);
    EXPECT_EQ(goldMatrix.getData(1,1), 0);
    EXPECT_EQ(goldMatrix.getData(1,2), 0);
    EXPECT_EQ(goldMatrix.getData(2,0), 0);
    EXPECT_EQ(goldMatrix.getData(2,1), 0);
    EXPECT_EQ(goldMatrix.getData(2,2), 1);

    double tol = 1e-6;

    EXPECT_EQ(9, goldMatrix.numEntries());

    for (int i=0;i<goldMatrix.numEntries();i++)
    {
        EXPECT_NEAR(goldMatrix.getData(i), matrix1.getData(i), tol) << " for index = " << i;
    }

    for (int row=0;row<3;row++)
    {
        for (int col=0;col<3;col++)
        {
            EXPECT_NEAR(goldMatrix.getData(row,col), matrix1.getData(row,col), tol) << " for (row,col) = (" << row << "," << col << ")";
        }
    }

    double coordinate[3] = { 1, 0, 0 };
    double result[3] = { 0, 0, 0 };
    double goldResult[3] = { 0, 1, 0 };

    matrix1.transformVec(coordinate, result);

    for (int i=0;i<3;i++)
    {
        EXPECT_NEAR(goldResult[i], result[i], tol) << " for index = " << i;
    }
}


TEST(PeriodicBoundarySearch, testRotationMatrixForRotationAboutY)
{
    double angleInRadians = PI/6;
    double axis[3] = { 0, 1, 0 };

    stk::mesh::matrix3x3 matrix1;
    stk::mesh::fillRotationMatrix(angleInRadians, axis[0], axis[1], axis[2], matrix1);

    stk::mesh::matrix3x3 goldMatrix;


    goldMatrix.setData(0, cos(angleInRadians));
    goldMatrix.setData(1, 0);
    goldMatrix.setData(2, sin(angleInRadians));

    goldMatrix.setData(3, 0);
    goldMatrix.setData(4, 1);
    goldMatrix.setData(5, 0);

    goldMatrix.setData(6, -sin(angleInRadians));
    goldMatrix.setData(7, 0);
    goldMatrix.setData(8, cos(angleInRadians));

    double tol = 1e-6;

    EXPECT_EQ(9, goldMatrix.numEntries());

    for (int i=0;i<goldMatrix.numEntries();i++)
    {
        EXPECT_NEAR(goldMatrix.getData(i), matrix1.getData(i), tol) << " for index = " << i;
    }
}

TEST(PeriodicBoundarySearch, testRotationMatrixForRotationAboutX)
{
    double angleInRadians = PI/3;

    double axis[3] = { 2, 0, 0 };

    stk::mesh::matrix3x3 matrix1;
    stk::mesh::fillRotationMatrix(angleInRadians, axis[0], axis[1], axis[2], matrix1);

    stk::mesh::matrix3x3 goldMatrix;


    goldMatrix.setData(0, 1);
    goldMatrix.setData(1, 0);
    goldMatrix.setData(2, 0);

    goldMatrix.setData(3, 0);
    goldMatrix.setData(4, cos(angleInRadians));
    goldMatrix.setData(5, -sin(angleInRadians));

    goldMatrix.setData(6, 0);
    goldMatrix.setData(7, sin(angleInRadians));
    goldMatrix.setData(8, cos(angleInRadians));

    double tol = 1e-6;

    EXPECT_EQ(9, goldMatrix.numEntries());

    for (int i=0;i<goldMatrix.numEntries();i++)
    {
        EXPECT_NEAR(goldMatrix.getData(i), matrix1.getData(i), tol) << " for index = " << i;
    }
}


