#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_mesh/fixtures/HexFixture.hpp>
#include <stk_search_util/stk_mesh/PeriodicBoundarySearch.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stk_mesh/fixtures/GearsFixture.hpp>
#include <stk_mesh/fixtures/Gear.hpp>
#include <stk_mesh/base/Field.hpp>


typedef stk::mesh::fixtures::HexFixture::CoordFieldType CoordFieldType;
typedef stk::mesh::GetCoordinates<CoordFieldType> CoordinateFunctor;
typedef stk::mesh::PeriodicBoundarySearch<CoordinateFunctor> PeriodicSearch;

namespace {

void expect_eq_for_shared_or_owned_node(stk::mesh::BulkData & bulk_data, stk::mesh::Entity node, const stk::mesh::Field<double> & theField, double expected_value )
{
  if (bulk_data.is_valid(node) && (bulk_data.bucket(node).owned() || bulk_data.bucket(node).shared() ) )
  {
    const double * const vol = bulk_data.field_data(theField, node);
    EXPECT_EQ(*vol, expected_value);
  }
}

void do_periodic_assembly(stk::mesh::BulkData & bulk_data, PeriodicSearch & pbc_search, stk::mesh::Field<double> & volField)
{
  //gather to domain node from possibly multiple ranges
  for (size_t i = 0; i < pbc_search.size(); ++i)
  {
    std::pair<stk::mesh::Entity, stk::mesh::Entity> entityPair = pbc_search.get_node_pair(i);
    ThrowRequire(bulk_data.is_valid(entityPair.first));
    ThrowRequire(bulk_data.is_valid(entityPair.second));
    double * domainField = bulk_data.field_data(volField, entityPair.first);
    double * rangeField = bulk_data.field_data(volField, entityPair.second);
    *domainField += *rangeField;
  }

  std::vector< stk::mesh::FieldBase const * > ghosted_fields;
  ghosted_fields.push_back(&volField);
  stk::mesh::communicate_field_data( pbc_search.get_ghosting(), ghosted_fields);

  //set ranges equal to domain value
  for (size_t i = 0; i < pbc_search.size(); ++i)
  {
    std::pair<stk::mesh::Entity, stk::mesh::Entity> entityPair = pbc_search.get_node_pair(i);
    ThrowRequire(bulk_data.is_valid(entityPair.first));
    ThrowRequire(bulk_data.is_valid(entityPair.second));
    double * domainField = bulk_data.field_data(volField, entityPair.first);
    double * rangeField = bulk_data.field_data(volField, entityPair.second);
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
        double * vol = bulk_data.field_data(volField, aNode[in]);
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
    stk::mesh::EntityId domain_node = search_results[i].first.ident.id();
    stk::mesh::EntityId range_node = search_results[i].second.ident.id();

    //entry in search is found in gold

    bool found = std::find(gold.begin(),gold.end(),std::make_pair(domain_node,range_node)) != gold.end();
    if (!found)
    {
      std::cout << "We can't find domain/range: " << search_results[i].first << "/" << search_results[i].second  << std::endl;
      EXPECT_TRUE(0);
    }
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
    stk::mesh::EntityId domain_node = search_results[i].first.ident.id();
    stk::mesh::EntityId range_node = search_results[i].second.ident.id();

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
    stk::mesh::EntityId domain_node = search_results[i].first.ident.id();
    stk::mesh::EntityId range_node = search_results[i].second.ident.id();

    EXPECT_TRUE((std::find(gold.begin(), gold.end(), std::make_pair(domain_node,range_node) ) ) != gold.end());
  }
}

void print_periodic_node_pairs(stk::mesh::BulkData & bulk_data,
    PeriodicSearch & pbc_search,
    CoordFieldType & coords_field)
{
  std::cout << "Periodic nodes identified:" << std::endl;
  for (size_t i = 0; i < pbc_search.size(); ++i)
  {
    std::pair<stk::mesh::Entity, stk::mesh::Entity> aPair = pbc_search.get_node_pair(i);
    stk::mesh::EntityId domainId = bulk_data.identifier(aPair.first);
    stk::mesh::EntityId rangeId = bulk_data.identifier(aPair.second);
    std::cout << "My Proc: " << std::setw(6) << bulk_data.parallel_rank() << std::setw(12) << domainId << ":";
    std::cout << std::setw(12) << rangeId << std::endl;
  }
}
}// namespace

STKUNIT_UNIT_TEST(CoarseSearch, PeriodicBC)
{
  const unsigned x = 3, y = 3, z = 3;


  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, x, y, z);

  stk::mesh::BulkData & bulk_data = fixture.m_bulk_data;
  stk::mesh::MetaData & meta_data = fixture.m_meta;
  CoordFieldType & coords_field = fixture.m_coord_field;

  stk::mesh::Part & side_0 = meta_data.declare_part("side_0", stk::topology::NODE_RANK);
  stk::mesh::Part & side_3 = meta_data.declare_part("side_3", stk::topology::NODE_RANK);

  stk::mesh::Field<double> & volField = meta_data.declare_field<stk::mesh::Field<double> >("volume");
  stk::mesh::put_field(volField, stk::topology::NODE_RANK, meta_data.universal_part());

  meta_data.commit();

  fixture.generate_mesh();

  // side 0 (master) is periodic with side 3 (slave)

  // add nodes to side 0 and 3
  stk::mesh::PartVector side_0_parts(1,&side_0);
  stk::mesh::PartVector side_3_parts(1,&side_3);

  bulk_data.modification_begin();
  for (unsigned i=0; i<y+1u; ++i) {
  for (unsigned j=0; j<z+1u; ++j) {
    stk::mesh::Entity node_0 = fixture.node(0,i,j);
    if (bulk_data.is_valid(node_0)  && bulk_data.bucket(node_0).owned()) {
      bulk_data.change_entity_parts( fixture.node(0,i,j), side_0_parts);
    }
    stk::mesh::Entity node_3 = fixture.node(3,i,j);
    if (bulk_data.is_valid(node_3)  && bulk_data.bucket(node_3).owned()) {
      bulk_data.change_entity_parts( fixture.node(3,i,j), side_3_parts);
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
  pbc_search.find_periodic_nodes(bulk_data.parallel(), true);


  check_gold(pbc_search.get_pairs() );

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

STKUNIT_UNIT_TEST(CoarseSearch, TwoWayMultiPeriodicBC)
{
  const unsigned x = 3, y = 3, z = 3;


  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, x, y, z);

  stk::mesh::BulkData & bulk_data = fixture.m_bulk_data;
  stk::mesh::MetaData & meta_data = fixture.m_meta;
  CoordFieldType & coords_field = fixture.m_coord_field;

  stk::mesh::Part & side_0 = meta_data.declare_part("side_0", stk::topology::NODE_RANK);
  stk::mesh::Part & side_2 = meta_data.declare_part("side_2", stk::topology::NODE_RANK);

  stk::mesh::Part & side_1 = meta_data.declare_part("side_1", stk::topology::NODE_RANK);
  stk::mesh::Part & side_3 = meta_data.declare_part("side_3", stk::topology::NODE_RANK);

  stk::mesh::Field<double> & volField = meta_data.declare_field<stk::mesh::Field<double> >("volume");
  stk::mesh::put_field(volField, stk::topology::NODE_RANK, meta_data.universal_part());

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
  pbc_search.find_periodic_nodes(bulk_data.parallel(), true);

  //check to see if re-entrant
//  pbc_search.find_periodic_nodes(bulk_data.parallel());

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

STKUNIT_UNIT_TEST(CoarseSearch, ThreeWayMultiPeriodicBC)
{
  const unsigned x = 3, y = 3, z = 3;

  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, x, y, z);

  stk::mesh::BulkData & bulk_data = fixture.m_bulk_data;
  stk::mesh::MetaData & meta_data = fixture.m_meta;
  CoordFieldType & coords_field = fixture.m_coord_field;

  stk::mesh::Part & side_0 = meta_data.declare_part("side_0", stk::topology::NODE_RANK);
  stk::mesh::Part & side_2 = meta_data.declare_part("side_2", stk::topology::NODE_RANK);

  stk::mesh::Part & side_1 = meta_data.declare_part("side_1", stk::topology::NODE_RANK);
  stk::mesh::Part & side_3 = meta_data.declare_part("side_3", stk::topology::NODE_RANK);

  stk::mesh::Part & side_4 = meta_data.declare_part("side_4", stk::topology::NODE_RANK);
  stk::mesh::Part & side_5 = meta_data.declare_part("side_5", stk::topology::NODE_RANK);

  stk::mesh::Field<double> & volField = meta_data.declare_field<stk::mesh::Field<double> >("volume");
  stk::mesh::put_field(volField, stk::topology::NODE_RANK, meta_data.universal_part());

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
  pbc_search.find_periodic_nodes(bulk_data.parallel(), true);

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



STKUNIT_UNIT_TEST(CoarseSearch, MultiPeriodicBCDisallowRotational)
{
  const unsigned x = 3, y = 3, z = 3;

  stk::mesh::fixtures::HexFixture fixture(MPI_COMM_WORLD, x, y, z);

  stk::mesh::BulkData & bulk_data = fixture.m_bulk_data;
  stk::mesh::MetaData & meta_data = fixture.m_meta;
  CoordFieldType & coords_field = fixture.m_coord_field;

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

  const double rotationAngle = TWO_PI/2.0;
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


template<typename SearchPairVector>
void check_gold_rotational_multiperiodic( const SearchPairVector & search_results )
  // check search result
{
  typedef std::vector<std::pair<stk::mesh::EntityId,stk::mesh::EntityId> > GoldVector;
  GoldVector gold;
  gold.push_back(std::make_pair(1, 181));
  gold.push_back(std::make_pair(2, 182));
  gold.push_back(std::make_pair(3, 183));
  gold.push_back(std::make_pair(4, 184));
  gold.push_back(std::make_pair(5, 185));
  gold.push_back(std::make_pair(6, 186));
  gold.push_back(std::make_pair(7, 187));
  gold.push_back(std::make_pair(8, 188));
  gold.push_back(std::make_pair(9, 189));
  gold.push_back(std::make_pair(11, 191));
  gold.push_back(std::make_pair(12, 192));
  gold.push_back(std::make_pair(13, 193));
  gold.push_back(std::make_pair(14, 194));
  gold.push_back(std::make_pair(15, 195));
  gold.push_back(std::make_pair(16, 196));
  gold.push_back(std::make_pair(17, 197));
  gold.push_back(std::make_pair(18, 198));
  gold.push_back(std::make_pair(19, 199));
  gold.push_back(std::make_pair(21, 201));
  gold.push_back(std::make_pair(22, 202));
  gold.push_back(std::make_pair(23, 203));
  gold.push_back(std::make_pair(24, 204));
  gold.push_back(std::make_pair(25, 205));
  gold.push_back(std::make_pair(26, 206));
  gold.push_back(std::make_pair(27, 207));
  gold.push_back(std::make_pair(28, 208));
  gold.push_back(std::make_pair(29, 209));
  gold.push_back(std::make_pair(31, 211));
  gold.push_back(std::make_pair(32, 212));
  gold.push_back(std::make_pair(33, 213));
  gold.push_back(std::make_pair(34, 214));
  gold.push_back(std::make_pair(35, 215));
  gold.push_back(std::make_pair(36, 216));
  gold.push_back(std::make_pair(37, 217));
  gold.push_back(std::make_pair(38, 218));
  gold.push_back(std::make_pair(39, 219));

  for (size_t i=0, size=search_results.size(); i<size; ++i) {
    stk::mesh::EntityId domain_node = search_results[i].first.ident.id();
    stk::mesh::EntityId range_node = search_results[i].second.ident.id();

    EXPECT_TRUE((std::find(gold.begin(), gold.end(), std::make_pair(domain_node,range_node) ) ) != gold.end());
  }
}

STKUNIT_UNIT_TEST(CoarseSearch, RotationalPeriodicBC)
{
  stk::mesh::fixtures::GearsFixture fixture(MPI_COMM_SELF, 1);
  stk::mesh::BulkData & bulk_data = fixture.bulk_data;
  stk::mesh::MetaData & meta_data = fixture.meta_data;
  CoordFieldType & coords_field = fixture.cartesian_coord_field;

  stk::mesh::Part & side_0 = meta_data.declare_part("side_0", stk::topology::NODE_RANK);
  stk::mesh::Part & side_1 = meta_data.declare_part("side_1", stk::topology::NODE_RANK);
  stk::mesh::PartVector side_0_vector(1, & side_0);
  stk::mesh::PartVector side_1_vector(1, & side_1);

  meta_data.commit();
  fixture.generate_mesh();

  //loop mesh entities and flag a wedge domain for putting in another part
  const stk::mesh::fixtures::Gear & aGear = fixture.get_gear(0);

  typedef std::pair<stk::mesh::Entity,stk::mesh::Entity> NodePairType;
  std::vector<NodePairType> nodePairs;

  bulk_data.modification_begin();
  //side_0 at angle 0
  const size_t domainAngleIndex = 0;
  const size_t rangeAngleIndex = 3;

  for ( size_t ir = 0 ; ir < aGear.rad_num -2 ; ++ir ) {
    for ( size_t iz = 0 ; iz < aGear.height_num -1 ; ++iz ) {
      stk::mesh::Entity domainNode = aGear.get_node(iz, ir, domainAngleIndex);
      if ( bulk_data.is_valid(domainNode) && bulk_data.bucket(domainNode).owned() )
      {
        bulk_data.change_entity_parts(domainNode, side_0_vector);
      }
      stk::mesh::Entity rangeNode = aGear.get_node(iz, ir, rangeAngleIndex);
      if ( bulk_data.is_valid(rangeNode) && bulk_data.bucket(rangeNode).owned() )
      {
        bulk_data.change_entity_parts(rangeNode, side_1_vector);
      }
      nodePairs.push_back(std::make_pair(domainNode, rangeNode));
    }
  }
  bulk_data.modification_end();

  //do periodic search
  typedef stk::mesh::GetCoordinates<CoordFieldType> CoordinateFunctor;
  typedef stk::mesh::PeriodicBoundarySearch<CoordinateFunctor> PeriodicSearch;
  PeriodicSearch pbc_search(bulk_data, CoordinateFunctor(bulk_data, coords_field));

  //element_size is in radians
  const double rotationAngle = double(rangeAngleIndex - domainAngleIndex)*TWO_PI/aGear.angle_num;
  const double rotationAxis[3] = {0.0, 0.0, 1.0};
  const double axisLocation[3] = {0.0, 0.0, 0.0};

  pbc_search.add_rotational_periodic_pair(side_0 & meta_data.locally_owned_part(),
                                          side_1 & meta_data.locally_owned_part(),
                                          rotationAngle,
                                          rotationAxis,
                                          axisLocation);

  pbc_search.find_periodic_nodes(MPI_COMM_SELF);

  check_gold_rotational_multiperiodic(pbc_search.get_pairs());
  if (bulk_data.parallel_size() == 1)
  {
    EXPECT_EQ(pbc_search.get_pairs().size(), 36u);
  }
  else
  {
    const int local_search_count = pbc_search.get_pairs().size();
    int global_search_count=0;
    stk::all_reduce_sum(bulk_data.parallel(), &local_search_count, &global_search_count, 1);

    EXPECT_GE(global_search_count, 36);
  }

  //now we ghost everything to do a local search
  bulk_data.modification_begin();
  pbc_search.create_ghosting("periodic_ghosts");
  bulk_data.modification_end();

  const stk::mesh::Ghosting & periodic_bc_ghosting = pbc_search.get_ghosting();

  std::vector< stk::mesh::FieldBase const * > ghosted_fields;
  ghosted_fields.push_back(&coords_field);
  stk::mesh::communicate_field_data( periodic_bc_ghosting, ghosted_fields);

  //lets do a local search to make sure the coords field and entities were properly ghosted
  PeriodicSearch pbc_local_search(bulk_data, CoordinateFunctor(bulk_data, coords_field));
  pbc_local_search.add_rotational_periodic_pair(side_0 & meta_data.locally_owned_part(),
                                          side_1 & meta_data.locally_owned_part(),
                                          rotationAngle,
                                          rotationAxis,
                                          axisLocation);

  pbc_local_search.find_periodic_nodes(MPI_COMM_SELF, true);

  //test to make sure it is re-entrant
  pbc_local_search.find_periodic_nodes(MPI_COMM_SELF, true);

  check_gold_rotational_multiperiodic(pbc_local_search.get_pairs());
  if (bulk_data.parallel_size() == 1)
  {
    EXPECT_EQ(pbc_local_search.get_pairs().size(), 36u);
  }
  else
  {
    const int local_search_count = pbc_local_search.get_pairs().size();
    int global_search_count=0;
    stk::all_reduce_sum(bulk_data.parallel(), &local_search_count, &global_search_count, 1);

    EXPECT_GE(global_search_count, 36);
  }
}


STKUNIT_UNIT_TEST(CoarseSearch, OffsetRotationalPeriodicBC)
{
  const double axisLocation[3] = {1.0, 0.0, 0.0};
  const stk::mesh::fixtures::GearMovement
      gear_movement_data(0.0, axisLocation[0], axisLocation[1], axisLocation[2]);

  stk::mesh::fixtures::GearsFixture fixture(MPI_COMM_SELF, 1);
  stk::mesh::BulkData & bulk_data = fixture.bulk_data;
  stk::mesh::MetaData & meta_data = fixture.meta_data;
  CoordFieldType & coords_field = fixture.cartesian_coord_field;
  CoordFieldType & disps_field = fixture.displacement_field;

  stk::mesh::Part & side_0 = meta_data.declare_part("side_0", stk::topology::NODE_RANK);
  stk::mesh::Part & side_1 = meta_data.declare_part("side_1", stk::topology::NODE_RANK);
  stk::mesh::PartVector side_0_vector(1, & side_0);
  stk::mesh::PartVector side_1_vector(1, & side_1);

  meta_data.commit();
  fixture.generate_mesh();

  //loop mesh entities and flag a wedge domain for putting in another part
  stk::mesh::fixtures::Gear & aGear = fixture.get_gear(0);
  aGear.move(gear_movement_data);

  typedef std::pair<stk::mesh::Entity,stk::mesh::Entity> NodePairType;
  std::vector<NodePairType> nodePairs;

  bulk_data.modification_begin();
  //side_0 at angle 0
  const size_t domainAngleIndex = 0;
  const size_t rangeAngleIndex = 3;

  for ( size_t ir = 0 ; ir < aGear.rad_num -2 ; ++ir ) {
    for ( size_t iz = 0 ; iz < aGear.height_num -1 ; ++iz ) {
      stk::mesh::Entity domainNode = aGear.get_node(iz, ir, domainAngleIndex);
      if ( bulk_data.is_valid(domainNode) && bulk_data.bucket(domainNode).owned() )
      {
        bulk_data.change_entity_parts(domainNode, side_0_vector);
      }
      stk::mesh::Entity rangeNode = aGear.get_node(iz, ir, rangeAngleIndex);
      if ( bulk_data.is_valid(rangeNode) && bulk_data.bucket(rangeNode).owned() )
      {
        bulk_data.change_entity_parts(rangeNode, side_1_vector);
      }
      nodePairs.push_back(std::make_pair(domainNode, rangeNode));
    }
  }
  bulk_data.modification_end();

  //do periodic search
  typedef stk::mesh::GetDisplacedCoordinates<CoordFieldType, CoordFieldType> DisplCoordsFunctor;
  typedef stk::mesh::PeriodicBoundarySearch<DisplCoordsFunctor> PeriodicSearch;
  PeriodicSearch pbc_search(bulk_data, DisplCoordsFunctor(bulk_data, coords_field, disps_field));

  //element_size is in radians
  const double rotationAngle = double(rangeAngleIndex - domainAngleIndex)*TWO_PI/aGear.angle_num;
  const double rotationAxis[3] = {0.0, 0.0, 1.0};

  pbc_search.add_rotational_periodic_pair(side_0 & meta_data.locally_owned_part(),
                                          side_1 & meta_data.locally_owned_part(),
                                          rotationAngle,
                                          rotationAxis,
                                          axisLocation);

  pbc_search.find_periodic_nodes(MPI_COMM_SELF);

  check_gold_rotational_multiperiodic(pbc_search.get_pairs());
  if (bulk_data.parallel_size() == 1)
  {
    EXPECT_EQ(pbc_search.get_pairs().size(), 36u);
  }
  else
  {
    const int local_search_count = pbc_search.get_pairs().size();
    int global_search_count=0;
    stk::all_reduce_sum(bulk_data.parallel(), &local_search_count, &global_search_count, 1);

    EXPECT_GE(global_search_count, 36);
  }

  //now we ghost everything to do a local search
  bulk_data.modification_begin();
  pbc_search.create_ghosting("periodic_ghosts");
  bulk_data.modification_end();

  const stk::mesh::Ghosting & periodic_bc_ghosting = pbc_search.get_ghosting();

  std::vector< stk::mesh::FieldBase const * > ghosted_fields;
  ghosted_fields.push_back(&coords_field);
  stk::mesh::communicate_field_data( periodic_bc_ghosting, ghosted_fields);

  //lets do a local search to make sure the coords field and entities were properly ghosted
  PeriodicSearch pbc_local_search(bulk_data, DisplCoordsFunctor(bulk_data, coords_field, disps_field));
  pbc_local_search.add_rotational_periodic_pair(side_0 & meta_data.locally_owned_part(),
                                          side_1 & meta_data.locally_owned_part(),
                                          rotationAngle,
                                          rotationAxis,
                                          axisLocation);

  pbc_local_search.find_periodic_nodes(MPI_COMM_SELF, true);

  //test to make sure it is re-entrant
  pbc_local_search.find_periodic_nodes(MPI_COMM_SELF, true);

  check_gold_rotational_multiperiodic(pbc_local_search.get_pairs());
  if (bulk_data.parallel_size() == 1)
  {
    EXPECT_EQ(pbc_local_search.get_pairs().size(), 36u);
  }
  else
  {
    const int local_search_count = pbc_local_search.get_pairs().size();
    int global_search_count=0;
    stk::all_reduce_sum(bulk_data.parallel(), &local_search_count, &global_search_count, 1);

    EXPECT_GE(global_search_count, 36);
  }
}
