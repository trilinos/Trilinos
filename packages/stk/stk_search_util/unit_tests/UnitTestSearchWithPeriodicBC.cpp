#include <stk_util/unit_test_support/stk_utest_macros.hpp>
#include <stk_search/CoarseSearch.hpp>
#include <stk_search/BoundingBox.hpp>
#include <stk_search/IdentProc.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/FieldParallel.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_mesh/fixtures/HexFixture.hpp>
#include <iostream>
#include <iomanip>
#include <sstream>

typedef double Scalar;
typedef stk::search::ident::IdentProc<stk::mesh::EntityId,unsigned> SearchId;
typedef stk::search::box::SphereBoundingBox<SearchId,Scalar,3> AABB;
typedef std::vector<AABB> AABBVector;
typedef std::vector<std::pair<SearchId,SearchId> > SearchPairVector;
typedef stk::mesh::fixtures::HexFixture::CoordFieldType CoordFieldType;

namespace {

size_t populate_search_vector( stk::mesh::BulkData & bulk_data
                              ,CoordFieldType & coords_field
                              ,stk::mesh::Selector side_selector
                              ,AABBVector & aabb_vector
                              ,double * centroid
                             )
{
  const int spatial_dimension = bulk_data.mesh_meta_data().spatial_dimension();
  const double radius = 1e-10;
  const unsigned parallel_rank = bulk_data.parallel_rank();

  stk::mesh::BucketVector buckets;
  stk::mesh::get_buckets( side_selector
                         ,bulk_data.buckets(stk::topology::NODE_RANK)
                         ,buckets
                        );

  size_t num_nodes = 0;

  for (size_t bindex = 0, num_buckets = buckets.size(); bindex < num_buckets; ++bindex) {
    stk::mesh::Bucket & b = *buckets[bindex];
    for (size_t ord =0, num_entities = b.size(); ord < num_entities; ++ord) {
      ++num_nodes;
      const double * const coords = bulk_data.field_data(coords_field, b, ord);
      SearchId search_id( bulk_data.identifier(b[ord]), parallel_rank);
      aabb_vector.push_back(AABB( coords, radius, search_id));
      for (int i=0; i<spatial_dimension; ++i) {
        centroid[i] += coords[i];
      }
    }
  }
  return num_nodes;
}

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

  for (size_t i=0, size=search_results.size(); i<size; ++i) {
    stk::mesh::EntityId domain_node = search_results[i].first.ident;
    stk::mesh::EntityId range_node = search_results[i].second.ident;

    EXPECT_TRUE((std::lower_bound(gold.begin(),gold.end(),std::make_pair(domain_node,range_node))) != gold.end());
  }
}

void find_periodic_nodes( stk::mesh::BulkData & bulk_data
                         ,stk::ParallelMachine parallel
                         ,CoordFieldType & coords_field
                         ,stk::mesh::Selector side1
                         ,stk::mesh::Selector side2
                         ,SearchPairVector & search_results
                        )
{
  AABBVector side_1_vector, side_2_vector;

  double local_centroid[6] = {0};
  size_t local_node_count[2];

  local_node_count[0] =
    populate_search_vector( bulk_data
                           ,coords_field
                           ,side1
                           ,side_1_vector
                           ,local_centroid
                          );

  local_node_count[1] =
    populate_search_vector( bulk_data
                           ,coords_field
                           ,side2
                           ,side_2_vector
                           ,local_centroid + 3
                           );

  size_t global_node_count[2] = {0};
  double global_centroid[6] = {0};

  stk::all_reduce_sum( parallel, local_node_count, global_node_count, 2);
  stk::all_reduce_sum( parallel, local_centroid, global_centroid, 6);

  // TODO pass in translation method...
  double translate[3];
  for (int i=0; i<3; ++i) {
    translate[i] =   (global_centroid[i+3]/global_node_count[1])
                   - (global_centroid[i]/global_node_count[0]);
  }

  // translate domain to range, i.e. master to slave
  for(size_t i=0, size=side_1_vector.size(); i<size; ++i) {
    for (int j=0; j<3; ++j) {
      side_1_vector[i].center[j] += translate[j];
    }
  }

  stk::search::FactoryOrder order;
  order.m_communicator = parallel;
  stk::search::coarse_search( search_results
                             ,side_2_vector
                             ,side_1_vector
                             ,order
                            );

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

  SearchPairVector search_results;

  find_periodic_nodes( bulk_data
                      ,bulk_data.parallel()
                      ,coords_field
                      ,side_0 & meta_data.locally_owned_part()
                      ,side_3 & meta_data.locally_owned_part()
                      ,search_results
                     );


  check_gold(search_results);


  const int parallel_rank = bulk_data.parallel_rank();
  size_t num_constraints = 0;
  std::vector<stk::mesh::EntityProc> send_nodes;
  for (size_t i=0, size=search_results.size(); i<size; ++i) {
      stk::mesh::Entity domain_node = bulk_data.get_entity(stk::topology::NODE_RANK, search_results[i].first.ident);
      int domain_proc = search_results[i].first.proc;
      stk::mesh::Entity range_node = bulk_data.get_entity(stk::topology::NODE_RANK, search_results[i].second.ident);
      int range_proc = search_results[i].second.proc;

      if (parallel_rank == domain_proc) ++num_constraints;

      if ((parallel_rank != domain_proc) && (parallel_rank == range_proc)) {
        send_nodes.push_back(stk::mesh::EntityProc(range_node, domain_proc));
      }
      else if ((parallel_rank == domain_proc) && (parallel_rank != range_proc)) {
        send_nodes.push_back(stk::mesh::EntityProc(domain_node, range_proc));
      }
      ThrowAssert((parallel_rank == domain_proc) || (parallel_rank == range_proc));
  }

  bulk_data.modification_begin();
  stk::mesh::Ghosting & periodic_bc_ghosting = bulk_data.create_ghosting( "periodic_bc");
  bulk_data.change_ghosting(periodic_bc_ghosting, send_nodes);

  bulk_data.modification_end();

  std::vector< stk::mesh::FieldBase const * > ghosted_fields;
  ghosted_fields.push_back(&coords_field);
  stk::mesh::communicate_field_data( periodic_bc_ghosting, ghosted_fields);

  SearchPairVector local_search_results;

  find_periodic_nodes( bulk_data
                      ,MPI_COMM_SELF
                      ,coords_field
                      ,side_0
                      ,side_3
                      ,local_search_results
                     );

  check_gold(local_search_results);

}

