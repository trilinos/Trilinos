/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <iostream>
#include <sstream>
#include <stdexcept>

#include <stk_util/unit_test_support/stk_utest_macros.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/EntityComm.hpp>
#include <stk_mesh/base/Comm.hpp>

#include <stk_mesh/fixtures/BoxFixture.hpp>
#include <stk_mesh/fixtures/RingFixture.hpp>

#include <unit_tests/UnitTestModificationEndWrapper.hpp>
#include <unit_tests/UnitTestRingFixture.hpp>

using stk_classic::mesh::Part;
using stk_classic::mesh::Bucket;
using stk_classic::mesh::PairIterRelation;
using stk_classic::mesh::PairIterEntityComm;
using stk_classic::mesh::MetaData;
using stk_classic::mesh::fem::FEMMetaData;
using stk_classic::mesh::BulkData;
using stk_classic::mesh::Selector;
using stk_classic::mesh::PartVector;
using stk_classic::mesh::BaseEntityRank;
using stk_classic::mesh::PairIterRelation;
using stk_classic::mesh::EntityProc;
using stk_classic::mesh::Entity;
using stk_classic::mesh::EntityId;
using stk_classic::mesh::EntityKey;
using stk_classic::mesh::EntityVector;
using stk_classic::mesh::EntityRank;
using stk_classic::mesh::fixtures::RingFixture;
using stk_classic::mesh::fixtures::BoxFixture;

namespace {
const EntityRank NODE_RANK = FEMMetaData::NODE_RANK;
} // empty namespace

void printEntity(std::ostringstream& msg, Entity *entity)
{
  msg << " :: " << print_entity_key(entity) << ":o[" << entity->owner_rank() << "]:l[" << entity->log_query()
      << "]:ec[";
  for ( PairIterEntityComm ec = entity->comm() ; ! ec.empty() ; ++ec ) {
    msg << "(" << ec->ghost_id << "," << ec->proc << ")";
  }
  msg << "]";
}

void printNode(std::ostringstream& msg, Entity *node)
{
  printEntity(msg, node);
  PairIterRelation rels = node->relations();
  for (unsigned i = 0; i < rels.size(); i++)
    {
      Entity *entity = rels[i].entity();
      if (entity->entity_rank() > node->entity_rank())
        printEntity(msg, entity);
    }
}

void printBuckets(std::ostringstream& msg, BulkData& mesh)
{
  const std::vector<Bucket*> & buckets = mesh.buckets(0);
  for (unsigned i=0; i < buckets.size(); i++)
    {
      const Bucket& bucket = *buckets[i];
      msg << " bucket[" << i << "] = ";
      size_t bucket_size = bucket.size();
      for (unsigned ie=0; ie < bucket_size; ie++)
        {
          msg << bucket[ie].identifier() << ", ";
        }
    }
}

static void checkBuckets( BulkData& mesh)
{
  const std::vector<Bucket*> & buckets = mesh.buckets(0);
  for (unsigned i=0; i < buckets.size(); i++)
    {
      Bucket* bucket = buckets[i];
      STKUNIT_ASSERT(bucket->assert_correct());
    }
}

STKUNIT_UNIT_TEST(UnitTestingOfBulkData, test_other_ghosting_2)
{
  //
  // testing if modification flags propagate properly for ghosted entities
  //
  // To test this, we focus on a single node shared on 2 procs, ghosted on others
  //

  /** 
   * 1D Mesh (node,owner)--[elem,owner]---(...)
   *
   * <---(70,0)--[500,0]--(41,1)--[301,1]---(42,2)---[402,2]---(70,0)--->
   *                              
   * <---(50,0)--[100,0]--(21,1)--[201,1]---(32,2)---[302,2]---(50,0)--->
   * 
   */

  // elem, node0, node1, owner
  EntityId elems_0[][4] = { {100, 21, 50, 0}, {201, 21, 32, 1}, {302, 32, 50, 2}, 
                            {500, 41, 70, 0}, {301, 41, 42, 1}, {402, 42, 70, 2}  };
  // node, owner
  EntityId nodes_0[][2] = { {21,1}, {50,0}, {32, 2}, {41, 1}, {42, 1}, {70, 0} };

  unsigned nelems = sizeof(elems_0)/4/sizeof(EntityId);
  unsigned nnodes = sizeof(nodes_0)/2/sizeof(EntityId);

  stk_classic::ParallelMachine pm = MPI_COMM_WORLD;

  // Set up meta and bulk data
  const unsigned spatial_dim = 2;

  std::vector<std::string> entity_rank_names = stk_classic::mesh::fem::entity_rank_names(spatial_dim);
  entity_rank_names.push_back("FAMILY_TREE");

  FEMMetaData meta_data(spatial_dim, entity_rank_names);
  //Part & part_tmp = meta_data.declare_part( "temp");
  
  meta_data.commit();
  unsigned max_bucket_size = 1;
  BulkData mesh(FEMMetaData::get_meta_data(meta_data), pm, max_bucket_size);
  //BulkData mesh(FEMMetaData::get_meta_data(meta_data), pm);
  unsigned p_rank = mesh.parallel_rank();
  unsigned p_size = mesh.parallel_size();

  if (p_size != 3) return;

  //
  // Begin modification cycle so we can create the entities and relations
  //

  // We're just going to add everything to the universal part
  stk_classic::mesh::PartVector empty_parts;

  // Create elements
  const EntityRank elem_rank = meta_data.element_rank();
  Entity * elem = 0;

  mesh.modification_begin();

  for (unsigned ielem=0; ielem < nelems; ielem++)
    {
      if (elems_0[ielem][3] == p_rank)
        {
          elem = &mesh.declare_entity(elem_rank, elems_0[ielem][0], empty_parts);

          EntityVector nodes;
          // Create node on all procs
          nodes.push_back( &mesh.declare_entity(NODE_RANK, elems_0[ielem][2], empty_parts) );
          nodes.push_back( &mesh.declare_entity(NODE_RANK, elems_0[ielem][1], empty_parts) );

          // Add relations to nodes
          mesh.declare_relation( *elem, *nodes[0], 0 );
          mesh.declare_relation( *elem, *nodes[1], 1 );

        }
    }

  mesh.modification_end();

  Entity* node1 = 0;

  // change node owners
  mesh.modification_begin();

  std::vector<EntityProc> change;

  for (unsigned inode=0; inode < nnodes; inode++)
    {
      node1 = mesh.get_entity(0, nodes_0[inode][0]);
      if (node1 && node1->owner_rank() == p_rank)
        {
          unsigned dest = nodes_0[inode][1];
          EntityProc eproc(node1, dest);
          change.push_back(eproc);
        }
    }

  mesh.change_entity_owner( change );

  mesh.modification_end();

  checkBuckets(mesh);

  MPI_Barrier(MPI_COMM_WORLD);


  // attempt to delete a node and its elems but on a ghosted proc
  mesh.modification_begin();

  if (p_rank == 2)
    {
      node1 = mesh.get_entity(0, 21);
      Entity *elem1 = mesh.get_entity(2, 201);
      Entity *elem2 = mesh.get_entity(2, 100);

      bool did_it_elem = mesh.destroy_entity(elem1);
      did_it_elem = did_it_elem & mesh.destroy_entity(elem2);
      STKUNIT_ASSERT(did_it_elem);
      bool did_it = mesh.destroy_entity(node1);
      STKUNIT_ASSERT(did_it);
    }

  mesh.modification_end();

  checkBuckets(mesh);

  // this node should no longer exist anywhere
  node1 = mesh.get_entity(0, 21);

  // uncomment to force failure of test
  // STKUNIT_ASSERT(node1 == 0);

}

