/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/
#include <use_cases/UseCase_Skinning.hpp>
#include <stk_mesh/fixtures/HexFixture.hpp>

#include <stk_mesh/base/BulkModification.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Selector.hpp>
#include <stk_mesh/base/GetBuckets.hpp>
#include <stk_mesh/base/EntityComm.hpp>

#include <stk_mesh/fem/EntityRanks.hpp>
#include <stk_mesh/fem/TopologyHelpers.hpp>
#include <stk_mesh/fem/BoundaryAnalysis.hpp>



#include <stk_util/parallel/ParallelReduce.hpp>

#include <stk_mesh/fem/SkinMesh.hpp>

namespace {

  unsigned count_skin_entities( stk::mesh::fixtures::HexFixture3x3x3 & fixture) {

    const unsigned mesh_rank = stk::mesh::Element;

    stk::mesh::BulkData & mesh = fixture.m_bulk_data;
    stk::mesh::MetaData & meta = fixture.m_meta_data;
    stk::mesh::Part & skin_part = fixture.m_skin_part;

    stk::mesh::Selector select_skin = skin_part & meta.locally_owned_part()  ;

    const std::vector<stk::mesh::Bucket*>& buckets = mesh.buckets( mesh_rank -1);

    return count_selected_entities( select_skin, buckets);
  }

}

bool skinning_use_case_1(stk::ParallelMachine pm)
{
  // TODO This test should be able to run with up to 27
  // processes, but it fails with corrupt communication
  // with 10 or more processes
  const unsigned p_size = stk::parallel_machine_size( pm );
  if ( 4 < p_size ) { return true; }

  //setup the mesh
  stk::mesh::fixtures::HexFixture3x3x3 fixture(pm);

  stk::mesh::MetaData& meta = fixture.m_meta_data;
  stk::mesh::BulkData& mesh = fixture.m_bulk_data;

  stk::mesh::Part & skin_part = fixture.m_skin_part;

  std::vector< stk::mesh::EntityId > elements_to_separate;

  //separate out the middle element
  elements_to_separate.push_back(fixture.m_elems_id[1][1][1]);

  separate_and_skin_mesh(
      meta,
      mesh,
      skin_part,
      elements_to_separate
      );


  // pointer to middle_element after mesh modification.
  stk::mesh::Entity * middle_element = mesh.get_entity(stk::mesh::Element,fixture.m_elems_id[1][1][1]);

  unsigned num_skin_entities = count_skin_entities(fixture);

  stk::all_reduce(pm, stk::ReduceSum<1>(&num_skin_entities));

  //there should be 66 faces in the skin part
  //54 on the outside
  //6 on the inside attached to the entire mesh
  //6 on the inside attected to the element that was detached
  bool correct_skin = ( num_skin_entities == 66 );
  bool correct_relations = true;
  bool correct_comm = true;

  //all nodes connected to the single element that has been broken off
  //should have relations.size() == 4 and comm.size() == 0
  if (middle_element != NULL && middle_element->owner_rank() == mesh.parallel_rank()) {

      stk::mesh::PairIterRelation relations = middle_element->relations(stk::mesh::Node);

      for (; relations.first != relations.second; ++relations.first) {
        stk::mesh::Entity * current_node = (relations.first->entity());
        //each node should be attached to only 1 element and 3 faces
        correct_relations &= ( current_node->relations().size() == 4 );
        //the entire closure of the element should exist on a single process
        correct_comm      &= ( current_node->comm().size() == 0 );
      }
  }

  return (correct_skin && correct_relations && correct_comm);
}
