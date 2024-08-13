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
#ifndef _UnitTestModificationEnd_hpp_
#define _UnitTestModificationEnd_hpp_

#include <gtest/gtest.h>
#include <stddef.h>  // for size_t
#include <stdlib.h>  // for exit

#include <exception>                   // for exception
#include <iostream>                    // for ostringstream, etc
#include <iterator>                    // for distance
#include <map>                         // for _Rb_tree_const_iterator, etc
#include <stdexcept>                   // for logic_error, runtime_error
#include <stk_mesh/base/BulkData.hpp>  // for BulkData, etc
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/EntityLess.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>   // for BulkData, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/Types.hpp>                 // for MeshIndex, EntityRank, etc
#include <stk_mesh/base/DumpMeshInfo.hpp>
#include <stk_mesh/baseImpl/BucketRepository.hpp>  // for BucketRepository
#include <stk_mesh/baseImpl/MeshImplUtils.hpp>
#include <stk_mesh/baseImpl/elementGraph/ElemElemGraph.hpp>
#include <stk_io/StkMeshIoBroker.hpp>
#include <stk_unit_test_utils/BulkDataTester.hpp>
#include <stk_unit_test_utils/CommandLineArgs.hpp>

namespace stk { namespace mesh { namespace unit_test {

void populateBulkDataWithFile(const std::string& exodusFileName, MPI_Comm communicator, stk::unit_test_util::BulkDataTester& bulkData);
void checkCommListAndMap(const stk::unit_test_util::BulkDataTester& stkMeshBulkData, bool isAfterIGMD);

void checkStatesOfEntities(std::vector<std::vector<stk::mesh::EntityState> > &nodeStates,
                           std::vector<std::vector<stk::mesh::EntityState> > &elementStates,
                           bool (&areNodesValid)[2][20], bool (&areElementsValid)[2][4],
stk::unit_test_util::BulkDataTester &stkMeshBulkData);
void checkThatMeshIsParallelConsistent(stk::unit_test_util::BulkDataTester& stkMeshBulkData);
void mark_element3_as_modified(stk::unit_test_util::BulkDataTester& stkMeshBulkData);
void makeSureEntityIsValidOnCommListAndBulkData(stk::unit_test_util::BulkDataTester& stkMeshBulkData, stk::mesh::EntityKey &entityKey);
void makeSureEntityIsValidOnCommListAndBut_NOT_BulkData(stk::unit_test_util::BulkDataTester& stkMeshBulkData, stk::mesh::EntityKey &entityKey);
void getMeshLineByLine(const stk::unit_test_util::BulkDataTester &stkMeshBulkData, std::vector<std::string> &output);
void checkCommMapsAndLists(stk::unit_test_util::BulkDataTester& stkMeshBulkData);
void destroy_element3_on_proc_1(stk::unit_test_util::BulkDataTester& stkMeshBulkData, stk::mesh::EntityKey &elementToDestroyKey);
void checkCommMapsAndListsAfterIRSMD(stk::unit_test_util::BulkDataTester& stkMeshBulkData);
void checkCommMapsAndListsAfterIRGMD(stk::unit_test_util::BulkDataTester& stkMeshBulkData);
void create_edges(stk::unit_test_util::BulkDataTester& stkMeshBulkData,
                  std::vector<stk::mesh::EntityId>& edgeIds,
                  std::vector<std::vector<stk::mesh::EntityId> > &nodeIdsForEdge,
                  std::vector<std::vector<stk::mesh::EntityId> > &elementRelations,
                  std::vector<stk::mesh::Entity> &edgeEntities, stk::mesh::Part& edge_part);

void checkResultsOfIRSMD_for_edges(stk::unit_test_util::BulkDataTester &stkMeshBulkData);
void checkResultsOfIRGMD_for_edges(stk::unit_test_util::BulkDataTester &stkMeshBulkData, std::vector<stk::mesh::Entity> &edgeEntities);
void checkItAllForThisCase(stk::unit_test_util::BulkDataTester &stkMeshBulkData);
void checkItAllForThisGhostedCase(stk::unit_test_util::BulkDataTester &stkMeshBulkData);


inline std::string getOption(const std::string& option, const std::string defaultString="no")
{
  std::string returnValue = defaultString;
  stk::unit_test_util::GlobalCommandLineArguments &args = stk::unit_test_util::GlobalCommandLineArguments::self();
   if (args.get_argv() != nullptr) {
    for (int i = 0; i < args.get_argc(); i++) {
      std::string input_argv(args.get_argv()[i]);
      if (option == input_argv) {
        if ((i + 1) < args.get_argc()) {
          returnValue = std::string(args.get_argv()[i + 1]);
        }
        break;
      }
    }
  }
  return returnValue;
}

// Write out vector of strings using proc id and label via an ostringstream
//void writeMesh(int myProcId, std::string label, const std::vector<std::string> &meshStart)
//{
//    std::ostringstream msg;
//    for (size_t i=0;i<meshStart.size();i++)
//    {
//        msg.str(std::string());
//        msg << "P[" << myProcId << "] " << label << "\t" << meshStart[i] << std::endl;
//        std::cerr << msg.str();
//    }
//}

void getMeshLineByLine(const stk::unit_test_util::BulkDataTester &stkMeshBulkData, std::vector<std::string> &output)
{
  std::ostringstream msg;
  impl::dump_all_mesh_info(stkMeshBulkData, msg);
  std::istringstream iss(msg.str());

  std::string s;
  while ( std::getline(iss, s) )
  {
    output.push_back(s);
  }
}

void populateBulkDataWithFile(const std::string& exodusFileName, MPI_Comm communicator, stk::unit_test_util::BulkDataTester& bulkData)
// STK IO module will be described in separate chapter.
// It is used here to read the mesh data from the Exodus file and populate an STK Mesh.
// The order of the following lines in {} are important
{
  stk::io::StkMeshIoBroker exodusFileReader(communicator);

  // Inform STK IO which STK Mesh objects to populate later
  exodusFileReader.set_bulk_data(bulkData);

  exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);

  // Populate the MetaData which has the descriptions of the Parts and Fields.
  exodusFileReader.create_input_mesh();

  // Populate entities in STK Mesh from Exodus file
  exodusFileReader.populate_bulk_data();
}

void checkCommListAndMap(const stk::unit_test_util::BulkDataTester& stkMeshBulkData, bool isAfterIRGMD)
{
  bool doesElementAppearInAuraCommMap[2][4] = { { false, true, true, false, },
                                                { false, true, true, false } };

  bool isElementValidInCommList[2][4] = { { false, true, true, false, },
                                          { false, true, true, false } };

  if ( isAfterIRGMD )
  {
    doesElementAppearInAuraCommMap[0][1] = false;
    doesElementAppearInAuraCommMap[0][2] = false;
    doesElementAppearInAuraCommMap[1][2] = false;
    isElementValidInCommList[0][2] = false;
    isElementValidInCommList[0][2] = false;
  }

  stk::mesh::EntityCommListInfoVector::const_iterator iter;
  int myProcId = stkMeshBulkData.parallel_rank();
  for (unsigned int element=1;element<=4;element++)
  {
    int elementoffset = element-1;
    stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
    iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), elementKey);
    if ( iter == stkMeshBulkData.my_internal_comm_list().end() || iter->key != elementKey )
    {
      EXPECT_EQ(isElementValidInCommList[myProcId][elementoffset], false) <<
                                                                             "Proc " << myProcId << " for element " << element;
    }
    else
    {
      EXPECT_EQ(isElementValidInCommList[myProcId][elementoffset], stkMeshBulkData.is_valid(iter->entity)) <<
                                                                                                              "proc " << myProcId << " using offset " << elementoffset;
    }

    const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
    EXPECT_EQ(doesElementAppearInAuraCommMap[myProcId][elementoffset], is_entity_ghosted ) <<
                                                                                              "proc " << myProcId << " using offset " << elementoffset;
  }
}


void checkStatesOfEntities(std::vector<std::vector<stk::mesh::EntityState> > &nodeStates,
                           std::vector<std::vector<stk::mesh::EntityState> > &elementStates,
                           bool (&areNodesValid)[2][20], bool (&areElementsValid)[2][4],
stk::unit_test_util::BulkDataTester &stkMeshBulkData)
{
  int myProcId = stkMeshBulkData.parallel_rank();
  for (unsigned nodeId=1;nodeId<=nodeStates[myProcId].size();nodeId++)
  {
    stk::mesh::EntityKey entity_key(stk::topology::NODE_RANK, nodeId);
    stk::mesh::Entity entity = stkMeshBulkData.get_entity(entity_key);
    ASSERT_EQ(areNodesValid[myProcId][nodeId-1], stkMeshBulkData.is_valid(entity)) <<
                                                                                      "Proc " << myProcId << " using node " << nodeId;
    if ( stkMeshBulkData.is_valid(entity) )
    {
      EXPECT_EQ(nodeStates[myProcId][nodeId-1], stkMeshBulkData.state(entity));
    }
  }

  for (unsigned elementId=1;elementId<=elementStates[myProcId].size();elementId++)
  {
    stk::mesh::EntityKey entity_key(stk::topology::ELEMENT_RANK, elementId);
    stk::mesh::Entity entity = stkMeshBulkData.get_entity(entity_key);
    ASSERT_EQ(areElementsValid[myProcId][elementId-1], stkMeshBulkData.is_valid(entity));
    if ( stkMeshBulkData.is_valid(entity) )
    {
      EXPECT_EQ(elementStates[myProcId][elementId-1], stkMeshBulkData.state(entity)) ;
    }
  }
}

void mark_element3_as_modified(stk::unit_test_util::BulkDataTester& stkMeshBulkData)
{
  int elementToModify = 3;
  stk::mesh::EntityKey elementToModifyKey(stk::topology::ELEMENT_RANK, elementToModify);
  stk::mesh::Entity entity = stkMeshBulkData.get_entity(elementToModifyKey);

  if ( stkMeshBulkData.parallel_rank() == 1 )
  {
    stkMeshBulkData.my_set_state(entity, stk::mesh::Modified);
  }
}

stk::mesh::EntityCommListInfoVector::const_iterator makeSureEntityIsValidOnCommList(stk::unit_test_util::BulkDataTester& stkMeshBulkData, stk::mesh::EntityKey &entityKey)
{
  stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), entityKey);
  EXPECT_TRUE(iter != stkMeshBulkData.my_internal_comm_list().end());
  EXPECT_EQ(entityKey,iter->key) << " looking for entityKey in comm_list. Did not find.";
  return iter;
}

stk::mesh::EntityCommListInfoVector::const_iterator makeSureEntityIs_NOT_ValidOnCommList(stk::unit_test_util::BulkDataTester& stkMeshBulkData, stk::mesh::EntityKey &entityKey)
{
  stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), entityKey);
  EXPECT_TRUE(iter == stkMeshBulkData.my_internal_comm_list().end() || entityKey != iter->key);
  return iter;
}

void makeSureEntityIsValidOnCommListAndBut_NOT_BulkData(stk::unit_test_util::BulkDataTester& stkMeshBulkData, stk::mesh::EntityKey &entityKey)
{
  stk::mesh::EntityCommListInfoVector::const_iterator iter = makeSureEntityIsValidOnCommList(stkMeshBulkData, entityKey);
  EXPECT_FALSE(stkMeshBulkData.is_valid(iter->entity));
}

void makeSureEntityIsValidOnCommListAndBulkData(stk::unit_test_util::BulkDataTester& stkMeshBulkData, stk::mesh::EntityKey &entityKey)
{
  stk::mesh::EntityCommListInfoVector::const_iterator iter = makeSureEntityIsValidOnCommList(stkMeshBulkData, entityKey);
  EXPECT_TRUE(stkMeshBulkData.is_valid(iter->entity));
}

void checkThatMeshIsParallelConsistent(stk::unit_test_util::BulkDataTester& stkMeshBulkData)
{
  std::ostringstream msg ;
  bool is_consistent = true;
  is_consistent = stkMeshBulkData.my_comm_mesh_verify_parallel_consistency( msg );
  EXPECT_TRUE(is_consistent) << msg.str();
}

void checkCommMapsAndLists(stk::unit_test_util::BulkDataTester& stkMeshBulkData)
{
  bool isNodeInCommList[2][20] = {
    {false, false, false, false,
     true,  true,  true,  true,
     true,  true,  true,  true,
     true,  true,  true,  true,
     false, false, false, false,},

    {false, false, false, false,
     true,  true,  true,  true,
     true,  true,  true,  true,
     true,  true,  true,  true,
     false, false, false, false}
  };

  bool isNodeInGhostedCommMap[2][20] = {
    {false, false, false, false,
     true,  true,  true,  true,
     false, false, false, false,
     true,  true,  true,  true,
     false, false, false, false,
    },
    {false, false, false, false,
     true,  true,  true,  true,
     false, false, false, false,
     true,  true,  true,  true,
     false, false, false, false,
    }
  };

  bool IsNodeInSharedCommMap[2][20] = {
    {false, false, false, false,
     false, false, false, false,
     true,  true,  true,  true,
     false, false, false, false,
     false, false, false, false,
    },
    {false, false, false, false,
     false, false, false, false,
     true,  true,  true,  true,
     false, false, false, false,
     false, false, false, false,
    }
  };

  for (unsigned int nodeId=1;nodeId<=20;nodeId++)
  {
    stk::mesh::EntityKey nodeKey(stk::topology::NODE_RANK,nodeId);
    stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), nodeKey);

    if ( iter != stkMeshBulkData.my_internal_comm_list().end() && iter->key == nodeKey )
    {
      EXPECT_TRUE(isNodeInCommList[stkMeshBulkData.parallel_rank()][nodeId-1]) <<
                                                                                  "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
    }
    else
    {
      EXPECT_FALSE(isNodeInCommList[stkMeshBulkData.parallel_rank()][nodeId-1]) <<
                                                                                   "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
    }

    const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(nodeKey, stkMeshBulkData.aura_ghosting()).empty();
    EXPECT_EQ( isNodeInGhostedCommMap[stkMeshBulkData.parallel_rank()][nodeId-1], is_entity_ghosted ) <<
                                                                                                         "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;

    const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(nodeKey, stkMeshBulkData.shared_ghosting()).empty();
    EXPECT_EQ( IsNodeInSharedCommMap[stkMeshBulkData.parallel_rank()][nodeId-1], is_entity_shared ) <<
                                                                                                       "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
  }
}

void destroy_element3_on_proc_1(stk::unit_test_util::BulkDataTester& stkMeshBulkData, stk::mesh::EntityKey &elementToDestroyKey)
{
  stkMeshBulkData.modification_begin();

  ASSERT_TRUE ( stkMeshBulkData.is_valid(stkMeshBulkData.get_entity(elementToDestroyKey)) );

  if ( stkMeshBulkData.parallel_rank() == 1 )
  {
    stkMeshBulkData.destroy_entity( stkMeshBulkData.get_entity(elementToDestroyKey) );
  }
}

void checkCommMapsAndListsAfterIRSMD(stk::unit_test_util::BulkDataTester& stkMeshBulkData)
{
  //============== Checking result of IRSMD


  /*
     * We deleted element 3 from mesh that looks like E1 : E2 : E3 : E4
     * Then we called internal_resolve_shared_modify_delete. Note,
     * nodes 9 thru 12 are on the boundary between E2 and E3. Also, proc 0
     * has E1 : E2 and proc 1 has E3 : E4. So the proc decomp looks like
     *
     * proc 0                               proc 1
     *    E1 : E2 : GE3                             GE2 : E3 : E4
     *
     * where GE2 and GE3 are ghosted elements.
     *
     * Result:
     *      For proc 0: nodes 9 thru 12 are still valid entities
     *      For proc 1: nodes 9 thru 12 are not valid anymore
     *
     *      Proc 0 still *thinks* it needes to communicate about elements 1, 2, and 3. Why?
     *      Proc 1 doesn't have any elements it needs to communicate with.
     *
     *      Proc 0 believes elements 2 and 3 are still needed for communication for aura.
     *      Proc 1 believes elements 2 and 3 are still needed for communication for aura.
     *
     * Question: Ghosting is NOT right, but is sharing right?
     *      According to test below, nothing is shared. So yes!?
     */

  bool isNodeValidInCommList[2][4] = { { true, true, true, true, },
                                       { false, false, false, false } };

  bool isElementValidInCommList[2][4] = { { false, true, true, false, },
                                          { false, false, false, false } };

  bool isElementInAuraCommMap[2][4] = { { false, true, true, false, },
                                        { false, true, true, false } };

  //============== Testing against the boolean vectors

  for (unsigned int node=9;node<=12;node++)
  {
    int nodeoffset = node-9;
    stk::mesh::EntityKey nodeEntityKey(stk::topology::NODE_RANK,node);
    stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), nodeEntityKey);
    if (iter == stkMeshBulkData.my_internal_comm_list().end() || iter->key != nodeEntityKey)
    {
      EXPECT_EQ(isNodeValidInCommList[stkMeshBulkData.parallel_rank()][nodeoffset], false);
    }
    else
    {
      EXPECT_EQ(isNodeValidInCommList[stkMeshBulkData.parallel_rank()][nodeoffset], stkMeshBulkData.is_valid(iter->entity));
    }

    const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(nodeEntityKey, stkMeshBulkData.aura_ghosting()).empty();
    const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(nodeEntityKey, stkMeshBulkData.shared_ghosting()).empty();
    EXPECT_FALSE( is_entity_shared );
    EXPECT_FALSE( is_entity_ghosted );
  }

  for (unsigned int element=1;element<=4;element++)
  {
    int elementoffset = element-1;
    stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
    stk::mesh::EntityCommListInfoVector::const_iterator  iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), elementKey);
    if (iter == stkMeshBulkData.my_internal_comm_list().end() || iter->key != elementKey)
    {
      EXPECT_EQ(isElementValidInCommList[stkMeshBulkData.parallel_rank()][elementoffset], false);
    }
    else
    {
      EXPECT_EQ(isElementValidInCommList[stkMeshBulkData.parallel_rank()][elementoffset], stkMeshBulkData.is_valid(iter->entity));
    }

    const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
    EXPECT_EQ(isElementInAuraCommMap[stkMeshBulkData.parallel_rank()][elementoffset], is_entity_ghosted);

    const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
    EXPECT_FALSE( is_entity_shared );
  }
}

void checkCommMapsAndListsAfterIRGMD(stk::unit_test_util::BulkDataTester& stkMeshBulkData)
{
  bool isElementValidInCommListAfterIRGMD[2][4] = { { false, true, false, false, },
                                                    { false, false, false, false } };

  //============== Check results against boolean above

  for (unsigned int element=1;element<=4;element++)
  {
    int elementoffset = element-1;
    stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
    stk::mesh::EntityCommListInfoVector::const_iterator  iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), elementKey);
    if (iter == stkMeshBulkData.my_internal_comm_list().end() || iter->key != elementKey)
    {
      EXPECT_EQ(isElementValidInCommListAfterIRGMD[stkMeshBulkData.parallel_rank()][elementoffset], false);
    }
    else
    {
      EXPECT_EQ(isElementValidInCommListAfterIRGMD[stkMeshBulkData.parallel_rank()][elementoffset], stkMeshBulkData.is_valid(iter->entity));
    }

    const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
    EXPECT_FALSE( is_entity_ghosted );

    const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
    EXPECT_FALSE( is_entity_shared );
  }

  //============== Check results of updating the comm list based on changes in the comm map

  stkMeshBulkData.my_update_comm_list_based_on_changes_in_comm_map();

  //============== No element is ghosted, No node is shared, and comm list is empty

  for (unsigned int element=1;element<=4;element++)
  {
    stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
    stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), elementKey);
    EXPECT_TRUE(iter == stkMeshBulkData.my_internal_comm_list().end());

    const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
    EXPECT_FALSE( is_entity_ghosted );

    const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
    EXPECT_FALSE( is_entity_shared );
  }

  // Element 3 has been delete so:
  for (unsigned int nodeId=1;nodeId<=20;nodeId++)
  {
    stk::mesh::EntityKey nodeKey(stk::topology::NODE_RANK,nodeId);
    stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), nodeKey);

    const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(nodeKey, stkMeshBulkData.aura_ghosting()).empty();
    if ( nodeId >=5 and nodeId <=8)
    {
      EXPECT_TRUE( iter != stkMeshBulkData.my_internal_comm_list().end() && iter->key == nodeKey );
      EXPECT_TRUE( is_entity_ghosted );
    }
    else
    {
      EXPECT_FALSE ( iter != stkMeshBulkData.my_internal_comm_list().end() && iter->key == nodeKey ) <<
                                                                                                        "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;

      EXPECT_FALSE( is_entity_ghosted ) <<
                                           "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
    }

    const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(nodeKey, stkMeshBulkData.shared_ghosting()).empty();
    EXPECT_FALSE( is_entity_shared ) <<
                                        "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
  }
}



void connectElementToEdge(stk::unit_test_util::BulkDataTester& stkMeshBulkData, stk::mesh::Entity element,
                          stk::mesh::Entity edge, const std::vector<stk::mesh::EntityId>& nodeIdsForEdge)
{
  std::vector<stk::mesh::Entity> nodes(nodeIdsForEdge.size());
  for (size_t i=0;i<nodeIdsForEdge.size();i++)
  {
    stk::mesh::EntityKey nodeKey(stk::topology::NODE_RANK, nodeIdsForEdge[i]);
    nodes[i] = stkMeshBulkData.get_entity(nodeKey);
  }

  stk::mesh::impl::connectUpwardEntityToEntity(stkMeshBulkData, element, edge, nodes.data());
}

void create_edges(stk::unit_test_util::BulkDataTester& stkMeshBulkData, std::vector<stk::mesh::EntityId>& edgeIds,
                  std::vector<std::vector<stk::mesh::EntityId> > &nodeIdsForEdge,
                  std::vector<std::vector<stk::mesh::EntityId> > &elementRelations,
                  std::vector<stk::mesh::Entity> &edgeEntities,
                  stk::mesh::Part& edge_part)
{
  //============== Create one edge between elements 2 and 3 (nodes 9 and 10) with edge id 100

  stk::mesh::OrdinalVector scratch1, scratch2, scratch3;

  stk::mesh::PartVector add_parts;
  add_parts.push_back( & stkMeshBulkData.mesh_meta_data().get_topology_root_part(stk::topology::LINE_2));
  add_parts.push_back(&edge_part);

  std::vector<bool> communicate_edge_for_ghosting(edgeIds.size(), false);

  for (size_t edge_index=0;edge_index<edgeIds.size();edge_index++)
  {
    stk::mesh::Entity edge = stkMeshBulkData.declare_edge(edgeIds[edge_index], add_parts);
    edgeEntities[edge_index] = edge;

    std::vector<stk::mesh::Entity> ghostedElements(10);
    stk::mesh::Permutation perm = stk::mesh::Permutation::INVALID_PERMUTATION;
    {
      std::vector<stk::mesh::Entity> nodes(2);
      ASSERT_EQ(2u, nodeIdsForEdge[edge_index].size());
      for (size_t n=0; n<nodeIdsForEdge[edge_index].size(); ++n)
      {
        stk::mesh::EntityKey nodeEntityKey(stk::topology::NODE_RANK,nodeIdsForEdge[edge_index][n]);
        stk::mesh::Entity node = stkMeshBulkData.get_entity(nodeEntityKey);
        ASSERT_TRUE(stkMeshBulkData.is_valid(node));
        stkMeshBulkData.declare_relation(edge, node,n, perm, scratch1, scratch2, scratch3);
        EXPECT_TRUE(stkMeshBulkData.bucket(node).member(edge_part));
        nodes[n] = node;
      }

      stk::mesh::Entity const * elemStartNode1 = stkMeshBulkData.begin_elements(nodes[0]);
      stk::mesh::Entity const * elemEndNode1 = stkMeshBulkData.end_elements(nodes[0]);
      stk::mesh::Entity const * elemStartNode2 = stkMeshBulkData.begin_elements(nodes[1]);
      stk::mesh::Entity const * elemEndNode2 = stkMeshBulkData.end_elements(nodes[1]);

      std::vector<stk::mesh::Entity> elems1(elemStartNode1, elemEndNode1);
      std::sort(elems1.begin(), elems1.end());
      std::vector<stk::mesh::Entity> elems2(elemStartNode2, elemEndNode2);
      std::sort(elems2.begin(), elems2.end());

      std::vector<stk::mesh::Entity>::iterator iter = std::set_intersection( elems1.begin(), elems1.end(),
                                                                             elems2.begin(), elems2.end(), ghostedElements.begin());

      ghostedElements.resize(iter-ghostedElements.begin());
    }

    for (size_t j=0;j<ghostedElements.size();j++)
    {
      if ( !stkMeshBulkData.bucket(ghostedElements[j]).owned() )
      {
        connectElementToEdge(stkMeshBulkData, ghostedElements[j], edge, nodeIdsForEdge[edge_index]);
      }
      else if ( stkMeshBulkData.is_ghosted_somewhere(stkMeshBulkData.entity_key(ghostedElements[j])) )
      {
        communicate_edge_for_ghosting[edge_index] = true;
      }
    }

    for (size_t j=0;j<elementRelations[edge_index].size();j++)
    {
      stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,elementRelations[edge_index][j]);
      stk::mesh::Entity element = stkMeshBulkData.get_entity(elementKey);
      connectElementToEdge(stkMeshBulkData, element, edge, nodeIdsForEdge[edge_index]);
    }
  }
}

void checkResultsOfIRSMD_for_edges(stk::unit_test_util::BulkDataTester &stkMeshBulkData)
{
  bool isNodeValidInCommList[2][4] = { { true, true, true, true, },
                                       { true, true, true, true } };

  //  bool isEdgeInCommList[2][1] = { {true}, {false} };
  //

  bool isElementValidInCommList[2][4] = { { false, true, true, false, },
                                          { false, true, true, false } };

  bool isElementInAuraCommMap[2][4] = { { false, true, true, false, },
                                        { false, true, true, false } };

  stk::mesh::EntityCommListInfoVector::const_iterator iter;

  //============== Testing against the boolean vectors

  for (unsigned int node=9;node<=12;node++)
  {
    int nodeoffset = node-9;
    stk::mesh::EntityKey nodeEntityKey(stk::topology::NODE_RANK,node);
    iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), nodeEntityKey);
    if (iter == stkMeshBulkData.my_internal_comm_list().end() || iter->key != nodeEntityKey)
    {
      EXPECT_EQ(isNodeValidInCommList[stkMeshBulkData.parallel_rank()][nodeoffset], false);
    }
    else
    {
      EXPECT_EQ(isNodeValidInCommList[stkMeshBulkData.parallel_rank()][nodeoffset], stkMeshBulkData.is_valid(iter->entity));
    }

    const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(nodeEntityKey, stkMeshBulkData.aura_ghosting()).empty();
    const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(nodeEntityKey, stkMeshBulkData.shared_ghosting()).empty();
    EXPECT_TRUE( is_entity_shared );
    EXPECT_FALSE( is_entity_ghosted );
  }

  for (unsigned int element=1;element<=4;element++)
  {
    int elementoffset = element-1;
    stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
    iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), elementKey);
    if (iter == stkMeshBulkData.my_internal_comm_list().end() || iter->key != elementKey)
    {
      EXPECT_EQ(isElementValidInCommList[stkMeshBulkData.parallel_rank()][elementoffset], false);
    }
    else
    {
      EXPECT_EQ(isElementValidInCommList[stkMeshBulkData.parallel_rank()][elementoffset], stkMeshBulkData.is_valid(iter->entity));
    }

    const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
    EXPECT_EQ(isElementInAuraCommMap[stkMeshBulkData.parallel_rank()][elementoffset], is_entity_ghosted);

    const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
    EXPECT_FALSE( is_entity_shared );
  }
}

void checkResultsOfIRGMD_for_edges(stk::unit_test_util::BulkDataTester &stkMeshBulkData, std::vector<stk::mesh::Entity> &edgeEntities)
{
  for (size_t i=0;i<edgeEntities.size();i++)
  {
    EXPECT_FALSE(stkMeshBulkData.bucket(edgeEntities[i]).shared());
  }

  bool isElementValidInCommListAfterIRGMD[2][4] = { { false, true, false, false, },
                                                    { false, false, true, false } };

  //============== Check results against boolean above

  for (unsigned int element=1;element<=4;element++)
  {
    int elementoffset = element-1;
    stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
    stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), elementKey);
    if (iter == stkMeshBulkData.my_internal_comm_list().end() || iter->key != elementKey)
    {
      EXPECT_EQ(isElementValidInCommListAfterIRGMD[stkMeshBulkData.parallel_rank()][elementoffset], false);
    }
    else
    {
      EXPECT_EQ(isElementValidInCommListAfterIRGMD[stkMeshBulkData.parallel_rank()][elementoffset], stkMeshBulkData.is_valid(iter->entity));
    }

    const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
    EXPECT_FALSE( is_entity_ghosted );

    const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
    EXPECT_FALSE( is_entity_shared );
  }

  //============== Check results of updating the comm list based on changes in the comm map

  stkMeshBulkData.my_update_comm_list_based_on_changes_in_comm_map();

  //============== No element is ghosted, No node is shared, and comm list is empty

  for (unsigned int element=1;element<=4;element++)
  {
    stk::mesh::EntityKey elementKey(stk::topology::ELEMENT_RANK,element);
    stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), elementKey);
    EXPECT_TRUE(iter == stkMeshBulkData.my_internal_comm_list().end());

    const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.aura_ghosting()).empty();
    EXPECT_FALSE( is_entity_ghosted );

    const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(elementKey, stkMeshBulkData.shared_ghosting()).empty();
    EXPECT_FALSE( is_entity_shared );
  }

  bool isNodeInCommList[2][20] = {
    {false, false, false, false,
     true,  true,  true,  true,
     true,  true,  true,  true,
     true,  true,  true,  true,
     false, false, false, false,},

    {false, false, false, false,
     true,  true,  true,  true,
     true,  true,  true,  true,
     true,  true,  true,  true,
     false, false, false, false}
  };

  bool isNodeInGhostedCommMap[2][20] = {
    {false, false, false, false,
     true,  true,  true,  true,
     false, false, false, false,
     true,  true,  true,  true,
     false, false, false, false,
    },
    {false, false, false, false,
     true,  true,  true,  true,
     false, false, false, false,
     true,  true,  true,  true,
     false, false, false, false,
    }
  };

  bool IsNodeInSharedCommMap[2][20] = {
    {false, false, false, false,
     false, false, false, false,
     true,  true,  true,  true,
     false, false, false, false,
     false, false, false, false,
    },
    {false, false, false, false,
     false, false, false, false,
     true,  true,  true,  true,
     false, false, false, false,
     false, false, false, false,
    }
  };

  for (unsigned int nodeId=1;nodeId<=20;nodeId++)
  {
    stk::mesh::EntityKey nodeKey(stk::topology::NODE_RANK,nodeId);
    stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.my_internal_comm_list().begin(), stkMeshBulkData.my_internal_comm_list().end(), nodeKey);

    if ( iter != stkMeshBulkData.my_internal_comm_list().end() && iter->key == nodeKey )
    {
      EXPECT_TRUE(isNodeInCommList[stkMeshBulkData.parallel_rank()][nodeId-1]) <<
                                                                                  "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
    }
    else
    {
      EXPECT_FALSE(isNodeInCommList[stkMeshBulkData.parallel_rank()][nodeId-1]) <<
                                                                                   "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
    }

    const bool is_entity_ghosted = !stkMeshBulkData.my_internal_entity_comm_map(nodeKey, stkMeshBulkData.aura_ghosting()).empty();
    EXPECT_EQ( isNodeInGhostedCommMap[stkMeshBulkData.parallel_rank()][nodeId-1], is_entity_ghosted ) <<
                                                                                                         "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;

    const bool is_entity_shared = !stkMeshBulkData.my_internal_entity_comm_map(nodeKey, stkMeshBulkData.shared_ghosting()).empty();
    EXPECT_EQ( IsNodeInSharedCommMap[stkMeshBulkData.parallel_rank()][nodeId-1], is_entity_shared ) <<
                                                                                                       "Proc " << stkMeshBulkData.parallel_rank() << " fails for node " << nodeId << std::endl;
  }
}

void checkResults(stk::unit_test_util::BulkDataTester& stkMeshBulkData,
                  const size_t numEntities,
                  const stk::mesh::Part& edge_part,
                  const std::vector<stk::mesh::EntityId>& entityIds,
                  stk::mesh::EntityRank rank,
                  bool *isEntityValidOnBulkData,
                  bool *isEntityValidOnCommList,
                  bool *isEntityOwned,
                  bool *isEntityGhosted,
                  bool *isEntityShared,
                  bool *isEntityOnEdgePart,
                  bool *isEntityOnAuraCommMap,
                  bool *isEntityOnSharedCommMap)
{
  int procId=-1;
  MPI_Comm_rank(MPI_COMM_WORLD, &procId);

  for (size_t i=0;i<numEntities;i++)
  {
    stk::mesh::EntityKey entity_key(rank, entityIds[i]);
    stk::mesh::Entity entity = stkMeshBulkData.get_entity(entity_key);
    EXPECT_EQ(isEntityValidOnBulkData[i], stkMeshBulkData.is_valid(entity)) << "P[" << procId << "] for rank: " << rank;
    if ( isEntityValidOnCommList[i] )
    {
      makeSureEntityIsValidOnCommList(stkMeshBulkData, entity_key);
    }
    else
    {
      makeSureEntityIs_NOT_ValidOnCommList(stkMeshBulkData, entity_key);
    }
    if ( isEntityValidOnBulkData[i] )
    {
      ASSERT_TRUE(stkMeshBulkData.bucket_ptr(entity) != 0);
      EXPECT_EQ(isEntityOwned[i], stkMeshBulkData.bucket(entity).owned())<< "P[" << procId << "] for rank: " << rank;
      EXPECT_EQ(isEntityGhosted[i], stkMeshBulkData.bucket(entity).member(stkMeshBulkData.mesh_meta_data().aura_part()))<< "P[" << procId << "] for rank: " << rank << " of index " << entity_key;
      EXPECT_EQ(isEntityShared[i], stkMeshBulkData.bucket(entity).shared())<< "P[" << procId << "] for rank: " << rank;
      EXPECT_EQ(isEntityOnEdgePart[i], stkMeshBulkData.bucket(entity).member(edge_part))<< "P[" << procId << "] for rank: " << rank << " with id " << entityIds[i];
    }
    else
    {
      EXPECT_FALSE(isEntityOwned[i])<< "P[" << procId << "] for rank: " << rank;
      EXPECT_FALSE(isEntityGhosted[i])<< "P[" << procId << "] for rank: " << rank;
      EXPECT_FALSE(isEntityShared[i])<< "P[" << procId << "] for rank: " << rank;
      EXPECT_FALSE(isEntityOnEdgePart[i])<< "P[" << procId << "] for rank: " << rank;
    }
    const bool is_node_on_aura_comm_map = !stkMeshBulkData.my_internal_entity_comm_map(entity_key, stkMeshBulkData.aura_ghosting()).empty();
    EXPECT_EQ( isEntityOnAuraCommMap[i], is_node_on_aura_comm_map )<< "P[" << procId << "] for rank: " << rank;

    const bool is_node_on_shared_comm_map = !stkMeshBulkData.my_internal_entity_comm_map(entity_key, stkMeshBulkData.shared_ghosting()).empty();
    EXPECT_EQ( isEntityOnSharedCommMap[i], is_node_on_shared_comm_map )<< "P[" << procId << "] for rank: " << rank;
  }
}

void checkEntityRelations(int procId, stk::unit_test_util::BulkDataTester& stkMeshBulkData)
{
  int counter=0;
  {
    stk::mesh::EntityId elementConn[24] = {
      1, 2, 4, 3, 5, 6, 8, 7,
      5, 6, 8, 7, 9, 10, 12, 11,
      9, 10, 12, 11, 13, 14, 16, 15
    };

    for (size_t i=0;i<3;i++)
    {
      stk::mesh::EntityKey entity_key(stk::topology::ELEMENT_RANK, i+1+procId);
      stk::mesh::Entity entity = stkMeshBulkData.get_entity(entity_key);
      stk::mesh::Entity const * nodesbegin = stkMeshBulkData.begin_nodes(entity);
      ASSERT_TRUE(nodesbegin!=0) << "for proc " << procId;
      stk::mesh::Entity const * nodesend = stkMeshBulkData.end_nodes(entity);
      for (stk::mesh::Entity const * node = nodesbegin; node != nodesend; ++node)
      {
        stk::mesh::EntityKey node_key(stk::topology::NODE_RANK, elementConn[counter]+4*procId);
        stk::mesh::Entity conn_node = stkMeshBulkData.get_entity(node_key);
        EXPECT_EQ(*node, conn_node) << "for proc " << procId;
        counter++;
      }
    }
  }

  stk::mesh::EntityKey edge_key(stk::topology::EDGE_RANK, 100);
  stk::mesh::Entity edge = stkMeshBulkData.get_entity(edge_key);

  stk::mesh::EntityKey element2_key(stk::topology::ELEMENT_RANK, 2);
  stk::mesh::Entity element2 = stkMeshBulkData.get_entity(element2_key);

  stk::mesh::EntityKey element3_key(stk::topology::ELEMENT_RANK, 3);
  stk::mesh::Entity element3 = stkMeshBulkData.get_entity(element3_key);

  {
    for (size_t i=2;i<=3;i++)
    {
      stk::mesh::EntityKey entity_key(stk::topology::ELEMENT_RANK, i);
      stk::mesh::Entity entity = stkMeshBulkData.get_entity(entity_key);
      stk::mesh::Entity const * edges_begin = stkMeshBulkData.begin_edges(entity);
      ASSERT_TRUE(edges_begin!=0) << "for proc " << procId << " against element " << entity_key;
      EXPECT_EQ( *edges_begin, edge) << "for proc " << procId << " against element " << entity_key;
    }
  }

  stk::mesh::EntityKey node9_key(stk::topology::NODE_RANK, 9);
  stk::mesh::Entity node9 = stkMeshBulkData.get_entity(node9_key);

  stk::mesh::EntityKey node10_key(stk::topology::NODE_RANK, 10);
  stk::mesh::Entity node10 = stkMeshBulkData.get_entity(node10_key);

  {
    stk::mesh::Entity const * entity = 0;
    entity = stkMeshBulkData.begin_nodes(edge);
    ASSERT_TRUE(entity!=0) << "for proc " << procId;

    EXPECT_EQ(*entity, node9) << "for proc " << procId;
    entity++;
    EXPECT_EQ(*entity, node10) << "for proc " << procId;
  }

  {
    stk::mesh::Entity const * elements_begin = stkMeshBulkData.begin_elements(edge);
    ASSERT_TRUE(elements_begin!=0) << "for proc " << procId;

    bool one_of_elements = (*elements_begin == element2) || (*elements_begin == element3);
    EXPECT_TRUE(one_of_elements) << "for proc " << procId;
    elements_begin++;
    one_of_elements = *elements_begin == element2 || *elements_begin == element3;
    EXPECT_TRUE(one_of_elements) << "for proc " << procId;
  }

  {
    stk::mesh::Entity const * entity = 0;

    entity = stkMeshBulkData.begin_edges(node9);
    ASSERT_TRUE(entity!=0) << "for proc " << procId;
    EXPECT_EQ(*entity, edge) << "for proc " << procId;

    entity = stkMeshBulkData.begin_edges(node10);
    ASSERT_TRUE(entity!=0) << "for proc " << procId;
    EXPECT_EQ(*entity, edge) << "for proc " << procId;

    entity = stkMeshBulkData.begin_elements(node9);
    ASSERT_TRUE(entity!=0) << "for proc " << procId;
    ASSERT_EQ(2u, stkMeshBulkData.num_elements(node9));
    bool has_elem2_and_elem3 = (entity[0] == element3 && entity[1] == element2)
        || (entity[0] == element2 && entity[1] == element3);
    EXPECT_TRUE(has_elem2_and_elem3) << "for proc " << procId;

    entity = stkMeshBulkData.begin_elements(node10);
    ASSERT_TRUE(entity!=0) << "for proc " << procId;
    ASSERT_EQ(2u, stkMeshBulkData.num_elements(node9));
    has_elem2_and_elem3 = (entity[0] == element3 && entity[1] == element2)
        || (entity[0] == element2 && entity[1] == element3);
    EXPECT_TRUE(has_elem2_and_elem3) << "for proc " << procId;
  }
}

void check_it_all_for_proc_0(stk::unit_test_util::BulkDataTester &stkMeshBulkData)
{
  size_t numNodes = 20;
  bool isNodeValidOnBulkData[20] = {
    true, true, true, true,
    true, true, true, true,
    true, true, true, true,
    true, true, true, true,
    false, false, false, false
  };
  bool isNodeValidOnCommList[20] = {
    false, false, false, false,
    true, true, true, true,
    true, true, true, true,
    true, true, true, true,
    false, false, false, false
  };
  bool isNodeOwned[20] = {
    true, true, true, true,
    true, true, true, true,
    true, true, true, true,
    false, false, false, false,
    false, false, false, false
  };
  bool isNodeShared[20] = {
    false, false, false, false,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    false, false, false, false
  };
  bool isNodeGhosted[20] = {
    false, false, false, false,
    false, false, false, false,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false
  };
  bool isNodeOnAuraCommMap[20] = {
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false
  };
  bool isNodeOnSharedCommMap[20] = {
    false, false, false, false,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    false, false, false, false
  };
  bool isNodeOnEdgePart[20] = {
    false, false, false, false,
    false, false, false, false,
    true, true, false, false,
    false, false, false, false,
    false, false, false, false
  };
  stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);
  std::vector<stk::mesh::EntityId> nodeIds(numNodes);
  for (size_t i=0;i<nodeIds.size();i++)
  {
    nodeIds[i] = i+1;
  }
  checkResults(stkMeshBulkData, numNodes, edge_part, nodeIds, stk::topology::NODE_RANK, isNodeValidOnBulkData, isNodeValidOnCommList, isNodeOwned,
               isNodeGhosted, isNodeShared, isNodeOnEdgePart, isNodeOnAuraCommMap, isNodeOnSharedCommMap);

  size_t numEdges = 1;
  bool isEdgeValidOnBulkData[1] = {
    true
  };
  bool isEdgeValidOnCommList[1] = {
    true
  };
  bool isEdgeOwned[1] = {
    true
  };
  bool isEdgeShared[1] = {
    true
  };
  bool isEdgeGhosted[1] = {
    false
  };
  bool isEdgeOnAuraCommMap[1] = {
    false
  };
  bool isEdgeOnSharedCommMap[1] = {
    true
  };
  bool isEdgeOnEdgePart[1] = {
    true
  };
  std::vector<stk::mesh::EntityId> edgeIds(numEdges);
  edgeIds[0] = 100;
  checkResults(stkMeshBulkData, numEdges, edge_part, edgeIds, stk::topology::EDGE_RANK, isEdgeValidOnBulkData, isEdgeValidOnCommList, isEdgeOwned,
               isEdgeGhosted, isEdgeShared, isEdgeOnEdgePart, isEdgeOnAuraCommMap, isEdgeOnSharedCommMap);

  size_t numElements = 4;
  bool isElementValidOnBulkData[4] = {
    true, true, true, false
  };
  bool isElementValidOnCommList[4] = {
    false, true, true, false
  };
  bool isElementOwned[4] = {
    true, true, false, false
  };
  bool isElementShared[4] = {
    false, false, false, false
  };
  bool isElementGhosted[4] = {
    false, false, true, false
  };
  bool isElementOnAuraCommMap[4] = {
    false, true, true, false
  };
  bool isElementOnSharedCommMap[4] = {
    false, false, false, false
  };
  bool isElementOnEdgePart[4] = {
    false, false, false, false
  };
  std::vector<stk::mesh::EntityId> elementIds(numElements);
  for (size_t i=0;i<elementIds.size();i++)
  {
    elementIds[i] = i+1;
  }
  checkResults(stkMeshBulkData, numElements, edge_part, elementIds, stk::topology::ELEMENT_RANK, isElementValidOnBulkData, isElementValidOnCommList, isElementOwned,
               isElementGhosted, isElementShared, isElementOnEdgePart, isElementOnAuraCommMap, isElementOnSharedCommMap);

  checkEntityRelations(0, stkMeshBulkData);
}

void check_it_all_for_proc_1(stk::unit_test_util::BulkDataTester &stkMeshBulkData)
{
  size_t numNodes = 20;
  bool isNodeValidOnBulkData[20] = {
    false, false, false, false,
    true, true, true, true,
    true, true, true, true,
    true, true, true, true,
    true, true, true, true
  };

  bool isNodeValidOnCommList[20] = {
    false, false, false, false,
    true, true, true, true,
    true, true, true, true,
    true, true, true, true,
    false, false, false, false
  };

  bool isNodeOwned[20] = {
    false, false, false, false,
    false, false, false, false,
    false, false, false, false,
    true, true, true, true,
    true, true, true, true
  };
  bool isNodeShared[20] = {
    false, false, false, false,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    false, false, false, false
  };
  bool isNodeGhosted[20] = {
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    false, false, false, false,
    false, false, false, false
  };
  bool isNodeOnAuraCommMap[20] = {
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false
  };
  bool isNodeOnSharedCommMap[20] = {
    false, false, false, false,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    false, false, false, false
  };
  bool isNodeOnEdgePart[20] = {
    false, false, false, false,
    false, false, false, false,
    true, true, false, false,
    false, false, false, false,
    false, false, false, false
  };
  stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);
  std::vector<stk::mesh::EntityId> nodeIds(numNodes);
  for (size_t i=0;i<nodeIds.size();i++)
  {
    nodeIds[i] = i+1;
  }
  checkResults(stkMeshBulkData, numNodes, edge_part, nodeIds, stk::topology::NODE_RANK, isNodeValidOnBulkData, isNodeValidOnCommList, isNodeOwned,
               isNodeGhosted, isNodeShared, isNodeOnEdgePart, isNodeOnAuraCommMap, isNodeOnSharedCommMap);

  size_t numEdges = 1;
  bool isEdgeValidOnBulkData[1] = {
    true
  };
  bool isEdgeValidOnCommList[1] = {
    true
  };
  bool isEdgeOwned[1] = {
    false
  };
  bool isEdgeShared[1] = {
    true
  };
  bool isEdgeGhosted[1] = {
    false
  };
  bool isEdgeOnAuraCommMap[1] = {
    false
  };
  bool isEdgeOnSharedCommMap[1] = {
    true
  };
  bool isEdgeOnEdgePart[1] = {
    true
  };
  std::vector<stk::mesh::EntityId> edgeIds(numEdges);
  edgeIds[0] = 100;
  checkResults(stkMeshBulkData, numEdges, edge_part, edgeIds, stk::topology::EDGE_RANK, isEdgeValidOnBulkData, isEdgeValidOnCommList, isEdgeOwned,
               isEdgeGhosted, isEdgeShared, isEdgeOnEdgePart, isEdgeOnAuraCommMap, isEdgeOnSharedCommMap);

  size_t numElements = 4;
  bool isElementValidOnBulkData[4] = {
    false, true, true, true
  };
  bool isElementValidOnCommList[4] = {
    false, true, true, false
  };
  bool isElementOwned[4] = {
    false, false, true, true
  };
  bool isElementShared[4] = {
    false, false, false, false
  };
  bool isElementGhosted[4] = {
    false, true, false, false
  };
  bool isElementOnAuraCommMap[4] = {
    false, true, true, false
  };
  bool isElementOnSharedCommMap[4] = {
    false, false, false, false
  };
  bool isElementOnEdgePart[4] = {
    false, false, false, false
  };
  std::vector<stk::mesh::EntityId> elementIds(numElements);
  for (size_t i=0;i<elementIds.size();i++)
  {
    elementIds[i] = i+1;
  }
  checkResults(stkMeshBulkData, numElements, edge_part, elementIds, stk::topology::ELEMENT_RANK, isElementValidOnBulkData, isElementValidOnCommList, isElementOwned,
               isElementGhosted, isElementShared, isElementOnEdgePart, isElementOnAuraCommMap, isElementOnSharedCommMap);

  checkEntityRelations(1, stkMeshBulkData);
}

void checkItAllForThisCase(stk::unit_test_util::BulkDataTester &stkMeshBulkData)
{
  checkThatMeshIsParallelConsistent(stkMeshBulkData);
  std::vector<size_t> globalCounts;
  stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
  EXPECT_EQ(20u, globalCounts[stk::topology::NODE_RANK]);
  EXPECT_EQ(1u, globalCounts[stk::topology::EDGE_RANK]);
  EXPECT_EQ(0u, globalCounts[stk::topology::FACE_RANK]);
  EXPECT_EQ(4u, globalCounts[stk::topology::ELEMENT_RANK]);

  if ( stkMeshBulkData.parallel_rank() == 0)
  {
    check_it_all_for_proc_0(stkMeshBulkData);
  }
  else
  {
    check_it_all_for_proc_1(stkMeshBulkData);
  }
}

void checkEntityRelationsGhosted(int procId, stk::mesh::BulkData& stkMeshBulkData)
{
  int counter=0;
  {
    stk::mesh::EntityId elementConn[24] = {
      1, 2, 4, 3, 5, 6, 8, 7,
      5, 6, 8, 7, 9, 10, 12, 11,
      9, 10, 12, 11, 13, 14, 16, 15
    };

    for (size_t i=0;i<3;i++)
    {
      stk::mesh::EntityKey entity_key(stk::topology::ELEMENT_RANK, i+1+procId);
      stk::mesh::Entity entity = stkMeshBulkData.get_entity(entity_key);
      stk::mesh::Entity const * nodesbegin = stkMeshBulkData.begin_nodes(entity);
      ASSERT_TRUE(nodesbegin!=0) << "for proc " << procId;
      stk::mesh::Entity const * nodesend = stkMeshBulkData.end_nodes(entity);
      for (stk::mesh::Entity const * node = nodesbegin; node != nodesend; ++node)
      {
        stk::mesh::EntityKey node_key(stk::topology::NODE_RANK, elementConn[counter]+4*procId);
        stk::mesh::Entity conn_node = stkMeshBulkData.get_entity(node_key);
        EXPECT_EQ(*node, conn_node) << "for proc " << procId;
        counter++;
      }
    }
  }

  stk::mesh::EntityKey edge_key1(stk::topology::EDGE_RANK, 100);
  stk::mesh::Entity edge1 = stkMeshBulkData.get_entity(edge_key1);

  stk::mesh::EntityKey edge_key2(stk::topology::EDGE_RANK, 101);
  stk::mesh::Entity edge2 = stkMeshBulkData.get_entity(edge_key2);

  std::vector<std::pair<int, stk::mesh::Entity> > element2edge;
  if ( procId == 0 )
  {
    element2edge.emplace_back(1, edge1);
    element2edge.emplace_back(2, edge1);
    element2edge.emplace_back(3, edge2);
  }
  else
  {
    element2edge.emplace_back(2, edge1);
    element2edge.emplace_back(3, edge2);
    element2edge.emplace_back(4, edge2);
  }

  {
    for (size_t i=0;i<element2edge.size();i++)
    {
      stk::mesh::EntityKey entity_key(stk::topology::ELEMENT_RANK, element2edge[i].first);
      stk::mesh::Entity entity = stkMeshBulkData.get_entity(entity_key);
      stk::mesh::Entity const * edges_begin = stkMeshBulkData.begin_edges(entity);
      ASSERT_TRUE(edges_begin!=0) << "for proc " << procId << " against element " << entity_key;
      EXPECT_EQ( *edges_begin, element2edge[i].second ) << "for proc " << procId << " against element " << entity_key;
    }
  }

  std::vector<int> edge_node_ids(2);
  if ( procId == 0 )
  {
    edge_node_ids[0] = 5;
    edge_node_ids[1] = 6;
  }
  else
  {
    edge_node_ids[0] = 13;
    edge_node_ids[1] = 14;
  }

  stk::mesh::EntityKey node1_key(stk::topology::NODE_RANK, edge_node_ids[0]);
  stk::mesh::Entity node1 = stkMeshBulkData.get_entity(node1_key);

  stk::mesh::EntityKey node2_key(stk::topology::NODE_RANK, edge_node_ids[1]);
  stk::mesh::Entity node2 = stkMeshBulkData.get_entity(node2_key);

  stk::mesh::Entity edgeToLocalNodes = edge1;
  if ( procId == 1)
  {
    edgeToLocalNodes = edge2;
  }

  {
    stk::mesh::Entity const * entity = 0;
    entity = stkMeshBulkData.begin_nodes(edgeToLocalNodes);
    ASSERT_TRUE(entity!=0) << "for proc " << procId;

    EXPECT_EQ(*entity, node1) << "for proc " << procId;
    entity++;
    EXPECT_EQ(*entity, node2) << "for proc " << procId;
  }

  std::vector<int> conn_elem_ids(2);
  if ( procId == 0 )
  {
    conn_elem_ids[0] = 1;
    conn_elem_ids[1] = 2;
  }
  else
  {
    conn_elem_ids[0] = 3;
    conn_elem_ids[1] = 4;
  }

  stk::mesh::EntityKey elementA_key(stk::topology::ELEMENT_RANK, conn_elem_ids[0]);
  stk::mesh::Entity first_element = stkMeshBulkData.get_entity(elementA_key);

  stk::mesh::EntityKey elementB_key(stk::topology::ELEMENT_RANK, conn_elem_ids[1]);
  stk::mesh::Entity second_element = stkMeshBulkData.get_entity(elementB_key);

  {
    stk::mesh::Entity const * elements_begin = stkMeshBulkData.begin_elements(edgeToLocalNodes);
    ASSERT_TRUE(elements_begin!=0) << "for proc " << procId;

    bool one_of_elements = (*elements_begin == first_element) || (*elements_begin == second_element);
    EXPECT_TRUE(one_of_elements) << "for proc " << procId;
    elements_begin++;
    one_of_elements = (*elements_begin == first_element) || (*elements_begin == second_element);
    EXPECT_TRUE(one_of_elements) << "for proc " << procId;
  }

  {
    stk::mesh::Entity const * entity = 0;

    entity = stkMeshBulkData.begin_edges(node1);
    ASSERT_TRUE(entity!=0) << "for proc " << procId;
    EXPECT_EQ(*entity, edgeToLocalNodes) << "for proc " << procId;

    entity = stkMeshBulkData.begin_edges(node2);
    ASSERT_TRUE(entity!=0) << "for proc " << procId;
    EXPECT_EQ(*entity, edgeToLocalNodes) << "for proc " << procId;

    entity = stkMeshBulkData.begin_elements(node1);
    ASSERT_TRUE(entity!=0) << "for proc " << procId;
    bool connected_to_valid_element = *entity == second_element || *entity == first_element;
    EXPECT_TRUE(connected_to_valid_element) << "for proc " << procId;
    entity++;
    connected_to_valid_element = *entity == second_element || *entity == first_element;
    EXPECT_TRUE(connected_to_valid_element) << "for proc " << procId;

    entity = stkMeshBulkData.begin_elements(node2);
    ASSERT_TRUE(entity!=0) << "for proc " << procId;
    connected_to_valid_element = *entity == second_element || *entity == first_element;
    EXPECT_TRUE(connected_to_valid_element) << "for proc " << procId;
    entity++;
    connected_to_valid_element = *entity == second_element || *entity == first_element;
    EXPECT_TRUE(connected_to_valid_element) << "for proc " << procId;
  }
}

void check_it_all_for_proc_0_ghosted(stk::unit_test_util::BulkDataTester &stkMeshBulkData)
{
  size_t numNodes = 20;
  bool isNodeValidOnBulkData[20] = {
    true, true, true, true,
    true, true, true, true,
    true, true, true, true,
    true, true, true, true,
    false, false, false, false
  };
  bool isNodeValidOnCommList[20] = {
    false, false, false, false,
    true, true, true, true,
    true, true, true, true,
    true, true, true, true,
    false, false, false, false
  };
  bool isNodeOwned[20] = {
    true, true, true, true,
    true, true, true, true,
    true, true, true, true,
    false, false, false, false,
    false, false, false, false
  };
  bool isNodeShared[20] = {
    false, false, false, false,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    false, false, false, false
  };
  bool isNodeGhosted[20] = {
    false, false, false, false,
    false, false, false, false,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false
  };
  bool isNodeOnAuraCommMap[20] = {
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false
  };
  bool isNodeOnSharedCommMap[20] = {
    false, false, false, false,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    false, false, false, false
  };
  bool isNodeOnEdgePart[20] = {
    false, false, false, false,
    true, true, false, false,
    false, false, false, false,
    true, true, false, false,
    false, false, false, false
  };
  stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);
  std::vector<stk::mesh::EntityId> nodeIds(numNodes);
  for (size_t i=0;i<nodeIds.size();i++)
  {
    nodeIds[i] = i+1;
  }
  checkResults(stkMeshBulkData, numNodes, edge_part, nodeIds, stk::topology::NODE_RANK, isNodeValidOnBulkData, isNodeValidOnCommList, isNodeOwned,
               isNodeGhosted, isNodeShared, isNodeOnEdgePart, isNodeOnAuraCommMap, isNodeOnSharedCommMap);

  size_t numEdges = 2;
  bool isEdgeValidOnBulkData[2] = {
    true, true
  };
  bool isEdgeValidOnCommList[2] = {
    true, true
  };
  bool isEdgeOwned[2] = {
    true, false
  };
  bool isEdgeShared[2] = {
    false, false
  };
  bool isEdgeGhosted[2] = {
    false, true
  };
  bool isEdgeOnAuraCommMap[2] = {
    true, true
  };
  bool isEdgeOnSharedCommMap[2] = {
    false, false
  };
  bool isEdgeOnEdgePart[2] = {
    true, true
  };
  std::vector<stk::mesh::EntityId> edgeIds(numEdges);
  edgeIds[0] = 100;
  edgeIds[1] = 101;

  checkResults(stkMeshBulkData, numEdges, edge_part, edgeIds, stk::topology::EDGE_RANK, isEdgeValidOnBulkData, isEdgeValidOnCommList, isEdgeOwned,
               isEdgeGhosted, isEdgeShared, isEdgeOnEdgePart, isEdgeOnAuraCommMap, isEdgeOnSharedCommMap);

  size_t numElements = 4;
  bool isElementValidOnBulkData[4] = {
    true, true, true, false
  };
  bool isElementValidOnCommList[4] = {
    false, true, true, false
  };
  bool isElementOwned[4] = {
    true, true, false, false
  };
  bool isElementShared[4] = {
    false, false, false, false
  };
  bool isElementGhosted[4] = {
    false, false, true, false
  };
  bool isElementOnAuraCommMap[4] = {
    false, true, true, false
  };
  bool isElementOnSharedCommMap[4] = {
    false, false, false, false
  };
  bool isElementOnEdgePart[4] = {
    false, false, false, false
  };
  std::vector<stk::mesh::EntityId> elementIds(numElements);
  for (size_t i=0;i<elementIds.size();i++)
  {
    elementIds[i] = i+1;
  }
  checkResults(stkMeshBulkData, numElements, edge_part, elementIds, stk::topology::ELEMENT_RANK, isElementValidOnBulkData, isElementValidOnCommList, isElementOwned,
               isElementGhosted, isElementShared, isElementOnEdgePart, isElementOnAuraCommMap, isElementOnSharedCommMap);

  checkEntityRelationsGhosted(0, stkMeshBulkData);
}

void check_it_all_for_proc_1_ghosted(stk::unit_test_util::BulkDataTester &stkMeshBulkData)
{
  size_t numNodes = 20;
  bool isNodeValidOnBulkData[20] = {
    false, false, false, false,
    true, true, true, true,
    true, true, true, true,
    true, true, true, true,
    true, true, true, true
  };

  bool isNodeValidOnCommList[20] = {
    false, false, false, false,
    true, true, true, true,
    true, true, true, true,
    true, true, true, true,
    false, false, false, false
  };

  bool isNodeOwned[20] = {
    false, false, false, false,
    false, false, false, false,
    false, false, false, false,
    true, true, true, true,
    true, true, true, true
  };
  bool isNodeShared[20] = {
    false, false, false, false,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    false, false, false, false
  };
  bool isNodeGhosted[20] = {
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    false, false, false, false,
    false, false, false, false
  };
  bool isNodeOnAuraCommMap[20] = {
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false
  };
  bool isNodeOnSharedCommMap[20] = {
    false, false, false, false,
    false, false, false, false,
    true, true, true, true,
    false, false, false, false,
    false, false, false, false
  };
  bool isNodeOnEdgePart[20] = {
    false, false, false, false,
    true, true, false, false,
    false, false, false, false,
    true, true, false, false,
    false, false, false, false
  };
  stk::mesh::Part& edge_part = stkMeshBulkData.mesh_meta_data().declare_part("edge_part", stk::topology::EDGE_RANK);
  std::vector<stk::mesh::EntityId> nodeIds(numNodes);
  for (size_t i=0;i<nodeIds.size();i++)
  {
    nodeIds[i] = i+1;
  }
  checkResults(stkMeshBulkData, numNodes, edge_part, nodeIds, stk::topology::NODE_RANK, isNodeValidOnBulkData, isNodeValidOnCommList, isNodeOwned,
               isNodeGhosted, isNodeShared, isNodeOnEdgePart, isNodeOnAuraCommMap, isNodeOnSharedCommMap);

  size_t numEdges = 2;
  bool isEdgeValidOnBulkData[2] = {
    true, true
  };
  bool isEdgeValidOnCommList[2] = {
    true, true
  };
  bool isEdgeOwned[2] = {
    false, true
  };
  bool isEdgeShared[2] = {
    false, false
  };
  bool isEdgeGhosted[2] = {
    true, false
  };
  bool isEdgeOnAuraCommMap[2] = {
    true, true
  };
  bool isEdgeOnSharedCommMap[2] = {
    false, false
  };
  bool isEdgeOnEdgePart[2] = {
    true, true
  };
  std::vector<stk::mesh::EntityId> edgeIds(numEdges);
  edgeIds[0] = 100;
  edgeIds[1] = 101;
  checkResults(stkMeshBulkData, numEdges, edge_part, edgeIds, stk::topology::EDGE_RANK, isEdgeValidOnBulkData, isEdgeValidOnCommList, isEdgeOwned,
               isEdgeGhosted, isEdgeShared, isEdgeOnEdgePart, isEdgeOnAuraCommMap, isEdgeOnSharedCommMap);

  size_t numElements = 4;
  bool isElementValidOnBulkData[4] = {
    false, true, true, true
  };
  bool isElementValidOnCommList[4] = {
    false, true, true, false
  };
  bool isElementOwned[4] = {
    false, false, true, true
  };
  bool isElementShared[4] = {
    false, false, false, false
  };
  bool isElementGhosted[4] = {
    false, true, false, false
  };
  bool isElementOnAuraCommMap[4] = {
    false, true, true, false
  };
  bool isElementOnSharedCommMap[4] = {
    false, false, false, false
  };
  bool isElementOnEdgePart[4] = {
    false, false, false, false
  };
  std::vector<stk::mesh::EntityId> elementIds(numElements);
  for (size_t i=0;i<elementIds.size();i++)
  {
    elementIds[i] = i+1;
  }
  checkResults(stkMeshBulkData, numElements, edge_part, elementIds, stk::topology::ELEMENT_RANK, isElementValidOnBulkData, isElementValidOnCommList, isElementOwned,
               isElementGhosted, isElementShared, isElementOnEdgePart, isElementOnAuraCommMap, isElementOnSharedCommMap);

  checkEntityRelationsGhosted(1, stkMeshBulkData);
}

void checkItAllForThisGhostedCase(stk::unit_test_util::BulkDataTester &stkMeshBulkData)
{
  checkThatMeshIsParallelConsistent(stkMeshBulkData);
  std::vector<size_t> globalCounts;
  stk::mesh::comm_mesh_counts(stkMeshBulkData, globalCounts);
  EXPECT_EQ(20u, globalCounts[stk::topology::NODE_RANK]);
  EXPECT_EQ(2u, globalCounts[stk::topology::EDGE_RANK]);
  EXPECT_EQ(0u, globalCounts[stk::topology::FACE_RANK]);
  EXPECT_EQ(4u, globalCounts[stk::topology::ELEMENT_RANK]);

  if ( stkMeshBulkData.parallel_rank() == 0)
  {
    check_it_all_for_proc_0_ghosted(stkMeshBulkData);
  }
  else
  {
    check_it_all_for_proc_1_ghosted(stkMeshBulkData);
  }
}

} } } // namespace stk mesh unit_test

#endif
