// Copyright (c) 2013, Sandia Corporation.
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
//     * Neither the name of Sandia Corporation nor the names of its
//       contributors may be used to endorse or promote products derived
//       from this software without specific prior written permission.
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

#include <stddef.h>                     // for size_t
#include <stdlib.h>                     // for exit
#include <exception>                    // for exception
#include <iostream>                     // for ostringstream, etc
#include <iterator>                     // for distance
#include <map>                          // for _Rb_tree_const_iterator, etc
#include <stdexcept>                    // for logic_error, runtime_error
#include <stk_mesh/base/BulkData.hpp>   // for BulkData, etc
#include <stk_mesh/base/FieldParallel.hpp>  // for communicate_field_data, etc
#include <stk_mesh/base/GetEntities.hpp>  // for count_entities, etc
#include <stk_mesh/fixtures/BoxFixture.hpp>  // for BoxFixture
#include <stk_mesh/fixtures/HexFixture.hpp>  // for HexFixture, etc
#include <stk_mesh/fixtures/QuadFixture.hpp>  // for QuadFixture
#include <stk_mesh/fixtures/RingFixture.hpp>  // for RingFixture
#include <stk_util/parallel/Parallel.hpp>  // for ParallelMachine, etc
#include <stk_util/parallel/ParallelReduce.hpp>  // for Reduce, ReduceSum, etc
#include <gtest/gtest.h>
#include <string>                       // for string, basic_string, etc
#include <unit_tests/UnitTestModificationEndWrapper.hpp>
#include <unit_tests/UnitTestRingFixture.hpp>  // for test_shift_ring
#include <unit_tests/Setup8Quad4ProcMesh.hpp>
#include <utility>                      // for pair
#include <vector>                       // for vector, etc
#include "stk_mesh/base/Bucket.hpp"     // for Bucket, has_superset
#include "stk_mesh/base/Entity.hpp"     // for Entity
#include "stk_mesh/base/EntityKey.hpp"  // for EntityKey
#include "stk_mesh/base/Field.hpp"      // for Field
#include "stk_mesh/base/FieldBase.hpp"  // for field_data, etc
#include "stk_mesh/base/Ghosting.hpp"   // for Ghosting
#include "stk_mesh/base/MetaData.hpp"   // for MetaData, entity_rank_names, etc
#include "stk_mesh/base/Part.hpp"       // for Part
#include "stk_mesh/base/Relation.hpp"
#include "stk_mesh/base/Selector.hpp"   // for Selector, operator|
#include "stk_mesh/base/Types.hpp"      // for EntityProc, EntityVector, etc
#include "stk_topology/topology.hpp"    // for topology, etc
#include "stk_util/util/PairIter.hpp"   // for PairIter
#include "stk_io/StkMeshIoBroker.hpp"
#include <stk_mesh/base/Comm.hpp>

namespace stk
{
namespace mesh
{
class FieldBase;
}
}

using stk::mesh::Part;
using stk::mesh::MetaData;
using stk::mesh::BulkData;
using stk::mesh::Selector;
using stk::mesh::PartVector;
using stk::mesh::BaseEntityRank;
using stk::mesh::PairIterRelation;
using stk::mesh::EntityProc;
using stk::mesh::Entity;
using stk::mesh::EntityId;
using stk::mesh::EntityKey;
using stk::mesh::EntityVector;
using stk::mesh::EntityRank;
using stk::mesh::fixtures::RingFixture;
using stk::mesh::fixtures::BoxFixture;

namespace
{

class FieldMgr
{
public:
    FieldMgr(stk::mesh::MetaData &stkMeshMetaData) :
        m_stkMeshMetaData(stkMeshMetaData),
        m_commListNodeFieldName("CommListNode"),
        m_commListNodeField(NULL),
        m_sharingCommMapNodeFieldName("SharingCommMapNode"),
        m_sharingCommMapNodeField(NULL),
        m_auraCommMapNodeFieldName("AuraCommMapNode"),
        m_auraCommMapNodeField(NULL),
        m_auraCommMapElementFieldName("AuraCommMapElement"),
        m_auraCommMapElementField(NULL)
    { }
    ~FieldMgr() {}

    void addCommListNodeField()
    {
        stk::mesh::Field<double> &field = m_stkMeshMetaData.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK,
                m_commListNodeFieldName, 1);
        m_commListNodeField = &field;
        stk::mesh::put_field(*m_commListNodeField, m_stkMeshMetaData.universal_part());
    }

    void addSharingCommMapNodeField()
    {
        stk::mesh::Field<double> &field = m_stkMeshMetaData.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK,
                m_sharingCommMapNodeFieldName, 1);
        m_sharingCommMapNodeField = &field;
        stk::mesh::put_field(*m_sharingCommMapNodeField, m_stkMeshMetaData.universal_part());
    }

    void addAuraCommMapNodeField()
    {
        stk::mesh::Field<double> &field = m_stkMeshMetaData.declare_field<stk::mesh::Field<double> >(stk::topology::NODE_RANK,
                m_auraCommMapNodeFieldName, 1);
        m_auraCommMapNodeField = &field;
        stk::mesh::put_field(*m_auraCommMapNodeField, m_stkMeshMetaData.universal_part());
    }

    void addAuraCommMapElementField()
    {
        stk::mesh::Field<double> &field = m_stkMeshMetaData.declare_field<stk::mesh::Field<double> >(stk::topology::ELEMENT_RANK,
                m_auraCommMapElementFieldName, 1);
        m_auraCommMapElementField = &field;
        stk::mesh::put_field(*m_auraCommMapElementField, m_stkMeshMetaData.universal_part());
    }

    void storeFieldPointers()
    {
        m_commListNodeField = static_cast<stk::mesh::Field<double> *>(m_stkMeshMetaData.get_field(stk::topology::NODE_RANK, m_commListNodeFieldName));
        m_sharingCommMapNodeField = static_cast<stk::mesh::Field<double> *>(m_stkMeshMetaData.get_field(stk::topology::NODE_RANK, m_sharingCommMapNodeFieldName));
        m_auraCommMapNodeField = static_cast<stk::mesh::Field<double > *>(m_stkMeshMetaData.get_field(stk::topology::NODE_RANK, m_auraCommMapNodeFieldName));
        m_auraCommMapElementField = static_cast<stk::mesh::Field<double > *>(m_stkMeshMetaData.get_field(stk::topology::ELEMENT_RANK, m_auraCommMapElementFieldName));
    }

    double getFieldValue(stk::mesh::Field<double>* field, stk::mesh::Entity entity)
    {
        double *fieldValue = stk::mesh::field_data(*field, entity);
        return *fieldValue;
    }

    stk::mesh::Field<double> * getCommListNodeField() { return m_commListNodeField; }
    stk::mesh::Field<double> * getSharingCommMapNodeField() { return m_sharingCommMapNodeField; }
    stk::mesh::Field<double> * getAuraCommMapNodeField() { return m_auraCommMapNodeField; }
    stk::mesh::Field<double> * getAuraCommMapElementField() { return m_auraCommMapElementField; }

private:
    stk::mesh::MetaData &m_stkMeshMetaData;

    std::string m_commListNodeFieldName;
    stk::mesh::Field<double> *m_commListNodeField;

    std::string m_sharingCommMapNodeFieldName;
    stk::mesh::Field<double> *m_sharingCommMapNodeField;

    std::string m_auraCommMapNodeFieldName;
    stk::mesh::Field<double> *m_auraCommMapNodeField;

    std::string m_auraCommMapElementFieldName;
    stk::mesh::Field<double> *m_auraCommMapElementField;
};

void testSubMesh(stk::mesh::BulkData &stkMeshBulkData, stk::mesh::Selector select, int elementToTestId, size_t goldNumberNodes, size_t goldNumberElements, FieldMgr &parallelFieldMgr);
void createSerialSubMesh(const stk::mesh::MetaData &meta, stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Selector subMeshSelector, stk::mesh::MetaData &newMeta, stk::mesh::BulkData &newBulkData);
void add_nodes_to_move(stk::mesh::BulkData& bulk,
                       stk::mesh::Entity elem,
                       int dest_proc,
                       std::vector<stk::mesh::EntityProc>& entities_to_move);
bool isEntityValidOnCommList(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity);
bool isEntityInGhostingCommMap(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity);
bool isEntityInSharingCommMap(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity);
void checkCommMaps(std::string message, stk::mesh::BulkData &stkMeshBulkData, int numElements, bool ownerOfElement[], bool isElementInAuraCommMap[], bool isElementValid[],
            int numNodes, bool ownerOfNode[], bool isNodeInSharedCommMap[], bool isNodeInAuraCommMap[], bool isNodeValid[]);
void putCommInfoDataOnFields(stk::mesh::BulkData &bulkData, FieldMgr &fieldMgr);

void writeCommInfoFields(stk::mesh::BulkData &bulkData, FieldMgr &fieldMgr, const std::string& filename, double time);

//////////////////////////////////////////////////
void moveElements2And3ToProc2(stk::mesh::BulkData &stkMeshBulkData);
void runProc0(stk::mesh::BulkData &stkMeshBulkData);
void runProc1(stk::mesh::BulkData &stkMeshBulkData);
void runProc2(stk::mesh::BulkData &stkMeshBulkData);


//////////////////////////////////////////////////////////////////////////////////////////////////////////////
TEST(UnitTestChangeEntityOwner, changeEntityOwnerCase1)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(comm, &numProcs);
    if(numProcs == 3)
    {
        std::string exodusFileName = "generated:1x1x6";
        const int spatialDim = 3;
        stk::mesh::MetaData stkMeshMetaData(spatialDim);
        stk::mesh::BulkData stkMeshBulkData(stkMeshMetaData, comm);

        stk::io::StkMeshIoBroker exodusFileReader(comm);

        exodusFileReader.set_bulk_data(stkMeshBulkData);
        exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);
        exodusFileReader.create_input_mesh();

        FieldMgr fieldMgr(stkMeshMetaData);
        fieldMgr.addCommListNodeField();
        fieldMgr.addSharingCommMapNodeField();
        fieldMgr.addAuraCommMapNodeField();
        fieldMgr.addAuraCommMapElementField();

        exodusFileReader.populate_bulk_data();

        putCommInfoDataOnFields(stkMeshBulkData, fieldMgr);

        double time = 0.5;
        writeCommInfoFields(stkMeshBulkData, fieldMgr, "testBefore.exo",  time);

        {
            stk::mesh::MetaData newMetaData(3);
            stk::mesh::BulkData nBulkData(newMetaData, MPI_COMM_SELF);

            createSerialSubMesh(stkMeshMetaData, stkMeshBulkData, stkMeshBulkData.mesh_meta_data().universal_part(), newMetaData, nBulkData);

            FieldMgr fieldMgrSubMesh(newMetaData);
            fieldMgrSubMesh.storeFieldPointers();

            std::ostringstream oss;
            oss << "testSingleBefore_" << stkMeshBulkData.parallel_rank() << ".exo";
            std::string filename = oss.str();

            writeCommInfoFields(nBulkData, fieldMgrSubMesh, filename,  time);

            unlink(filename.c_str());
        }

        moveElements2And3ToProc2(stkMeshBulkData);

        putCommInfoDataOnFields(stkMeshBulkData, fieldMgr);

        writeCommInfoFields(stkMeshBulkData, fieldMgr, "testAfter.exo",  time);

        MPI_Barrier(MPI_COMM_WORLD);

        {
            stk::mesh::MetaData newMetaData(3);
            stk::mesh::BulkData nBulkData(newMetaData, MPI_COMM_SELF);

            createSerialSubMesh(stkMeshMetaData, stkMeshBulkData, stkMeshBulkData.mesh_meta_data().universal_part(), newMetaData, nBulkData);

            FieldMgr fieldMgrSubMesh(newMetaData);
            fieldMgrSubMesh.storeFieldPointers();

            std::ostringstream oss;
            oss << "testSingleAfter_" << stkMeshBulkData.parallel_rank() << ".exo";
            std::string filename = oss.str();

            writeCommInfoFields(nBulkData, fieldMgrSubMesh, filename,  time);

            unlink(filename.c_str());
        }

        if(stkMeshBulkData.parallel_rank() == 0)
        {
            unlink("testBefore.exo.3.0");
            unlink("testBefore.exo.3.1");
            unlink("testBefore.exo.3.2");
            unlink("testAfter.exo.3.0");
            unlink("testAfter.exo.3.1");
            unlink("testAfter.exo.3.2");
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

TEST(UnitTestChangeEntityOwner, testCreateSubMesh)
{
    MPI_Comm comm = MPI_COMM_WORLD;
    int numProcs = -1;
    MPI_Comm_size(comm, &numProcs);
    if(numProcs == 2)
    {
        std::string exodusFileName = "generated:1x1x4";
        const int spatialDim = 3;
        stk::mesh::MetaData stkMeshMetaData(spatialDim);
        stk::mesh::BulkData stkMeshBulkData(stkMeshMetaData, comm);

        stk::io::StkMeshIoBroker exodusFileReader(comm);

        exodusFileReader.set_bulk_data(stkMeshBulkData);
        exodusFileReader.add_mesh_database(exodusFileName, stk::io::READ_MESH);
        exodusFileReader.create_input_mesh();

        FieldMgr fieldMgr(stkMeshMetaData);
        fieldMgr.addCommListNodeField();
        fieldMgr.addSharingCommMapNodeField();
        fieldMgr.addAuraCommMapNodeField();
        fieldMgr.addAuraCommMapElementField();

        exodusFileReader.populate_bulk_data();
        putCommInfoDataOnFields(stkMeshBulkData, fieldMgr);

        {
            stk::mesh::Selector sel = stkMeshMetaData.universal_part();
            int elementToTestId = 3;
            size_t goldNumberNodes = 16;
            size_t goldNumberElements = 3;

            testSubMesh(stkMeshBulkData, sel, elementToTestId,  goldNumberNodes, goldNumberElements, fieldMgr);
        }

        {
            stk::mesh::Selector sel = stkMeshMetaData.locally_owned_part() | stkMeshMetaData.globally_shared_part();
            int elementToTestId = 2 + stkMeshBulkData.parallel_rank();
            size_t goldNumberNodes = 12;
            size_t goldNumberElements = 2;

            testSubMesh(stkMeshBulkData, sel, elementToTestId,  goldNumberNodes, goldNumberElements, fieldMgr);
        }

        {
            stk::mesh::Part *part = stkMeshMetaData.get_part("block_1");
            stk::mesh::Selector sel = *part;
            int elementToTestId = 2 + stkMeshBulkData.parallel_rank();
            size_t goldNumberNodes = 16;
            size_t goldNumberElements = 3;

            testSubMesh(stkMeshBulkData, sel, elementToTestId,  goldNumberNodes, goldNumberElements, fieldMgr);
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void writeCommInfoFields(stk::mesh::BulkData &bulkData, FieldMgr &fieldMgr, const std::string& filename, double time)
{
    stk::io::StkMeshIoBroker ioBroker(bulkData.parallel());
    ioBroker.set_bulk_data(bulkData);

    int indexAfter = ioBroker.create_output_mesh(filename, stk::io::WRITE_RESULTS);
    ioBroker.add_field(indexAfter, *fieldMgr.getCommListNodeField());
    ioBroker.add_field(indexAfter, *fieldMgr.getSharingCommMapNodeField());
    ioBroker.add_field(indexAfter, *fieldMgr.getAuraCommMapNodeField());
    ioBroker.add_field(indexAfter, *fieldMgr.getAuraCommMapElementField());
    ioBroker.write_output_mesh(indexAfter);
    ioBroker.begin_output_step(indexAfter, time);
    ioBroker.write_defined_output_fields(indexAfter);
    ioBroker.end_output_step(indexAfter);
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void runProc0(stk::mesh::BulkData &stkMeshBulkData)
{
    {
        int numElements = 6;
        bool ownerOfElement[] = { true, true, false, false, false, false };
        bool isElementInAuraCommMap[] = { false, true, true, false, false, false };
        bool isElementValid[] = { true, true, true, false, false, false };

        int numNodes = 28;
        bool ownerOfNode[] = {
                true, true, true, true,
                true, true, true, true,
                true, true, true, true,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false
        };

        bool isNodeInSharedCommMap[] = {
                false, false, false, false,
                false, false, false, false,
                true, true, true, true,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false
        };

        bool isNodeInAuraCommMap[] = {
                false, false, false, false,
                true, true, true, true,
                false, false, false, false,
                true, true, true, true,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false
        };

        bool isNodeValid[] = {
               true, true, true, true,
               true, true, true, true,
               true, true, true, true,
               true, true, true, true,
               false, false, false, false,
               false, false, false, false,
               false, false, false, false
        };

        checkCommMaps("Before on Proc 0", stkMeshBulkData, numElements, ownerOfElement, isElementInAuraCommMap, isElementValid,
                numNodes, ownerOfNode, isNodeInSharedCommMap, isNodeInAuraCommMap, isNodeValid);
    }

    stk::mesh::Entity entity = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, 2);
    std::vector<EntityProc> entitiesToMove;
    int proc2=2;
    entitiesToMove.push_back(std::make_pair(entity, proc2));

    add_nodes_to_move(stkMeshBulkData, entity, proc2, entitiesToMove);

    stkMeshBulkData.change_entity_owner(entitiesToMove);

    {
        int numElements = 6;
        bool ownerOfElement[] = { true, false, false, false, false, false };
        bool isElementInAuraCommMap[] = { true, true, false, false, false, false };
        bool isElementValid[] = { true, true, false, false, false, false };

        int numNodes = 28;
        bool ownerOfNode[] = {
                true, true, true, true,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false
        };

        bool isNodeInSharedCommMap[] = {
                false, false, false, false,
                true, true, true, true,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false
        };

        bool isNodeInAuraCommMap[] = {
                true, true, true, true,
                false, false, false, false,
                true, true, true, true,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false
        };

        bool isNodeValid[] = {
               true, true, true, true,
               true, true, true, true,
               true, true, true, true,
               false, false, false, false,
               false, false, false, false,
               false, false, false, false,
               false, false, false, false
        };

        checkCommMaps("After on Proc0", stkMeshBulkData, numElements, ownerOfElement, isElementInAuraCommMap, isElementValid,
                numNodes, ownerOfNode, isNodeInSharedCommMap, isNodeInAuraCommMap, isNodeValid);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void runProc1(stk::mesh::BulkData &stkMeshBulkData)
{
    {
        int numElements = 6;
        bool ownerOfElement[] = { false, false, true, true, false, false };
        bool isElementInAuraCommMap[] = { false, true, true, true, true, false };
        bool isElementValid[] = { false, true, true, true, true, false };

        int numNodes = 28;
        bool ownerOfNode[] = {
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                true, true, true, true,
                true, true, true, true,
                false, false, false, false,
                false, false, false, false
        };

        bool isNodeInSharedCommMap[] = {
                false, false, false, false,
                false, false, false, false,
                true, true, true, true,
                false, false, false, false,
                true, true, true, true,
                false, false, false, false,
                false, false, false, false
        };

        bool isNodeInAuraCommMap[] = {
                false, false, false, false,
                true, true, true, true,
                false, false, false, false,
                true, true, true, true,
                false, false, false, false,
                true, true, true, true,
                false, false, false, false
        };

        bool isNodeValid[] = {
                false, false, false, false,
                true, true, true, true,
                true, true, true, true,
                true, true, true, true,
                true, true, true, true,
                true, true, true, true,
               false, false, false, false
        };

        checkCommMaps("Before on Proc 1", stkMeshBulkData, numElements, ownerOfElement, isElementInAuraCommMap, isElementValid,
                numNodes, ownerOfNode, isNodeInSharedCommMap, isNodeInAuraCommMap, isNodeValid);
    }

    stk::mesh::Entity entity = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, 3);
    std::vector<EntityProc> entitiesToMove;
    int proc2=2;
    entitiesToMove.push_back(std::make_pair(entity, proc2));

    add_nodes_to_move(stkMeshBulkData, entity, proc2, entitiesToMove);

    stkMeshBulkData.change_entity_owner(entitiesToMove);

    {
            int numElements = 6;
            bool ownerOfElement[] = { false, false, false, true, false, false };
            bool isElementInAuraCommMap[] = { false, true, true, true, true, false };
            bool isElementValid[] = { false, false, true, true, true, false };

            int numNodes = 28;
            bool ownerOfNode[] = {
                    false, false, false, false,
                    false, false, false, false,
                    false, false, false, false,
                    false, false, false, false,
                    true, true, true, true,
                    false, false, false, false,
                    false, false, false, false
            };

            bool isNodeInSharedCommMap[] = {
                    false, false, false, false,
                    false, false, false, false,
                    false, false, false, false,
                    true, true, true, true,
                    true, true, true, true,
                    false, false, false, false,
                    false, false, false, false
            };

            bool isNodeInAuraCommMap[] = {
                    false, false, false, false,
                    false, false, false, false,
                    true, true, true, true,
                    false, false, false, false,
                    false, false, false, false,
                    true, true, true, true,
                    false, false, false, false
            };

            bool isNodeValid[] = {
                    false, false, false, false,
                    false, false, false, false,
                    true, true, true, true,
                    true, true, true, true,
                    true, true, true, true,
                    true, true, true, true,
                   false, false, false, false
            };

            checkCommMaps("After on Proc 1", stkMeshBulkData, numElements, ownerOfElement, isElementInAuraCommMap, isElementValid,
                    numNodes, ownerOfNode, isNodeInSharedCommMap, isNodeInAuraCommMap, isNodeValid);
        }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void runProc2(stk::mesh::BulkData &stkMeshBulkData)
{
    {
        int numElements = 6;
        bool ownerOfElement[] = { false, false, false, false, true, true };
        bool isElementInAuraCommMap[] = { false, false, false, true, true, false };
        bool isElementValid[] = { false, false, false, true, true, true };

        int numNodes = 28;
        bool ownerOfNode[] = {
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                true, true, true, true,
                true, true, true, true
        };

        bool isNodeInSharedCommMap[] = {
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                true, true, true, true,
                false, false, false, false,
                false, false, false, false
        };

        bool isNodeInAuraCommMap[] = {
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                true, true, true, true,
                false, false, false, false,
                true, true, true, true,
                false, false, false, false
        };

        bool isNodeValid[] = {
                false, false, false, false,
                false, false, false, false,
                false, false, false, false,
                true, true, true, true,
                true, true, true, true,
                true, true, true, true,
                true, true, true, true
        };

        checkCommMaps("Before on Proc 2", stkMeshBulkData, numElements, ownerOfElement, isElementInAuraCommMap, isElementValid,
                numNodes, ownerOfNode, isNodeInSharedCommMap, isNodeInAuraCommMap, isNodeValid);
    }

    std::vector<EntityProc> entitiesToMove;

    stkMeshBulkData.change_entity_owner(entitiesToMove);

    {
        int numElements = 6;
        bool ownerOfElement[] = { false, true, true, false, true, true };
        bool isElementInAuraCommMap[] = { true, true, true, true, true, false };
        bool isElementValid[] = { true, true, true, true, true, true };

        int numNodes = 28;
        bool ownerOfNode[] = {
                false, false, false, false,
                true, true, true, true,
                true, true, true, true,
                true, true, true, true,
                false, false, false, false,
                true, true, true, true,
                true, true, true, true
        };

        bool isNodeInSharedCommMap[] = {
                false, false, false, false,
                true, true, true, true,
                false, false, false, false,
                true, true, true, true,
                true, true, true, true,
                false, false, false, false,
                false, false, false, false
        };

        bool isNodeInAuraCommMap[] = {
                true, true, true, true,               // left of element 1
                false, false, false, false,           // element 1 : element 2
                true, true, true, true,               // element 2 : element 3
                false, false, false, false,           // element 3 : element 4
                false, false, false, false,           // element 4 : element 5
                true, true, true, true,               // element 5 : element 6
                false, false, false, false            // right of element 6
        };

        bool isNodeValid[] = {
                true, true, true, true,
                true, true, true, true,
                true, true, true, true,
                true, true, true, true,
                true, true, true, true,
                true, true, true, true,
                true, true, true, true
        };

        checkCommMaps("After on Proc 2", stkMeshBulkData, numElements, ownerOfElement, isElementInAuraCommMap, isElementValid,
                numNodes, ownerOfNode, isNodeInSharedCommMap, isNodeInAuraCommMap, isNodeValid);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void moveElements2And3ToProc2(stk::mesh::BulkData &stkMeshBulkData)
{
    if ( stkMeshBulkData.parallel_rank() == 0 )
    {
        runProc0(stkMeshBulkData);
    }
    else if ( stkMeshBulkData.parallel_rank() == 1 )
    {
        runProc1(stkMeshBulkData);
    }
    else if ( stkMeshBulkData.parallel_rank() == 2 )
    {
        runProc2(stkMeshBulkData);
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void createSerialSubMesh(const stk::mesh::MetaData &meta, stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Selector subMeshSelector, stk::mesh::MetaData &newMeta, stk::mesh::BulkData &newBulkData)
{
    const stk::mesh::PartVector &allparts = meta.get_mesh_parts();
    for(size_t i = 0; i < allparts.size(); i++)
    {
        stk::mesh::Part *part = NULL;
        if(allparts[i]->topology() != stk::topology::INVALID_TOPOLOGY)
        {
            part = &newMeta.declare_part_with_topology(allparts[i]->name(), allparts[i]->topology());
        }
        else
        {
            part = &newMeta.declare_part(allparts[i]->name());
        }
        stk::io::put_io_part_attribute(*part);
    }

    const stk::mesh::FieldVector &fields = meta.get_fields();
    for(size_t i = 0; i < fields.size(); i++)
    {
        stk::mesh::FieldBase* newField = newMeta.declare_field_base(fields[i]->name(), fields[i]->entity_rank(), fields[i]->data_traits(),
                fields[i]->field_array_rank(), fields[i]->dimension_tags(), fields[i]->number_of_states());

        stk::mesh::Selector selectFieldParts = stk::mesh::selectField(*fields[i]);
        stk::mesh::PartVector oldParts;
        selectFieldParts.get_parts(oldParts);
        if(!oldParts.empty())
        {
            stk::mesh::PartVector newParts(oldParts.size());
            for(size_t k=0; k<oldParts.size(); k++)
            {
                newParts[k] = &newMeta.get_part(oldParts[k]->mesh_meta_data_ordinal());
            }
            stk::mesh::Selector selectNewParts = stk::mesh::selectUnion(newParts);
            stk::mesh::put_field(*newField, selectNewParts, newMeta.spatial_dimension());
        }
    }

    newMeta.commit();

    std::map<stk::mesh::Entity, stk::mesh::Entity> oldToNewEntityMap;

    newBulkData.modification_begin();
    for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEMENT_RANK; rank++)
    {
        const stk::mesh::BucketVector &buckets = stkMeshBulkData.get_buckets(rank, subMeshSelector);
        for(size_t i = 0; i < buckets.size(); i++)
        {
            stk::mesh::Bucket &bucket = *buckets[i];
            {
                for(size_t j = 0; j < bucket.size(); j++)
                {
                    const stk::mesh::PartVector &oldParts = bucket.supersets();
                    stk::mesh::PartVector newParts(oldParts.size(), 0);
                    for(size_t k = 0; k < oldParts.size(); k++)
                    {
                        newParts[k] = &newMeta.get_part(oldParts[k]->mesh_meta_data_ordinal());
                    }
                    stk::mesh::Entity newEntity = newBulkData.declare_entity(rank, stkMeshBulkData.identifier(bucket[j]), newParts);
                    oldToNewEntityMap[bucket[j]] = newEntity;

                    if ( rank != stk::topology::NODE_RANK )
                    {
                        unsigned numNodes = stkMeshBulkData.num_nodes(bucket[j]);
                        const stk::mesh::Entity *nodes = stkMeshBulkData.begin_nodes(bucket[j]);
                        for(unsigned nodeOrder = 0; nodeOrder < numNodes; nodeOrder++)
                        {
                            newBulkData.declare_relation(newEntity, oldToNewEntityMap[nodes[nodeOrder]], nodeOrder);
                        }
                    }
                }
            }
        }
    }
    newBulkData.modification_end();

    const stk::mesh::FieldVector &newFields = newMeta.get_fields();
    for(stk::mesh::EntityRank rank = stk::topology::NODE_RANK; rank <= stk::topology::ELEMENT_RANK; rank++)
    {
        const stk::mesh::BucketVector &buckets = stkMeshBulkData.get_buckets(rank, subMeshSelector);
        for(size_t i = 0; i < buckets.size(); i++)
        {
            stk::mesh::Bucket &bucket = *buckets[i];
            {
                for(size_t k = 0; k < fields.size(); k++)
                {
                    if(bucket.field_data_is_allocated(*fields[k]))
                    {
                        for(size_t j = 0; j < bucket.size(); j++)
                        {
                            stk::mesh::Entity oldEntity = bucket[j];
                            stk::mesh::Entity newEntity = oldToNewEntityMap[oldEntity];
                            void *oldData = stk::mesh::field_data(*fields[k], oldEntity);
                            void *newData = stk::mesh::field_data(*newFields[k], newEntity);
                            memcpy(newData, oldData, stk::mesh::field_bytes_per_entity(*fields[k], oldEntity));
                        }
                    }
                }
            }
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void testSubMesh(stk::mesh::BulkData &stkMeshBulkData, stk::mesh::Selector select, int elementToTestId, size_t goldNumberNodes, size_t goldNumberElements, FieldMgr &parallelFieldMgr)
{
    const stk::mesh::MetaData &stkMeshMetaData = stkMeshBulkData.mesh_meta_data();
    stk::mesh::MetaData newMetaData(3);
    stk::mesh::BulkData nBulkData(newMetaData, MPI_COMM_SELF);

    // Need to test io_part attribute item
    // Need to make sure field data is converted properly

    createSerialSubMesh(stkMeshMetaData, stkMeshBulkData, select, newMetaData, nBulkData);

    std::vector<unsigned> entityCounts;
    stk::mesh::count_entities(newMetaData.universal_part(), nBulkData, entityCounts);
    EXPECT_EQ(goldNumberElements, entityCounts[stk::topology::ELEMENT_RANK]);
    EXPECT_EQ(goldNumberNodes, entityCounts[stk::topology::NODE_RANK]);
    EXPECT_EQ(stkMeshMetaData.get_parts().size(), newMetaData.get_parts().size());
    EXPECT_EQ(stkMeshMetaData.get_fields().size(), newMetaData.get_fields().size());

    stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, elementToTestId);
    stk::mesh::Entity elementSerial = nBulkData.get_entity(stk::topology::ELEMENT_RANK, elementToTestId);

    EXPECT_EQ(stkMeshBulkData.identifier(element), nBulkData.identifier(elementSerial));

    std::set<stk::mesh::EntityId> parallelNodeIds;
    std::set<stk::mesh::EntityId> serialNodeIds;

    unsigned numNodesParallel = stkMeshBulkData.num_nodes(element);
    unsigned numNodesSerial = nBulkData.num_nodes(elementSerial);

    ASSERT_EQ(numNodesParallel, numNodesSerial);

    const stk::mesh::Entity *nodesParallel = stkMeshBulkData.begin_nodes(element);
    const stk::mesh::Entity *nodesSerial = nBulkData.begin_nodes(elementSerial);

    for(unsigned i = 0; i < numNodesParallel; i++)
    {
        parallelNodeIds.insert(stkMeshBulkData.identifier(nodesParallel[i]));
        serialNodeIds.insert(nBulkData.identifier(nodesSerial[i]));
    }

    EXPECT_TRUE(parallelNodeIds == serialNodeIds);

    FieldMgr serialFieldMgr(newMetaData);
    serialFieldMgr.storeFieldPointers();

    std::set<stk::mesh::EntityId>::iterator iter = serialNodeIds.begin();
    for (; iter != serialNodeIds.end(); iter++)
    {
        stk::mesh::Entity serialNode = nBulkData.get_entity(stk::topology::NODE_RANK, *iter);
        stk::mesh::Entity parallelNode = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, *iter);
        EXPECT_EQ(serialFieldMgr.getFieldValue(serialFieldMgr.getAuraCommMapNodeField(), serialNode),
                parallelFieldMgr.getFieldValue(parallelFieldMgr.getAuraCommMapNodeField(), parallelNode) );
        EXPECT_EQ(serialFieldMgr.getFieldValue(serialFieldMgr.getCommListNodeField(), serialNode),
                parallelFieldMgr.getFieldValue(parallelFieldMgr.getCommListNodeField(), parallelNode) );
        EXPECT_EQ(serialFieldMgr.getFieldValue(serialFieldMgr.getSharingCommMapNodeField(), serialNode),
                parallelFieldMgr.getFieldValue(parallelFieldMgr.getSharingCommMapNodeField(), parallelNode) );
    }

    EXPECT_EQ(serialFieldMgr.getFieldValue(serialFieldMgr.getAuraCommMapElementField(), elementSerial),
            parallelFieldMgr.getFieldValue(parallelFieldMgr.getAuraCommMapElementField(), element) );
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool isEntityValidOnCommList(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity)
{
    stk::mesh::EntityKey entityKey = stkMeshBulkData.entity_key(entity);
    stk::mesh::EntityCommListInfoVector::const_iterator iter = std::lower_bound(stkMeshBulkData.comm_list().begin(), stkMeshBulkData.comm_list().end(), entityKey);
    return iter != stkMeshBulkData.comm_list().end() && entityKey == iter->key;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool isEntityInGhostingCommMap(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity)
{
    stk::mesh::EntityKey entityKey = stkMeshBulkData.entity_key(entity);
    bool is_entity_in_aura_comm_map = !stkMeshBulkData.entity_comm_map(entityKey, stkMeshBulkData.aura_ghosting()).empty();
    return is_entity_in_aura_comm_map;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool isEntityInSharingCommMap(stk::mesh::BulkData& stkMeshBulkData, stk::mesh::Entity entity)
{
    stk::mesh::EntityKey entityKey = stkMeshBulkData.entity_key(entity);
    bool is_entity_in_shared_comm_map = !stkMeshBulkData.entity_comm_map(entityKey, stkMeshBulkData.shared_ghosting()).empty();
    return is_entity_in_shared_comm_map;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void checkCommMaps(std::string message, stk::mesh::BulkData &stkMeshBulkData, int numElements, bool ownerOfElement[], bool isElementInAuraCommMap[], bool isElementValid[],
            int numNodes, bool ownerOfNode[], bool isNodeInSharedCommMap[], bool isNodeInAuraCommMap[], bool isNodeValid[])
{
    for (int i=0;i<numElements;i++)
    {
        stk::mesh::Entity element = stkMeshBulkData.get_entity(stk::topology::ELEMENT_RANK, i+1);
        EXPECT_EQ(isElementValid[i], stkMeshBulkData.is_valid(element)) << message << " for element " << i+1;
        if ( isElementValid[i] && stkMeshBulkData.is_valid(element) )
        {
            bool amIOwner = ( stkMeshBulkData.parallel_owner_rank(element) == stkMeshBulkData.parallel_rank() );
            EXPECT_EQ(ownerOfElement[i], amIOwner) << message << " for element " << i+1;
            EXPECT_EQ(isElementInAuraCommMap[i], isEntityInGhostingCommMap(stkMeshBulkData, element)) << message << " for element " << i+1;
        }
    }

    for (int i=0;i<numNodes;i++)
    {
        stk::mesh::Entity node = stkMeshBulkData.get_entity(stk::topology::NODE_RANK, i+1);
        EXPECT_EQ(isNodeValid[i], stkMeshBulkData.is_valid(node)) << message << " for node " << i+1;
        if ( isNodeValid[i] && stkMeshBulkData.is_valid(node) )
        {
            bool amIOwner = ( stkMeshBulkData.parallel_owner_rank(node) == stkMeshBulkData.parallel_rank() );
            EXPECT_EQ(ownerOfNode[i], amIOwner) << message << " for node " << i+1;
            EXPECT_EQ(isNodeInSharedCommMap[i], isEntityInSharingCommMap(stkMeshBulkData, node)) << message << " for node " << i+1;
            EXPECT_EQ(isNodeInAuraCommMap[i], isEntityInGhostingCommMap(stkMeshBulkData, node)) << message << " for node " << i+1;
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void add_nodes_to_move(stk::mesh::BulkData& bulk,
                       stk::mesh::Entity elem,
                       int dest_proc,
                       std::vector<stk::mesh::EntityProc>& entities_to_move)
{
    const stk::mesh::Entity* nodes = bulk.begin_nodes(elem);
    for(unsigned i = 0; i < bulk.num_nodes(elem); ++i)
    {
        if(bulk.parallel_owner_rank(nodes[i]) == bulk.parallel_rank())
        {
            entities_to_move.push_back(stk::mesh::EntityProc(nodes[i], dest_proc));
        }
    }
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////

void putCommInfoDataOnFields(stk::mesh::BulkData &bulkData, FieldMgr &fieldMgr)
{
    const stk::mesh::BucketVector &nodeBuckets = bulkData.buckets(stk::topology::NODE_RANK);

    int ghostedToThisProcValue = 10;
    int ownedAndNotShareEntityValue = 0;
    int ghostedToAnotherProcValue = 8;
    int ownedAndSharedNodeValue = 2;
    int sharedOnlyNodeValue = 5;

    for(size_t i = 0; i < nodeBuckets.size(); i++)
    {
        const stk::mesh::Bucket &bucket = *nodeBuckets[i];
        //(bucket.owned())
        {
            for(size_t j = 0; j < bucket.size(); j++)
            {
                double *commListNode = stk::mesh::field_data(*fieldMgr.getCommListNodeField(), bucket[j]);
                double *sharedCommMapNode = stk::mesh::field_data(*fieldMgr.getSharingCommMapNodeField(), bucket[j]);
                double *auraCommMapNode = stk::mesh::field_data(*fieldMgr.getAuraCommMapNodeField(), bucket[j]);

                // 0: not ghosted
                // 1: ghosted else where
                // 2: ghosted here

                if(isEntityInGhostingCommMap(bulkData, bucket[j]))
                {
                    if ( bucket.in_aura() )
                    {
                        *auraCommMapNode = 2;
                    }
                    else
                    {
                        *auraCommMapNode = 1;
                    }
                }
                else
                {
                    *auraCommMapNode = 0;
                }

                // 0: owned only
                // 1: owned and shared
                // 2: shared only
                // 3: ghosted to me

                if(isEntityInSharingCommMap(bulkData, bucket[j]))
                {
                    if ( bucket.owned() )
                    {
                        *sharedCommMapNode = ownedAndSharedNodeValue;
                    }
                    else
                    {
                        *sharedCommMapNode = sharedOnlyNodeValue;
                    }
                }
                else
                {
                    if (bucket.in_aura() )
                    {
                        *sharedCommMapNode = ghostedToThisProcValue;
                    }
                    else
                    {
                        *sharedCommMapNode = ownedAndNotShareEntityValue;
                    }
                }

                if(isEntityValidOnCommList(bulkData, bucket[j]))
                {
                    *commListNode = 1;
                }
                else
                {
                    *commListNode = 0;
                }
            }
        }
    }

    const stk::mesh::BucketVector &elementBuckets = bulkData.buckets(stk::topology::ELEMENT_RANK);

    for(size_t i = 0; i < elementBuckets.size(); i++)
    {
        const stk::mesh::Bucket &bucket = *elementBuckets[i];
        //if(bucket.owned())
        {
            for(size_t j = 0; j < bucket.size(); j++)
            {
                double *auraCommMap = stk::mesh::field_data(*fieldMgr.getAuraCommMapElementField(), bucket[j]);

                // 0: Not ghosted (owned)
                // 1: ghosted to other proc
                // 2: ghosted here from another proc

                if(isEntityInGhostingCommMap(bulkData, bucket[j]))
                {
                    if ( bucket.in_aura() )
                    {
                        *auraCommMap = ghostedToThisProcValue;
                    }
                    else
                    {
                        *auraCommMap = ghostedToAnotherProcValue;
                    }
                }
                else
                {
                    *auraCommMap = ownedAndNotShareEntityValue;
                }
            }
        }
    }
}
//////////////////////////////////////////////////////////////////////////////////////////////////////////////

}
