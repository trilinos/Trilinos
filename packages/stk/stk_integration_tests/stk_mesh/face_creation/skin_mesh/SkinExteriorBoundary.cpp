/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <Ioss_IOFactory.h>             // for IOFactory
#include <Ioss_Region.h>                // for Region
#include <init/Ionit_Initializer.h>     // for Initializer
#include <stddef.h>                     // for size_t, nullptr
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD
#include "stk_io/DatabasePurpose.hpp"
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture

//#define VERBOSE_OUTPUT

namespace
{

/////////////////////////////////////////////////////////////////////////////////////////////////////////////

typedef std::map<std::string, unsigned int> TestCaseData;
typedef TestCaseData::value_type TestCaseDatum;

TestCaseData get_test_cases()
{

//The Magical Alphabet of Hexes, Shells & Sidesets
//
// A = hex in block A
// B = hex in block B
// e = shell in block E
// f = shell in block F
// L = sideset associated with the side on the left
// R = "        "           "   "  "     "   " right
// D = sideset containing 2 sides, one associated to left and one to right
// X = sideset associated with all sides on this surface
// J = two hexes in block A connected to the same 8 nodes
//
// .e = the language of our Patron Saint Exodus
//
// RR = pronounced like a pirate
// RRR = roll the R

    static TestCaseData test_cases = {
            {"AA.e",      10u},
            {"AB.e",      10u},
            {"ADA.e",     10u},
            {"ADB.e",     10u},
            {"ADDA.e",    10u},
            {"ADDB.e",    10u},
            {"ADeDA.e",   10u},
            {"ADeDB.e",   10u},
            {"ADe.e",      6u},
            {"ADeLA.e",   10u},
            {"ADeLB.e",   10u},
            {"ADeRA.e",   10u},
            {"ADeRB.e",   10u},
            {"A.e",        6u},
            {"AeA.e",     10u},
            {"AeB.e",     10u},
            {"Ae.e",       6u},
            {"AefA.e",    10u},
            {"AefB.e",    10u},
            {"Aef.e",      6u},
            {"ALA.e",     10u},
            {"ALB.e",     10u},
            {"AL.e",       6u},
            {"ALeDA.e",   10u},
            {"ALeDB.e",   10u},
            {"ALeDfRA.e", 10u},
            {"ALeDfRB.e", 10u},
            {"ALe.e",      6u},
            {"ALeLA.e",   10u},
            {"ALeLB.e",   10u},
            {"ALeRA.e",   10u},
            {"ALeRB.e",   10u},
            {"ALRA.e",    10u},
            {"ALRB.e",    10u},
            {"ARA.e",     10u},
            {"ARB.e",     10u},
            {"AReDA.e",   10u},
            {"AReDB.e",   10u},
            {"ARe.e",      6u},
            {"AReLA.e",   10u},
            {"AReLB.e",   10u},
            {"AReRA.e",   10u},
            {"AReRB.e",   10u},
            {"e.e",        2u},
            {"eL.e",       2u},
            {"ALeXfRA.e", 10u},
            {"ADReA.e",   10u},
            {"ADReB.e",   10u},
            {"ALReA.e",   10u},
            {"ALReB.e",   10u},
            {"ARReA.e",   10u},
            {"ARReB.e",   10u},
            {"Re.e",       2u},
            {"ReL.e",      2u},
            {"ALefRA.e",  10u},
            {"ARefLA.e",  10u},
            {"AeDfA.e",   10u},
            {"ALJ.e",     10u}
    };

    return test_cases;
}

class SkinnedMesh: public stk::unit_test_util::MeshTestFixture
{
protected:
    void input_from_file(const std::string &meshSpec, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        initialize_mesh();
        setup_empty_mesh(auraOption);

#if defined(VERBOSE_OUTPUT)
        if(get_bulk().parallel_rank() == 0) std::cerr << "Reading " << meshSpec << std::endl;
#endif

        stk::unit_test_util::read_from_serial_file_and_decompose(meshSpec, get_bulk(), "cyclic");
    }

    stk::mesh::Part & run_skin_mesh(stk::mesh::Selector blocksToSkin)
    {
        stk::mesh::Part &skin = get_meta().declare_part("skin", get_meta().side_rank());
        EXPECT_NO_FATAL_FAILURE(stk::mesh::create_exposed_boundary_sides(get_bulk(), blocksToSkin, skin));
        return skin;
    }

    void test_read_file(const std::string &meshSpec, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        EXPECT_NO_FATAL_FAILURE(input_from_file(meshSpec, auraOption));
    }

    void test_skin_file(const TestCaseDatum& testCase, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        EXPECT_NO_FATAL_FAILURE(input_from_file(testCase.first, auraOption));
        stk::mesh::Part &skin = run_skin_mesh(get_skin_block_selector());
        EXPECT_EQ(testCase.second, get_global_number_of_skinned_faces(skin));
        expect_aura_correct(auraOption, testCase.second, skin);
        EXPECT_TRUE(stk::mesh::check_exposed_boundary_sides(get_bulk(), get_skin_block_selector(), skin));
    }

    stk::mesh::Selector get_skin_block_selector()
    {
        return get_meta().universal_part();
    }

private:
    std::vector<unsigned> get_entity_count(const stk::mesh::Part &skin)
    {
        std::vector<unsigned> skin_counts;
        stk::mesh::count_entities(skin & get_meta().locally_owned_part(), get_bulk(), skin_counts);
        return skin_counts;
    }

    unsigned get_global_number_of_skinned_faces(const stk::mesh::Part &skin)
    {
        std::vector<unsigned> skin_counts = get_entity_count(skin);
        unsigned localSkinnedCount = skin_counts[get_meta().side_rank()], globalSkinnedCount = 0;
        stk::all_reduce_sum<unsigned>( get_comm(), &localSkinnedCount, &globalSkinnedCount , 1 );
        return globalSkinnedCount;
    }

    void expect_ghost_sides_connected_to_ghost_elements(stk::mesh::Part& skin)
    {
        stk::mesh::EntityVector ghostSides;
        stk::mesh::get_selected_entities((get_meta().aura_part() & skin), get_bulk().buckets(get_meta().side_rank()), ghostSides);
        for(stk::mesh::Entity ghostSide : ghostSides)
            expect_ghost_side_connected(ghostSide);
    }

    void expect_aura_correct(stk::mesh::BulkData::AutomaticAuraOption auraOption,
                             unsigned numGlobalFaces,
                             stk::mesh::Part& skin)
    {
        if(auraOption == stk::mesh::BulkData::AUTO_AURA && stk::mesh::count_selected_entities(get_meta().universal_part(),
                                                                                              get_bulk().buckets(stk::topology::ELEM_RANK)))
        {
            EXPECT_EQ(numGlobalFaces, stk::mesh::count_selected_entities(skin, get_bulk().buckets(get_meta().side_rank())));
            expect_ghost_sides_connected_to_ghost_elements(skin);
        }
    }

    void expect_ghost_side_connected(stk::mesh::Entity ghostSide)
    {
        unsigned numConnectedElems = get_bulk().num_elements(ghostSide);
        EXPECT_TRUE(numConnectedElems > 0);
        const stk::mesh::Entity* connectedElems = get_bulk().begin_elements(ghostSide);
        for(unsigned i = 0; i < numConnectedElems; i++)
            EXPECT_TRUE(get_bulk().bucket(connectedElems[i]).in_aura());
    }
};

class ReadMesh: public SkinnedMesh
{
protected:
    virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        TestCaseData test_cases = get_test_cases();
        for(const TestCaseDatum& testCase : test_cases)
        {
            test_read_file(testCase.first, auraOption);
        }
    }
};

class BasicSkinnedMesh: public SkinnedMesh
{
protected:
    virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        TestCaseData test_cases = get_test_cases();
        for(const TestCaseDatum& testCase : test_cases)
        {
            test_skin_file(testCase, auraOption);
        }
    }
};

class SkinnedMeshWithModifiedSkinPart: public SkinnedMesh
{
protected:
    virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        EXPECT_NO_FATAL_FAILURE(input_from_file("ARA.e", auraOption));
        stk::mesh::Part &skin = run_skin_mesh(get_skin_block_selector());
        run_modification(skin);
        EXPECT_FALSE(stk::mesh::check_exposed_boundary_sides(get_bulk(), get_skin_block_selector(), skin));
    }

    virtual void run_modification(stk::mesh::Part &skin) = 0;

    stk::mesh::EntityVector get_faces(stk::mesh::Selector selector)
    {
        stk::mesh::EntityVector faces;
        stk::mesh::get_selected_entities(selector, get_bulk().buckets(stk::topology::FACE_RANK), faces);
        return faces;
    }
};

class SkinnedMeshWithExtraFace: public SkinnedMeshWithModifiedSkinPart
{
protected:
    virtual void run_modification(stk::mesh::Part &skin)
    {
        get_bulk().modification_begin();
        add_extra_face_to_skin(skin);
        get_bulk().modification_end();
    }

private:
    void add_extra_face_to_skin(stk::mesh::Part &skin)
    {
        if(get_bulk().parallel_rank() == 1) {
            stk::mesh::EntityVector notSkinFaces = get_faces(!skin);
            ASSERT_EQ(1u, notSkinFaces.size());
            get_bulk().change_entity_parts(notSkinFaces[0], {&skin});
        }
    }
};

class SkinnedMeshWithMissingFace: public SkinnedMeshWithModifiedSkinPart
{
protected:
    virtual void run_modification(stk::mesh::Part &skin)
    {
        get_bulk().modification_begin();
        remove_face_from_skin(skin);
        get_bulk().modification_end();
    }

private:
    void remove_face_from_skin(stk::mesh::Part &skin)
    {
        stk::mesh::EntityVector skinFaces = get_faces(skin);
        ASSERT_EQ(5u, skinFaces.size());
        get_bulk().change_entity_parts(skinFaces[0], {}, {&skin});
    }
};

TEST_F(ReadMesh, read_all_files_aura)
{
    run_test(stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(ReadMesh, read_all_files_no_aura)
{
    run_test(stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST_F(BasicSkinnedMesh, skin_all_files_aura)
{
    run_test(stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(BasicSkinnedMesh, skin_all_files_no_aura)
{
    run_test(stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST_F(SkinnedMeshWithExtraFace, skin_no_aura)
{
    run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}

TEST_F(SkinnedMeshWithMissingFace, skin_no_aura)
{
    run_test_on_num_procs(2, stk::mesh::BulkData::NO_AUTO_AURA);
}

} //namespace
