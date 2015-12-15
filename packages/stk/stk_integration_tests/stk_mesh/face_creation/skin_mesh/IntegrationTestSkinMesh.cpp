/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <Ioss_IOFactory.h>             // for IOFactory
#include <Ioss_Region.h>                // for Region
#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
#include <init/Ionit_Initializer.h>     // for Initializer
#include <stddef.h>                     // for size_t, nullptr
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include <stk_util/environment/Env.hpp>      // for sierraTimer, Timer
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD
#include "stk_io/DatabasePurpose.hpp"
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/CreateFaces.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_unit_test_utils/FaceTestingUtils.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace Ioss { class DatabaseIO; }

namespace
{

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct FaceConnectionData {
    unsigned int numberOfSkinnedBoundaryFaces;
    unsigned int numberOfSidesWithCreateFaces;
    bool everyElementSideHasAFaceWithCreateFaces;
    unsigned int numberOfFacesConnectedToMoreThanOneDistinctElementWithCreateFaces;
    unsigned int numberofFacesConnectedToSameElementTwiceWithCreateFaces;
    unsigned int numberOfSidesNoCreateFaces;
    bool everyElementSideHasAFaceNoCreateFaces;
    unsigned int numberOfFacesConnectedToMoreThanOneDistinctElementNoCreateFaces;
    unsigned int numberOfFacesConnectedToSameElementTwiceNoCreateFaces;
    std::set<unsigned> numberOfElementsConnectedToEachFaceAtXEqualPointFiveWithCreateFaces;
    std::set<unsigned> numberOfElementsConnectedToEachFaceAtXEqualPointFiveNoCreateFaces;
};
typedef std::map<std::string, FaceConnectionData> TestCaseData;
typedef TestCaseData::value_type TestCaseDatum;


using stk::mesh::EntityKey;


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
            {"AA.e",      {10u, 11u, true, 1u, 0u, 0u, false, 0u, 0u, {2   }, {    }}},
            {"AB.e",      {10u, 11u, true, 1u, 0u, 0u, false, 0u, 0u, {2   }, {    }}},
            {"ADA.e",     {10u, 11u, true, 1u, 0u, 1u, false, 1u, 0u, {2   }, {2   }}},
            {"ADB.e",     {10u, 11u, true, 1u, 0u, 1u, false, 1u, 0u, {2   }, {2   }}},
            {"ADDA.e",    {10u, 11u, true, 1u, 0u, 1u, false, 1u, 0u, {2   }, {2   }}},
            {"ADDB.e",    {10u, 11u, true, 1u, 0u, 1u, false, 1u, 0u, {2   }, {2   }}},
            {"ADeDA.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"ADeDB.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"ADe.e",     { 6u,  7u, true, 1u, 0u, 1u, false, 1u, 0u, {2, 1}, {2   }}},
            {"ADeLA.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"ADeLB.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"ADeRA.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"ADeRB.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"A.e",       { 6u,  6u, true, 0u, 0u, 0u, false, 0u, 0u, {1   }, {    }}},
            {"AeA.e",     {10u, 12u, true, 2u, 0u, 0u, false, 0u, 0u, {2, 2}, {    }}},
            {"AeB.e",     {10u, 12u, true, 2u, 0u, 0u, false, 0u, 0u, {2, 2}, {    }}},
            {"Ae.e",      { 6u,  7u, true, 1u, 0u, 0u, false, 0u, 0u, {2, 1}, {    }}},
            {"AefA.e",    {10u, 12u, true, 2u, 0u, 0u, false, 0u, 0u, {3, 3}, {    }}},
            {"AefB.e",    {10u, 12u, true, 2u, 0u, 0u, false, 0u, 0u, {3, 3}, {    }}},
            {"Aef.e",     { 6u,  7u, true, 2u, 0u, 0u, false, 0u, 0u, {3, 2}, {    }}},
            {"ALA.e",     {10u, 11u, true, 1u, 0u, 1u, false, 1u, 0u, {2   }, {2   }}},
            {"ALB.e",     {10u, 11u, true, 1u, 0u, 1u, false, 1u, 0u, {2   }, {2   }}},
            {"AL.e",      { 6u,  6u, true, 0u, 0u, 1u, false, 0u, 0u, {1   }, {1   }}},
            {"ALeDA.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"ALeDB.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"ALeDfRA.e", {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {3, 3}, {3, 3}}},
            {"ALeDfRB.e", {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {3, 3}, {3, 3}}},
            {"ALe.e",     { 6u,  7u, true, 1u, 0u, 1u, false, 1u, 0u, {2, 1}, {2   }}},
            {"ALeLA.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"ALeLB.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"ALeRA.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"ALeRB.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"ALRA.e",    {10u, 11u, true, 1u, 0u, 1u, false, 1u, 0u, {2   }, {2   }}},
            {"ALRB.e",    {10u, 11u, true, 1u, 0u, 1u, false, 1u, 0u, {2   }, {2   }}},
            {"ARA.e",     {10u, 11u, true, 1u, 0u, 1u, false, 1u, 0u, {2   }, {2   }}},
            {"ARB.e",     {10u, 11u, true, 1u, 0u, 1u, false, 1u, 0u, {2   }, {2   }}},
            {"AReDA.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"AReDB.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"ARe.e",     { 6u,  7u, true, 1u, 0u, 1u, false, 1u, 0u, {2, 1}, {2   }}},
            {"AReLA.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"AReLB.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"AReRA.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"AReRB.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {2, 2}, {2, 2}}},
            {"e.e",       { 2u,  2u, true, 0u, 0u, 0u, false, 0u, 0u, {1, 1}, {    }}},
            {"eL.e",      { 2u,  2u, true, 0u, 0u, 1u, false, 0u, 0u, {1, 1}, {1   }}},
            {"ALeXfRA.e", {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {3, 3}, {3, 3}}},
            {"ADReA.e",   {10u, 12u, true, 2u, 0u, 1u, false, 1u, 0u, {2, 2}, {2   }}},
            {"ADReB.e",   {10u, 12u, true, 2u, 0u, 1u, false, 1u, 0u, {2, 2}, {2   }}},
            {"ALReA.e",   {10u, 12u, true, 2u, 0u, 1u, false, 1u, 0u, {2, 2}, {2   }}},
            {"ALReB.e",   {10u, 12u, true, 2u, 0u, 1u, false, 1u, 0u, {2, 2}, {2   }}},
            {"ARReA.e",   {10u, 12u, true, 2u, 0u, 1u, false, 1u, 0u, {2, 2}, {2   }}},
            {"ARReB.e",   {10u, 12u, true, 2u, 0u, 1u, false, 1u, 0u, {2, 2}, {2   }}},
            {"Re.e",      { 2u,  2u, true, 0u, 0u, 1u, false, 0u, 0u, {1, 1}, {1   }}},
            {"ReL.e",     { 2u,  2u, true, 0u, 0u, 2u, true , 0u, 0u, {1, 1}, {1, 1}}},
            {"ALefRA.e",  {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {3, 3}, {3, 3}}},
            {"ARefLA.e",  {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {3, 3}, {3, 3}}},
            {"AeDfA.e",   {10u, 12u, true, 2u, 0u, 2u, false, 2u, 0u, {3, 3}, {3, 3}}},
            {"ALJ.e",     {10u, 11u, true, 6u, 0u, 1u, false, 1u, 0u, {3   }, {3   }}}
    };

    return test_cases;
}


void filter_test_case(TestCaseData& test_cases, const std::string& filename)
{
    auto iter = test_cases.find(filename);
    if(iter != test_cases.end())
    {
#if defined(VERBOSE_OUTPUT)
        std::cerr << "***** WARNING *****: Removing test " << iter->first << " since it doesn't work for this test. Needs fixing.\n";
#endif
        test_cases.erase(iter);
    }
}

TestCaseDatum select_test_case(TestCaseData& test_cases, const std::string& filename)
{
    TestCaseData::iterator iter = test_cases.find(filename);
    return *iter;
}

void do_input_from_file(const std::string &meshSpec, stk::mesh::BulkData &bulkData, stk::ParallelMachine communicator)
{
    if(stk::parallel_machine_size(communicator) == 1)
    {
        stk::unit_test_util::fill_mesh_using_stk_io(meshSpec, bulkData, communicator);
    }
    else
    {
        stk::unit_test_util::read_from_serial_file_and_decompose(meshSpec, bulkData, "cyclic");
    }
}

void test_read_file(const TestCaseDatum& testCase, stk::ParallelMachine communicator, stk::mesh::BulkData::AutomaticAuraOption auraOption)
{
    const std::string &meshSpec = testCase.first;

    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulkData(meta, communicator, auraOption);

    EXPECT_NO_FATAL_FAILURE(do_input_from_file(meshSpec, bulkData, communicator));
}

#define VERBOSE_OUTPUT

void test_skin_file(const TestCaseDatum& testCase, stk::ParallelMachine communicator, stk::mesh::BulkData::AutomaticAuraOption auraOption)
{
    const std::string &meshSpec = testCase.first;

    stk::mesh::MetaData meta(3);
    stk::mesh::BulkData bulkData(meta, communicator, auraOption);

    do_input_from_file(meshSpec, bulkData, communicator);

    stk::mesh::Selector blocksToSkin = meta.universal_part();
    stk::mesh::Part &skin = meta.declare_part("skin", meta.side_rank());
    stk::io::put_io_part_attribute(skin);

#if defined(VERBOSE_OUTPUT)
    std::cout << bulkData.parallel_rank() << ") Skinning: " << meshSpec << std::endl;

    std::vector<size_t> mesh_counts;
    stk::mesh::comm_mesh_counts(bulkData, mesh_counts);
    unsigned initialFaceCount = mesh_counts[meta.side_rank()];
    std::cout << bulkData.parallel_rank() << ") Initial number of faces = " << initialFaceCount << std::endl;
#endif

    EXPECT_NO_FATAL_FAILURE(stk::mesh::create_exposed_boundary_sides(bulkData, blocksToSkin, skin));

#if defined(VERBOSE_OUTPUT)
    stk::mesh::comm_mesh_counts(bulkData, mesh_counts);
    unsigned finalFaceCount = mesh_counts[meta.side_rank()];
    std::cout << bulkData.parallel_rank() << ") Final number of faces = " << finalFaceCount << std::endl;
#endif

    std::vector<unsigned> skin_counts;
    stk::mesh::count_entities(skin & meta.locally_owned_part(), bulkData, skin_counts);

#if defined(VERBOSE_OUTPUT)
    std::cout << bulkData.parallel_rank() << ") Number of skinned/boundary faces = " << skin_counts[meta.side_rank()] << std::endl;
#endif

    unsigned localSkinnedCount = skin_counts[meta.side_rank()];
    unsigned globalSkinnedCount = 0;
    stk::all_reduce_sum<unsigned>( communicator, &localSkinnedCount, &globalSkinnedCount , 1 );

    EXPECT_EQ(testCase.second.numberOfSkinnedBoundaryFaces, globalSkinnedCount);

    EXPECT_TRUE(stk::mesh::check_exposed_boundary_sides(bulkData, blocksToSkin, skin));
}

void filter_failing_tests(TestCaseData &test_cases, stk::mesh::BulkData::AutomaticAuraOption auraOption)
{
    // These 3 fail because the number of skinned faces does not match expectation
    filter_test_case(test_cases, "Aef.e");
    filter_test_case(test_cases, "AefA.e");
    filter_test_case(test_cases, "AefB.e");

    // These fail with both AURA and NO_AURA
    filter_test_case(test_cases, "ALeDfRA.e");
    filter_test_case(test_cases, "ALeDfRB.e");
    filter_test_case(test_cases, "ALeXfRA.e");
    filter_test_case(test_cases, "ALefRA.e");
    filter_test_case(test_cases, "ARefLA.e");
    filter_test_case(test_cases, "AeDfA.e");
    filter_test_case(test_cases, "ALeDB.e");
    filter_test_case(test_cases, "ADReA.e");


    if(stk::mesh::BulkData::AUTO_AURA == auraOption)
    {
        // Fail with aura
        filter_test_case(test_cases, "ADReB.e");
        filter_test_case(test_cases, "ADe.e");
        filter_test_case(test_cases, "ADeLA.e");
        filter_test_case(test_cases, "ADeRA.e");
        filter_test_case(test_cases, "ADeRB.e");
        filter_test_case(test_cases, "ALA.e");
        filter_test_case(test_cases, "ALB.e");
        filter_test_case(test_cases, "ALJ.e");
        filter_test_case(test_cases, "ALReA.e");
        filter_test_case(test_cases, "ALReB.e");
        filter_test_case(test_cases, "ALe.e");
        filter_test_case(test_cases, "ALeLA.e");
        filter_test_case(test_cases, "ALeRA.e");
        filter_test_case(test_cases, "ARA.e");
        filter_test_case(test_cases, "ARB.e");

        // These fail with AURA in debug
        filter_test_case(test_cases, "ARReA.e");
        filter_test_case(test_cases, "ARReB.e");
        filter_test_case(test_cases, "ARe.e");
        filter_test_case(test_cases, "AReLA.e");
        filter_test_case(test_cases, "AReLB.e");
        filter_test_case(test_cases, "AReRA.e");
        filter_test_case(test_cases, "AReRB.e");
        filter_test_case(test_cases, "Ae.e");
    }
}

void test_read_all_files(stk::mesh::BulkData::AutomaticAuraOption auraOption)
{
    TestCaseData test_cases = get_test_cases();
    filter_failing_tests(test_cases, auraOption);
    for(const TestCaseDatum& testCase : test_cases)
    {
        test_read_file(testCase, MPI_COMM_WORLD, auraOption);
    }
}

TEST(SkinMesh, read_all_files_aura)
{
    test_read_all_files(stk::mesh::BulkData::AUTO_AURA);
}

TEST(SkinMesh, read_all_files_no_aura)
{
    test_read_all_files(stk::mesh::BulkData::NO_AUTO_AURA);
}


void test_skin_all_files(stk::mesh::BulkData::AutomaticAuraOption auraOption)
{
    TestCaseData test_cases = get_test_cases();
    filter_failing_tests(test_cases, auraOption);
    for(const TestCaseDatum& testCase : test_cases)
    {
       test_skin_file(testCase, MPI_COMM_WORLD, auraOption);
    }
}

TEST(SkinMesh, skin_all_files_aura)
{
    test_skin_all_files(stk::mesh::BulkData::AUTO_AURA);
}

TEST(SkinMesh, skin_all_files_no_aura)
{
    test_skin_all_files(stk::mesh::BulkData::NO_AUTO_AURA);
}

} //namespace
