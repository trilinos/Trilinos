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
#include <string>                       // for string
#include "mpi.h"                        // for MPI_COMM_WORLD
#include <stk_io/StkMeshIoBroker.hpp>   // for StkMeshIoBroker
#include "stk_io/DatabasePurpose.hpp"
#include <stk_io/IossBridge.hpp>
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
struct MeshInfoData {
    unsigned int numberOfSkinnedFaces;
    unsigned int numberOfIOElementBlocks;
};

typedef std::map<std::string, MeshInfoData> TestCaseData;
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
            {"AA.e",      {10u, 1u}},
            {"AB.e",      {10u, 2u}},
            {"ADA.e",     {10u, 1u}},
            {"ADB.e",     {10u, 2u}},
            {"ADDA.e",    {10u, 1u}},
            {"ADDB.e",    {10u, 2u}},
            {"ADeDA.e",   {10u, 2u}},
            {"ADeDB.e",   {10u, 3u}},
            {"ADe.e",     { 6u, 2u}},
            {"ADeLA.e",   {10u, 2u}},
            {"ADeLB.e",   {10u, 3u}},
            {"ADeRA.e",   {10u, 2u}},
            {"ADeRB.e",   {10u, 3u}},
            {"A.e",       { 6u, 1u}},
            {"AeA.e",     {10u, 2u}},
            {"AeB.e",     {10u, 3u}},
            {"Ae.e",      { 6u, 2u}},
            {"AefA.e",    {10u, 3u}},
            {"AefB.e",    {10u, 4u}},
            {"Aef.e",     { 6u, 3u}},
            {"ALA.e",     {10u, 1u}},
            {"ALB.e",     {10u, 2u}},
            {"AL.e",      { 6u, 1u}},
            {"ALeDA.e",   {10u, 2u}},
            {"ALeDB.e",   {10u, 3u}},
            {"ALeDfRA.e", {10u, 3u}},
            {"ALeDfRB.e", {10u, 4u}},
            {"ALe.e",     { 6u, 2u}},
            {"ALeLA.e",   {10u, 2u}},
            {"ALeLB.e",   {10u, 3u}},
            {"ALeRA.e",   {10u, 2u}},
            {"ALeRB.e",   {10u, 3u}},
            {"ALRA.e",    {10u, 1u}},
            {"ALRB.e",    {10u, 2u}},
            {"ARA.e",     {10u, 1u}},
            {"ARB.e",     {10u, 2u}},
            {"AReDA.e",   {10u, 2u}},
            {"AReDB.e",   {10u, 3u}},
            {"ARe.e",     { 6u, 2u}},
            {"AReLA.e",   {10u, 2u}},
            {"AReLB.e",   {10u, 3u}},
            {"AReRA.e",   {10u, 2u}},
            {"AReRB.e",   {10u, 3u}},
            {"e.e",       { 2u, 1u}},
            {"eL.e",      { 2u, 1u}},
            {"ALeXfRA.e", {10u, 3u}},
            {"ADReA.e",   {10u, 2u}},
            {"ADReB.e",   {10u, 3u}},
            {"ALReA.e",   {10u, 2u}},
            {"ALReB.e",   {10u, 3u}},
            {"ARReA.e",   {10u, 2u}},
            {"ARReB.e",   {10u, 3u}},
            {"Re.e",      { 2u, 1u}},
            {"ReL.e",     { 2u, 1u}},
            {"ALefRA.e",  {10u, 3u}},
            {"ARefLA.e",  {10u, 3u}},
            {"AeDfA.e",   {10u, 3u}},
            {"ALJ.e",     {10u, 1u}}
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

    bool is_element_block(const stk::mesh::Part &part)
    {
        return ((part.primary_entity_rank() == stk::topology::ELEMENT_RANK) &&
                part.topology().is_valid()) && stk::io::is_part_io_part(part);
    }

    void get_element_blocks(stk::mesh::PartVector &element_blocks)
    {
        const stk::mesh::PartVector partVector = get_meta().get_mesh_parts();
        for(stk::mesh::Part *part : partVector) {
            if(is_element_block(*part))
                element_blocks.push_back(part);
        }
    }

    void get_element_blocks_for_entity(stk::mesh::Entity entity, stk::mesh::PartVector &element_blocks)
    {
        const stk::mesh::PartVector &partVector = get_bulk().bucket(entity).supersets();
        for(stk::mesh::Part *part : partVector) {
            if(is_element_block(*part))
                element_blocks.push_back(part);
        }
    }

    stk::mesh::EntityVector get_all_elements()
    {
        stk::mesh::EntityVector elements;
        stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::ELEMENT_RANK), elements);
        return elements;
    }
};

class ElementBlockPartTest: public SkinnedMesh
{
protected:
    virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        TestCaseData test_cases = get_test_cases();
        for(const TestCaseDatum& testCase : test_cases)
        {
            test_element_blocks(testCase, auraOption);
        }
    }

    void test_element_blocks(const TestCaseDatum& testCase, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        std::vector<stk::mesh::Part*> element_blocks;
        EXPECT_NO_FATAL_FAILURE(input_from_file(testCase.first, auraOption));
        get_element_blocks(element_blocks);
        EXPECT_EQ(testCase.second.numberOfIOElementBlocks, element_blocks.size());
    }
};

class TwoElementMesh: public SkinnedMesh
{
public:
    TwoElementMesh(std::string inputFile)
    : m_inputFile(inputFile)
    {}

protected:
    virtual void run_test(stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        EXPECT_NO_FATAL_FAILURE(input_from_file(m_inputFile, auraOption));
        stk::mesh::EntityVector allElements = get_all_elements();
        EXPECT_EQ(2u, allElements.size());
        run_difference_test(allElements);
    }

    virtual void check_difference(stk::mesh::PartVector &blocksForElement1, stk::mesh::PartVector &blocksForElement2) = 0;

    stk::mesh::PartVector get_difference(stk::mesh::PartVector &blocksForElement1, stk::mesh::PartVector &blocksForElement2)
    {
        stk::mesh::PartVector diff(blocksForElement1.size() + blocksForElement2.size());
        std::sort(blocksForElement1.begin(), blocksForElement1.end());
        std::sort(blocksForElement2.begin(), blocksForElement2.end());
        stk::mesh::PartVector::iterator it=std::set_difference (blocksForElement1.begin(), blocksForElement1.end(), blocksForElement2.begin(), blocksForElement2.end(), diff.begin());
        diff.resize(it-diff.begin());
        return diff;
    }

private:
    void run_difference_test(stk::mesh::EntityVector &allElements)
    {
        stk::mesh::PartVector blocksForElement1, blocksForElement2;
        get_element_blocks_for_entity(allElements[0], blocksForElement1);
        get_element_blocks_for_entity(allElements[1], blocksForElement2);
        check_difference(blocksForElement1, blocksForElement2);
    }

    std::string m_inputFile;
};

class TwoElementsInSameBlock: public TwoElementMesh
{
public:
    TwoElementsInSameBlock()
       : TwoElementMesh("AA.e")
    { }

protected:
    virtual void check_difference(stk::mesh::PartVector &blocksForElement1, stk::mesh::PartVector &blocksForElement2)
    {
        stk::mesh::PartVector diff12 = get_difference(blocksForElement1, blocksForElement2);
        stk::mesh::PartVector diff21 = get_difference(blocksForElement2, blocksForElement1);
        EXPECT_EQ(0u, diff12.size());
        EXPECT_EQ(0u, diff21.size());
    }
};

class TwoElementsInDifferentBlocks: public TwoElementMesh
{
public:
    TwoElementsInDifferentBlocks()
       : TwoElementMesh("AB.e")
    { }

protected:
    void check_difference(stk::mesh::PartVector &blocksForElement1, stk::mesh::PartVector &blocksForElement2)
    {
        stk::mesh::PartVector diff12 = get_difference(blocksForElement1, blocksForElement2);
        stk::mesh::PartVector diff21 = get_difference(blocksForElement2, blocksForElement1);
        EXPECT_EQ(1u, diff12.size());
        EXPECT_EQ(1u, diff21.size());
    }
};

TEST_F(ElementBlockPartTest, is_element_block_no_aura)
{
    run_test(stk::mesh::BulkData::NO_AUTO_AURA);
}


TEST_F(ElementBlockPartTest, is_element_block_aura)
{
    run_test(stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(TwoElementsInSameBlock, test_difference_with_aura)
{
    run_test_on_num_procs_or_less(2, stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(TwoElementsInDifferentBlocks, test_difference_with_aura)
{
    run_test_on_num_procs_or_less(2, stk::mesh::BulkData::AUTO_AURA);
}

} //namespace
