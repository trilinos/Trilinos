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
#include <stk_mesh/base/GetEntities.hpp>
#include <stk_mesh/base/Comm.hpp>
#include <stk_mesh/base/SkinBoundary.hpp>
#include <stk_unit_test_utils/ioUtils.hpp>
#include <stk_unit_test_utils/MeshFixture.hpp>  // for MeshTestFixture
#include <stk_mesh/baseImpl/EquivalentEntityBlocks.hpp>

//#define VERBOSE_OUTPUT

namespace
{

/////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef std::map<std::string, unsigned> TestCaseData;
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
// Z = degenerate hex in block Z
// Y = degenerate hex in block Y
// T = tet in block T
// g = degenerate shell (usually attached to face of tet T)
// .e = the language of our Patron Saint Exodus
//
// RR = pronounced like a pirate
// RRR = roll the R

    static TestCaseData test_cases = {
        //   Filename      Number Element Blocks (GOLD)
            {"AA.e",       1u},
            {"AB.e",       2u},
            {"ADA.e",      1u},
            {"ADB.e",      2u},
            {"ADDA.e",     1u},
            {"ADDB.e",     2u},
            {"ADeDA.e",    2u},
            {"ADeDB.e",    3u},
            {"ADe.e",      2u},
            {"ADeLA.e",    2u},
            {"ADeLB.e",    3u},
            {"ADeRA.e",    2u},
            {"ADeRB.e",    3u},
            {"A.e",        1u},
            {"AeA.e",      2u},
            {"AeB.e",      3u},
            {"Ae.e",       2u},
            {"AefA.e",     3u},
            {"AefB.e",     4u},
            {"Aef.e",      3u},
            {"ALA.e",      1u},
            {"ALB.e",      2u},
            {"AL.e",       1u},
            {"ALeDA.e",    2u},
            {"ALeDB.e",    3u},
            {"ALeDfRA.e",  3u},
            {"ALeDfRB.e",  4u},
            {"ALe.e",      2u},
            {"ALeLA.e",    2u},
            {"ALeLB.e",    3u},
            {"ALeRA.e",    2u},
            {"ALeRB.e",    3u},
            {"ALRA.e",     1u},
            {"ALRB.e",     2u},
            {"ARA.e",      1u},
            {"ARB.e",      2u},
            {"AReDA.e",    2u},
            {"AReDB.e",    3u},
            {"ARe.e",      2u},
            {"AReLA.e",    2u},
            {"AReLB.e",    3u},
            {"AReRA.e",    2u},
            {"AReRB.e",    3u},
            {"e.e",        1u},
            {"eL.e",       1u},
            {"ALeXfRA.e",  3u},
            {"ADReA.e",    2u},
            {"ADReB.e",    3u},
            {"ALReA.e",    2u},
            {"ALReB.e",    3u},
            {"ARReA.e",    2u},
            {"ARReB.e",    3u},
            {"Re.e",       1u},
            {"ReL.e",      1u},
            {"ALefRA.e",   3u},
            {"ARefLA.e",   3u},
            {"AeDfA.e",    3u},
            {"ALJ.e",      1u},
            {"TDg.e",      2u},
            {"ZDZ.e",      1u},
            {"Tg.e",       2u},
            {"ZY.e",       2u},
            {"ef.e",       2u},
            {"AB_doubleKissing.e",      2u},
            {"ADDA_doubleKissing.e",    1u},
            {"ALRA_doubleKissing.e",    1u},
            {"ALLA_doubleKissing.e",    1u},
            {"ARRA_doubleKissing.e",    1u}
    };

    return test_cases;
}

class LoadMesh: public stk::unit_test_util::MeshTestFixture
{
public:
    virtual ~LoadMesh() {}

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

    stk::mesh::EntityVector get_all_elements()
    {
        stk::mesh::EntityVector elements;
        stk::mesh::get_selected_entities(get_meta().universal_part(), get_bulk().buckets(stk::topology::ELEMENT_RANK), elements);
        return elements;
    }

};

class TwoElementMesh: public LoadMesh
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
        check_difference(allElements[0], allElements[1]);
    }

    virtual void check_difference(stk::mesh::Entity element1, stk::mesh::Entity element2) = 0;

private:
    std::string m_inputFile;
};

class TwoElementsInSameBlock: public TwoElementMesh
{
public:
    TwoElementsInSameBlock()
       : TwoElementMesh("AA.e")
    { }

protected:
    virtual void check_difference(stk::mesh::Entity element1, stk::mesh::Entity element2)
    {
        EXPECT_TRUE(stk::mesh::impl::are_entity_element_blocks_equivalent(get_bulk(), element1, element2));
    }
};

class TwoElementsInDifferentBlocks: public TwoElementMesh
{
public:
    TwoElementsInDifferentBlocks()
       : TwoElementMesh("AB.e")
    { }

protected:
    virtual void check_difference(stk::mesh::Entity element1, stk::mesh::Entity element2)
    {
        EXPECT_FALSE(stk::mesh::impl::are_entity_element_blocks_equivalent(get_bulk(), element1, element2));
    }
};


TEST_F(TwoElementsInSameBlock, test_difference_with_aura)
{
    run_test_on_num_procs_or_less(2, stk::mesh::BulkData::AUTO_AURA);
}

TEST_F(TwoElementsInDifferentBlocks, test_difference_with_aura)
{
    run_test_on_num_procs_or_less(2, stk::mesh::BulkData::AUTO_AURA);
}

} //namespace
