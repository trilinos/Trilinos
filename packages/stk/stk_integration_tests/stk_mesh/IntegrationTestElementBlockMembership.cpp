/*--------------------------------------------------------------------*/
/*    Copyright 2009 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <gtest/gtest.h>                // for AssertHelper, EXPECT_EQ, etc
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

class LoadMesh: public stk::unit_test_util::MeshTestFixture
{
public:
    virtual ~LoadMesh() {}

protected:
    void input_from_file(const std::string &meshSpec, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        reset_mesh();
        setup_empty_mesh(auraOption);

#if defined(VERBOSE_OUTPUT)
        if(get_bulk().parallel_rank() == 0) std::cout << "Reading " << meshSpec << std::endl;
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
