#include <gtest/gtest.h>
#include <stk_topology/topology.hpp>
#include <Ioss_Utils.h>
#include <Ioss_ElementTopology.h>
// #include <Ioss_Initializer.h>
#include <init/Ionit_Initializer.h>
#include <stk_io/IossBridge.hpp>
#include <stk_mesh/base/MetaData.hpp>

#include <vector>

namespace {

struct TopologyMapper
{
    std::string exodusName;
    int exodusNumNodes;
    std::string iossTopologyName;
    stk::topology stkTopology;
    stk::mesh::CellTopology shardsTopology;
    TopologyMapper(const std::string &exodusName_, const int exodusNumNodes_, const std::string& iossTopologyName_,
            stk::topology stkTopology_, stk::mesh::CellTopology shardsTopology_) :
                exodusName(exodusName_), exodusNumNodes(exodusNumNodes_), iossTopologyName(iossTopologyName_),
                stkTopology(stkTopology_), shardsTopology(shardsTopology_)
    {}
};

void createIossElementRegistryForKnownElementTopologies()
{
    Ioss::Init::Initializer io;
}

//Begin documentation test here

void setUpMappingsToTest(std::vector<TopologyMapper>& topologyMappings)
{
    std::string exodusName;
    int exodusNumNodes=-1;
    std::string iossTopology;
    stk::topology stkTopology;
    stk::mesh::CellTopology shardsTopology;

    exodusName="sphere";
    exodusNumNodes=1;
    iossTopology="sphere";
    stkTopology=stk::topology::PARTICLE;
    shardsTopology=stk::mesh::CellTopology(shards::getCellTopologyData< shards::Particle >());
    topologyMappings.push_back(TopologyMapper(exodusName, exodusNumNodes, iossTopology, stkTopology, shardsTopology));

    exodusName="BEam";
    exodusNumNodes=3;
    iossTopology="bar3";
    stkTopology=stk::topology::BEAM_3;
    shardsTopology=stk::mesh::CellTopology(shards::getCellTopologyData< shards::Beam<3> >());
    topologyMappings.push_back(TopologyMapper(exodusName, exodusNumNodes, iossTopology, stkTopology, shardsTopology));

    exodusName="Tri";
    exodusNumNodes=3;
    iossTopology="trishell3";
    stkTopology=stk::topology::SHELL_TRIANGLE_3;
    shardsTopology=stk::mesh::CellTopology(shards::getCellTopologyData< shards::ShellTriangle<3> >());
    topologyMappings.push_back(TopologyMapper(exodusName, exodusNumNodes, iossTopology, stkTopology, shardsTopology));

    exodusName="hex";
    exodusNumNodes=20;
    iossTopology="hex20";
    stkTopology=stk::topology::HEXAHEDRON_20;
    shardsTopology=stk::mesh::CellTopology(shards::getCellTopologyData< shards::Hexahedron<20> >());
    topologyMappings.push_back(TopologyMapper(exodusName, exodusNumNodes, iossTopology, stkTopology, shardsTopology));
}

TEST(Understanding, sierra_topologies)
{
    int spatialDim = 3;
    std::vector<TopologyMapper> topologyMappings;
    setUpMappingsToTest(topologyMappings);

    size_t numMappings = topologyMappings.size();

    createIossElementRegistryForKnownElementTopologies();

    for (size_t i=0;i<numMappings;i++)
    {
        TopologyMapper goldValues = topologyMappings[i];

        std::string fixedExodusName = Ioss::Utils::fixup_type(topologyMappings[i].exodusName, topologyMappings[i].exodusNumNodes, spatialDim);
        Ioss::ElementTopology *iossTopology = Ioss::ElementTopology::factory(fixedExodusName, true);
        ASSERT_TRUE(iossTopology != NULL);
        EXPECT_EQ(goldValues.iossTopologyName, iossTopology->name());

        stk::topology mappedStkTopologyFromIossTopology = stk::io::map_ioss_topology_to_stk(iossTopology);
        EXPECT_EQ(goldValues.stkTopology, mappedStkTopologyFromIossTopology);

        stk::mesh::CellTopology mappedShardsTopologyFromStkTopology = stk::mesh::get_cell_topology(mappedStkTopologyFromIossTopology);
        EXPECT_EQ(goldValues.shardsTopology, mappedShardsTopologyFromStkTopology);

        stk::topology mappedStkTopologyFromShards = stk::mesh::get_topology(mappedShardsTopologyFromStkTopology, spatialDim);
        EXPECT_EQ(goldValues.stkTopology, mappedStkTopologyFromShards);
    }
}

//End documentation test here

}
