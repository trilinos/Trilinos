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

#include <stk_util/stk_config.h>

#ifdef STK_HAVE_STKIO

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
    TopologyMapper(const std::string &exodusName_, const int exodusNumNodes_,
                   const std::string& iossTopologyName_, stk::topology stkTopology_)
      : exodusName(exodusName_), exodusNumNodes(exodusNumNodes_),
        iossTopologyName(iossTopologyName_), stkTopology(stkTopology_)
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
    std::string iossTopologyName;
    stk::topology stkTopology;

    exodusName="sphere";
    exodusNumNodes=1;
    iossTopologyName="sphere";
    stkTopology=stk::topology::PARTICLE;
    topologyMappings.push_back(TopologyMapper(exodusName, exodusNumNodes, iossTopologyName, stkTopology));

    exodusName="BEam";
    exodusNumNodes=3;
    iossTopologyName="bar3";
    stkTopology=stk::topology::BEAM_3;
    topologyMappings.push_back(TopologyMapper(exodusName, exodusNumNodes, iossTopologyName, stkTopology));

    exodusName="Tri";
    exodusNumNodes=3;
    iossTopologyName="trishell3";
    stkTopology=stk::topology::SHELL_TRIANGLE_3;
    topologyMappings.push_back(TopologyMapper(exodusName, exodusNumNodes, iossTopologyName, stkTopology));

    exodusName="hex";
    exodusNumNodes=20;
    iossTopologyName="hex20";
    stkTopology=stk::topology::HEXAHEDRON_20;
    topologyMappings.push_back(TopologyMapper(exodusName, exodusNumNodes, iossTopologyName, stkTopology));

    exodusName="quad";
    exodusNumNodes=6;
    iossTopologyName="quad6";
    stkTopology=stk::topology::QUAD_6;
    topologyMappings.push_back(TopologyMapper(exodusName, exodusNumNodes, iossTopologyName, stkTopology));

    exodusName="wedge";
    exodusNumNodes=12;
    iossTopologyName="wedge12";
    stkTopology=stk::topology::WEDGE_12;
    topologyMappings.push_back(TopologyMapper(exodusName, exodusNumNodes, iossTopologyName, stkTopology));
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

        stk::topology mappedStkTopologyFromIossTopology = stk::io::map_ioss_topology_to_stk(iossTopology, spatialDim);
        EXPECT_EQ(goldValues.stkTopology, mappedStkTopologyFromIossTopology);
    }
}

//End documentation test here

}

#endif

