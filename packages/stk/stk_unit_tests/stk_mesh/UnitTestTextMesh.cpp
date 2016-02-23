#include "gtest/gtest.h"
#include <mpi.h>
#include <string>
#include <sstream>
#include <vector>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>

namespace {

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems) {
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.push_back(item);
    }
    return elems;
}


std::vector<std::string> split(const std::string &s, char delim) {
    std::vector<std::string> elems;
    split(s, delim, elems);
    return elems;
}

std::string remove_spaces(const std::string& data)
{
    std::string str(data);
    std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
    str.erase(end_pos, str.end());
    return str;
}

std::string make_upper_case(const std::string& data)
{
    std::string str(data);
    std::transform(str.begin(), str.end(), str.begin(), ::toupper);
    return str;
}

stk::topology get_topology_by_name(const std::string& name)
{
    std::map<std::string,stk::topology> nameToTopology =
    {
       {  "NODE"          , stk::topology::NODE         },
       {  "LINE_2"        , stk::topology::LINE_2       },
       {  "LINE_3"        , stk::topology::LINE_3       },
       {  "TRI_3"         , stk::topology::TRI_3        },
       {  "TRI_4"         , stk::topology::TRI_4        },
       {  "TRI_6"         , stk::topology::TRI_6        },
       {  "QUAD_4"        , stk::topology::QUAD_4       },
       {  "QUAD_8"        , stk::topology::QUAD_8       },
       {  "QUAD_9"        , stk::topology::QUAD_9       },
       {  "PARTICLE"      , stk::topology::PARTICLE     },
       {  "LINE_2_1D"     , stk::topology::LINE_2_1D    },
       {  "LINE_3_1D"     , stk::topology::LINE_3_1D    },
       {  "BEAM_2"        , stk::topology::BEAM_2       },
       {  "BEAM_3"        , stk::topology::BEAM_3       },
       {  "SHELL_LINE_2"  , stk::topology::SHELL_LINE_2 },
       {  "SHELL_LINE_3"  , stk::topology::SHELL_LINE_3 },
       {  "TRI_3_2D"      , stk::topology::TRI_3_2D     },
       {  "TRI_4_2D"      , stk::topology::TRI_4_2D     },
       {  "TRI_6_2D"      , stk::topology::TRI_6_2D     },
       {  "QUAD_4_2D"     , stk::topology::QUAD_4_2D    },
       {  "QUAD_8_2D"     , stk::topology::QUAD_8_2D    },
       {  "QUAD_9_2D"     , stk::topology::QUAD_9_2D    },
       {  "SHELL_TRI_3"   , stk::topology::SHELL_TRI_3  },
       {  "SHELL_TRI_4"   , stk::topology::SHELL_TRI_4  },
       {  "SHELL_TRI_6"   , stk::topology::SHELL_TRI_6  },
       {  "SHELL_QUAD_4"  , stk::topology::SHELL_QUAD_4 },
       {  "SHELL_QUAD_8"  , stk::topology::SHELL_QUAD_8 },
       {  "SHELL_QUAD_9"  , stk::topology::SHELL_QUAD_9 },
       {  "TET_4"         , stk::topology::TET_4        },
       {  "TET_8"         , stk::topology::TET_8        },
       {  "TET_10"        , stk::topology::TET_10       },
       {  "TET_11"        , stk::topology::TET_11       },
       {  "PYRAMID_5"     , stk::topology::PYRAMID_5    },
       {  "PYRAMID_13"    , stk::topology::PYRAMID_13   },
       {  "PYRAMID_14"    , stk::topology::PYRAMID_14   },
       {  "WEDGE_6"       , stk::topology::WEDGE_6      },
       {  "WEDGE_15"      , stk::topology::WEDGE_15     },
       {  "WEDGE_18"      , stk::topology::WEDGE_18     },
       {  "HEX_8"         , stk::topology::HEX_8        },
       {  "HEX_20"        , stk::topology::HEX_20       },
       {  "HEX_27"        , stk::topology::HEX_27       },
    };

    stk::topology topology = stk::topology::INVALID_TOPOLOGY;
    std::map<std::string,stk::topology>::iterator iter = nameToTopology.find(name);
    if (iter != nameToTopology.end())
    {
        topology = iter->second;
    }
    return topology;

}

struct ElementData
{
    int proc;
    stk::topology topology;
    stk::mesh::EntityId identifier;
    stk::mesh::EntityIdVector nodeIds;
};

struct MeshData
{
    int spatialDim;
    std::vector<ElementData> elementDataVec;
};

class TextMesh
{
public:
    typedef std::map<int,std::set<stk::mesh::EntityId>> ProcToNodeIdsMap;

    TextMesh(const std::string& meshDescription, stk::mesh::BulkData::AutomaticAuraOption auraOption)
    {
        MeshData meshData = parse_input(meshDescription);
        determine_number_of_procs(meshData);
        meta = new stk::mesh::MetaData(meshData.spatialDim);
        bulk = new stk::mesh::BulkData(*meta, MPI_COMM_WORLD, auraOption);
        setup_mesh(meshData);
    }
    ~TextMesh() { delete bulk ; delete meta;}
    bool num_procs_ok() { return get_bulk().parallel_size() == numberOfProcs; }
    stk::mesh::BulkData& get_bulk() { return *bulk; }
    stk::mesh::MetaData& get_meta() { return *meta; }

private:

    void setup_mesh(const MeshData& meshData)
    {
        get_bulk().modification_begin();
        for(const ElementData& elementData : meshData.elementDataVec)
        {
            if(get_bulk().parallel_rank() == elementData.proc)
            {
                stk::mesh::PartVector topologyParts = {&get_meta().get_topology_root_part(elementData.topology)};
                stk::mesh::declare_element(get_bulk(),
                                           topologyParts,
                                           elementData.identifier,
                                           elementData.nodeIds);
            }
        }
        setup_node_sharing(meshData);
        get_bulk().modification_end();
    }

    void determine_number_of_procs(MeshData& meshData)
    {
        int maxProc = 0;
        for (const ElementData& elementData : meshData.elementDataVec)
            maxProc = std::max(maxProc,elementData.proc);
        numberOfProcs = maxProc+1;
    }

    MeshData parse_input(const std::string& meshDescription)
    {
        MeshData data;
        int spatialDim = -1;
        int spatialDimLine = -1;
        const std::vector<std::string> lines = split(meshDescription,'\n');
        for (size_t lineI=0 ; lineI<lines.size() ; ++lineI)
        {
            size_t userLineNumber = lineI+1;
            std::string line = make_upper_case(remove_spaces(lines[lineI]));
            const std::vector<std::string> tokens = split(line,',');
            ThrowRequireMsg(tokens.size()>=4, "Error!  Each line must contain the following fields (with at least one node):  Processor, Element Topology, GlobalId, NodeIds.  Error on line " << userLineNumber << ".");
            ElementData elementData;
            elementData.proc = std::stoi(tokens[0]);
            elementData.identifier = static_cast<stk::mesh::EntityId>(std::stoi(tokens[1]));
            elementData.topology = get_topology_by_name(tokens[2]);
            ThrowRequireMsg(elementData.topology != stk::topology::INVALID_TOPOLOGY, "Error!  Topology = >>" << tokens[2] << "<< is invalid from line " << userLineNumber << ".");
            if (-1 == spatialDim)
            {
                spatialDim = elementData.topology.dimension();
                spatialDimLine = userLineNumber;
            }
            else
                ThrowRequireMsg(elementData.topology.defined_on_spatial_dimension(spatialDim), "Error!  Topology = " << elementData.topology << " is not defined on spatial dimension = " << spatialDim << " that was set on line " << spatialDimLine << ".  Error on line " << userLineNumber << ".");
            unsigned numNodes = elementData.topology.num_nodes();
            ThrowRequireMsg(tokens.size() == numNodes+3u, "Error!  The input line contains " << tokens.size()-3 << " nodes, but the topology " << elementData.topology.name() << " needs " << numNodes << " nodes on line " << userLineNumber << ".");
            elementData.nodeIds.resize(numNodes);
            for (unsigned i=0 ; i < numNodes ; ++i)
            {
                elementData.nodeIds[i] = static_cast<stk::mesh::EntityId>(std::stoi(tokens[3+i]));
            }
            data.elementDataVec.push_back(elementData);
        }
        ThrowRequireMsg(spatialDim>1, "Error!  Spatial dimension not defined to be 2 or 3!");
        data.spatialDim = spatialDim;
        return data;
    }


    ProcToNodeIdsMap get_proc_to_node_ids_map(const MeshData& meshData)
    {
        ProcToNodeIdsMap procToNodeIdsMap;
        for (const ElementData& elementData : meshData.elementDataVec)
        {
            procToNodeIdsMap[elementData.proc].insert(
                    elementData.nodeIds.begin(),
                    elementData.nodeIds.end()
                    );
        }
        return procToNodeIdsMap;
    }

    void setup_node_sharing(const MeshData& parsedData)
    {
        const ProcToNodeIdsMap procToNodeIdsMap = get_proc_to_node_ids_map(parsedData);
        const ProcToNodeIdsMap::const_iterator mapIterForMyNodeIds = procToNodeIdsMap.find(get_bulk().parallel_rank());
        if (mapIterForMyNodeIds != procToNodeIdsMap.end())
        {
            for (const ProcToNodeIdsMap::value_type& procSetPair : procToNodeIdsMap)
            {
                if (get_bulk().parallel_rank() == procSetPair.first)
                    continue;
                std::set<stk::mesh::EntityId> intersectionNodeIds = find_intersecting_nodes(procSetPair.second, mapIterForMyNodeIds->second);
                add_node_sharing_for_intersecting_nodes(intersectionNodeIds, procSetPair.first);
            }
        }
    }

    std::set<stk::mesh::EntityId> find_intersecting_nodes(const std::set<stk::mesh::EntityId>& otherProcNodeIds,
                                                          const std::set<stk::mesh::EntityId>& myNodeIds)
    {
        std::set<stk::mesh::EntityId> intersectionNodeIds;
        std::set_intersection(otherProcNodeIds.begin(),
                              otherProcNodeIds.end(),
                              myNodeIds.begin(),
                              myNodeIds.end(),
                              std::insert_iterator<std::set<stk::mesh::EntityId> >(intersectionNodeIds,
                                                                                   intersectionNodeIds.begin()));
        return intersectionNodeIds;
    }

    void add_node_sharing_for_intersecting_nodes(const std::set<stk::mesh::EntityId>& intersectionNodeIds, const int proc)
    {
        for(stk::mesh::EntityId nodeId : intersectionNodeIds)
        {
            stk::mesh::Entity node = get_bulk().get_entity(stk::topology::NODE_RANK, nodeId);
            get_bulk().add_node_sharing(node, proc);
        }
    }

    stk::mesh::BulkData* bulk;
    stk::mesh::MetaData* meta;
    int numberOfProcs = -1;
};

void verify_nodes(const stk::mesh::BulkData& bulk, stk::mesh::Entity element, const stk::mesh::EntityIdVector& goldNodeIds)
{
    stk::mesh::EntityVector nodes(bulk.begin_nodes(element), bulk.end_nodes(element));
    stk::mesh::EntityIdVector nodeIds(nodes.size());
    for(size_t i = 0; i < nodes.size(); ++i)
    {
        nodeIds[i] = bulk.identifier(nodes[i]);
    }
    EXPECT_EQ(goldNodeIds, nodeIds);
}

void verify_shared_nodes(const stk::mesh::BulkData& bulk, const stk::mesh::EntityIdVector& nodeIds , int sharingProc)
{
    for (stk::mesh::EntityId nodeId : nodeIds)
        EXPECT_TRUE(bulk.in_shared(stk::mesh::EntityKey(stk::topology::NODE_RANK,nodeId),sharingProc));
}

void verify_single_element(const stk::mesh::BulkData& bulk, stk::mesh::EntityId elemId, stk::topology topology, const stk::mesh::EntityIdVector& nodeIds)
{
    stk::mesh::Entity element = bulk.get_entity(stk::topology::ELEM_RANK, elemId);
    EXPECT_TRUE(bulk.is_valid(element));
    EXPECT_EQ(topology, bulk.bucket(element).topology());
    verify_nodes(bulk, element, nodeIds);
}

TEST(TextMesh, singlHex)
{
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8";
    TextMesh textMesh(meshDesc, stk::mesh::BulkData::NO_AUTO_AURA);
    if (!textMesh.num_procs_ok())
        return;
    verify_single_element(textMesh.get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
}

TEST(TextMesh, singleQuad)
{
    std::string meshDesc = "0,1,QUAD_4_2D,1,2,3,4";
    TextMesh textMesh(meshDesc, stk::mesh::BulkData::NO_AUTO_AURA);
    if (!textMesh.num_procs_ok())
        return;
    verify_single_element(textMesh.get_bulk(), 1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1,2,3,4});
}

TEST(TextMesh, mixedSpatialDim)
{
    std::string meshDesc =
"0,1,HEX_8,1,2,3,4,5,6,7,8\n\
0,2,QUAD_4_2D,5,6,7,8";
    EXPECT_THROW(new TextMesh(meshDesc, stk::mesh::BulkData::NO_AUTO_AURA), std::logic_error);
}

TEST(TextMesh, singlHexWithSpaces)
{
    std::string meshDesc = "0, 1, HEX_8, 1, 2, 3, 4, 5, 6, 7, 8";
    TextMesh textMesh(meshDesc, stk::mesh::BulkData::NO_AUTO_AURA);
    if (!textMesh.num_procs_ok())
        return;
    verify_single_element(textMesh.get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
}

TEST(TextMesh, singlHexWithLowerCase)
{
    std::string meshDesc = "0,1,Hex_8,1,2,3,4,5,6,7,8";
    TextMesh textMesh(meshDesc, stk::mesh::BulkData::NO_AUTO_AURA);
    if (!textMesh.num_procs_ok())
        return;
    verify_single_element(textMesh.get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
}

TEST(TextMesh, tooFewNodes)
{
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7";
    EXPECT_THROW(new TextMesh(meshDesc, stk::mesh::BulkData::NO_AUTO_AURA), std::logic_error);
}

TEST(TextMesh, tooManyNodes)
{
    std::string meshDesc = "0,1,HEX_8,1,2,3,4,5,6,7,8,9";
    EXPECT_THROW(new TextMesh(meshDesc, stk::mesh::BulkData::NO_AUTO_AURA), std::logic_error);
}

TEST(TextMesh, tooLittleData)
{
    std::string meshDesc = "0,1,";
    EXPECT_THROW(new TextMesh(meshDesc, stk::mesh::BulkData::NO_AUTO_AURA), std::logic_error);
}

TEST(TextMesh, invalidTopology)
{
    std::string meshDesc = "0,1,invalid,1";
    EXPECT_THROW(new TextMesh(meshDesc, stk::mesh::BulkData::NO_AUTO_AURA), std::logic_error);
}


TEST(TextMesh, twoHexesSerial)
{
    std::string meshDesc =
"0,1,HEX_8,1,2,3,4,5,6,7,8\n\
0,2,HEX_8,5,6,7,8,9,10,11,12";
    TextMesh textMesh(meshDesc, stk::mesh::BulkData::NO_AUTO_AURA);
    if (!textMesh.num_procs_ok())
        return;
    EXPECT_EQ(1,textMesh.get_bulk().parallel_size());
    verify_single_element(textMesh.get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
    verify_single_element(textMesh.get_bulk(), 2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{5,6,7,8,9,10,11,12});
}

TEST(TextMesh, twoHexesParallel)
{
    std::string meshDesc =
"0,1,HEX_8,1,2,3,4,5,6,7,8\n\
1,2,HEX_8,5,6,7,8,9,10,11,12";
    TextMesh textMesh(meshDesc, stk::mesh::BulkData::NO_AUTO_AURA);
    if (!textMesh.num_procs_ok())
        return;
    EXPECT_EQ(2,textMesh.get_bulk().parallel_size());
    if (textMesh.get_bulk().parallel_rank() == 0)
    {
        verify_single_element(textMesh.get_bulk(), 1u, stk::topology::HEX_8, stk::mesh::EntityIdVector{1,2,3,4,5,6,7,8});
        verify_shared_nodes(textMesh.get_bulk(), {5,6,7,8}, 1);
    }
    else
    {
        verify_single_element(textMesh.get_bulk(), 2u, stk::topology::HEX_8, stk::mesh::EntityIdVector{5,6,7,8,9,10,11,12});
        verify_shared_nodes(textMesh.get_bulk(), {5,6,7,8}, 0);
    }
}

TEST(TextMesh, DISABLED_twoQuadOneShellParallel)
{
    std::string meshDesc =
"0,1,QUAD_4_2D,1,2,3,4\n\
1,2,QUAD_4_2D,3,4,5,6\n\
2,3,SHELL_LINE_2,3,4";
    TextMesh textMesh(meshDesc, stk::mesh::BulkData::NO_AUTO_AURA);
    if (!textMesh.num_procs_ok())
        return;
    EXPECT_EQ(3,textMesh.get_bulk().parallel_size());
    if (textMesh.get_bulk().parallel_rank() == 0)
    {
        verify_single_element(textMesh.get_bulk(), 1u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{1,2,3,4});
        verify_shared_nodes(textMesh.get_bulk(), {3,4}, 1);
        verify_shared_nodes(textMesh.get_bulk(), {3,4}, 2);
    }
    else if (textMesh.get_bulk().parallel_rank() == 1)
    {
        verify_single_element(textMesh.get_bulk(), 2u, stk::topology::QUAD_4_2D, stk::mesh::EntityIdVector{3,4,5,6});
        verify_shared_nodes(textMesh.get_bulk(), {3,4}, 0);
        verify_shared_nodes(textMesh.get_bulk(), {3,4}, 2);
    }
    else
    {
        verify_single_element(textMesh.get_bulk(), 3u, stk::topology::SHELL_LINE_2, stk::mesh::EntityIdVector{3,4});
        verify_shared_nodes(textMesh.get_bulk(), {3,4}, 0);
        verify_shared_nodes(textMesh.get_bulk(), {3,4}, 1);

    }
}

} // namespace
