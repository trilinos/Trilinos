#include "TextMesh.hpp"
#include <mpi.h>
#include <string>
#include <sstream>
#include <vector>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/FEMHelpers.hpp>
#include <stk_mesh/base/GetEntities.hpp>
#include "stk_util/environment/ReportHandler.hpp"

namespace stk
{
namespace unit_test_util
{

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

typedef std::map<int,std::set<stk::mesh::EntityId>> ProcToNodeIdsMap;

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

void add_node_sharing_for_intersecting_nodes(stk::mesh::BulkData &bulkData, const std::set<stk::mesh::EntityId>& intersectionNodeIds, const int proc)
{
    for(stk::mesh::EntityId nodeId : intersectionNodeIds)
    {
        stk::mesh::Entity node = bulkData.get_entity(stk::topology::NODE_RANK, nodeId);
        bulkData.add_node_sharing(node, proc);
    }
}

void setup_node_sharing(const MeshData& parsedData, stk::mesh::BulkData &bulkData)
{
    const ProcToNodeIdsMap procToNodeIdsMap = get_proc_to_node_ids_map(parsedData);
    const ProcToNodeIdsMap::const_iterator mapIterForMyNodeIds = procToNodeIdsMap.find(bulkData.parallel_rank());
    if (mapIterForMyNodeIds != procToNodeIdsMap.end())
    {
        for (const ProcToNodeIdsMap::value_type& procSetPair : procToNodeIdsMap)
        {
            if (bulkData.parallel_rank() == procSetPair.first)
                continue;
            std::set<stk::mesh::EntityId> intersectionNodeIds = find_intersecting_nodes(procSetPair.second, mapIterForMyNodeIds->second);
            add_node_sharing_for_intersecting_nodes(bulkData, intersectionNodeIds, procSetPair.first);
        }
    }
}

void setup_mesh(const MeshData& meshData, stk::mesh::BulkData &bulkData)
{
    bulkData.modification_begin();
    for(const ElementData& elementData : meshData.elementDataVec)
    {
        if(bulkData.parallel_rank() == elementData.proc || bulkData.parallel_size() == 1)
        {
            stk::mesh::PartVector topologyParts = {&bulkData.mesh_meta_data().get_topology_root_part(elementData.topology)};
            stk::mesh::declare_element(bulkData,
                                       topologyParts,
                                       elementData.identifier,
                                       elementData.nodeIds);
        }
    }
    setup_node_sharing(meshData, bulkData);
    bulkData.modification_end();
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




void fill_mesh_using_text_mesh(const std::string &meshDesc, stk::mesh::BulkData &bulkData)
{
    MeshData meshData = parse_input(meshDesc);
    setup_mesh(meshData, bulkData);
}

} // namespace unit_test_util
} // namespace stk
