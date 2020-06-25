#include "balanceUtils.hpp"
#include "mpi.h"
#include "search_tolerance/FaceSearchTolerance.hpp"
#include "stk_mesh/base/Field.hpp"  // for field_data
#include "stk_mesh/base/FieldBase.hpp"  // for field_data
#include "stk_util/diag/StringUtil.hpp"
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_topology/topology.hpp>

namespace stk
{
namespace balance
{

//////////////////////////////////////////////////////////////////////////

BalanceSettings::BalanceSettings()
  : m_initialDecompMethod("RIB")
{}

size_t BalanceSettings::getNumNodesRequiredForConnection(stk::topology element1Topology, stk::topology element2Topology) const
{
    return 1;
}

double BalanceSettings::getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const
{
    return 1;
}

int BalanceSettings::getGraphVertexWeight(stk::topology type) const
{
    return 1;
}

double BalanceSettings::getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index) const
{
    return 1;
}

BalanceSettings::GraphOption BalanceSettings::getGraphOption() const
{
    return BalanceSettings::LOAD_BALANCE;
}

bool BalanceSettings::includeSearchResultsInGraph() const
{
    return false;
}

void BalanceSettings::setIncludeSearchResultsInGraph(bool doContactSearch) 
{
}

double BalanceSettings::getToleranceForFaceSearch(const stk::mesh::BulkData & mesh,
                                                  const stk::mesh::FieldBase & coordField,
                                                  const stk::mesh::Entity * faceNodes,
                                                  const unsigned numFaceNodes) const
{
    return 0.0;
}

void BalanceSettings::setToleranceFunctionForFaceSearch(std::shared_ptr<stk::balance::FaceSearchTolerance> faceSearchTolerance)
{
}

bool BalanceSettings::isConstantFaceSearchTolerance() const
{
    return true;
}

double BalanceSettings::getToleranceForParticleSearch() const
{
    return 0.0;
}

double BalanceSettings::getAbsoluteToleranceForParticleSearch(stk::mesh::Entity particle) const
{
    return getParticleRadius(particle) * getToleranceForParticleSearch();
}

double BalanceSettings::getGraphEdgeWeightForSearch() const
{
    return 1.0;
}

bool BalanceSettings::getEdgesForParticlesUsingSearch() const
{
    return false;
}

double BalanceSettings::getVertexWeightMultiplierForVertexInSearch() const
{
    return 15;
}

bool BalanceSettings::isIncrementalRebalance() const
{
    return false;
}

bool BalanceSettings::isMultiCriteriaRebalance() const
{
    return false;
}

bool BalanceSettings::areVertexWeightsProvidedInAVector() const
{
    return false;
}

bool BalanceSettings::areVertexWeightsProvidedViaFields() const
{
    return false;
}

std::vector<double> BalanceSettings::getVertexWeightsViaVector() const
{
    return std::vector<double>();
}

double BalanceSettings::getImbalanceTolerance() const
{
    return 1.01;
}

void BalanceSettings::setDecompMethod(const std::string& method)
{
}

std::string BalanceSettings::getDecompMethod() const
{
    return std::string("parmetis");
}

std::string BalanceSettings::getInitialDecompMethod() const
{
    return m_initialDecompMethod;
}

std::string BalanceSettings::getCoordinateFieldName() const
{
    return std::string("coordinates");
}

bool BalanceSettings::shouldPrintMetrics() const
{
    return false;
}

int BalanceSettings::getNumCriteria() const
{
    return 1;
}

void BalanceSettings::modifyDecomposition(DecompositionChangeList & decomp) const
{}

double BalanceSettings::getParticleRadius(stk::mesh::Entity particle) const
{
    return 0.5;
}

bool BalanceSettings::setVertexWeightsBasedOnNumberAdjacencies() const
{
    return false;
}

// For graph based methods (parmetis) only
bool BalanceSettings::allowModificationOfVertexWeightsForSmallMeshes() const
{
    return true;
}

// For graph based methods (parmetis) only
bool BalanceSettings::shouldFixMechanisms() const
{
    return false;
}

bool BalanceSettings::shouldFixSpiders() const
{
    return false;
}

std::string BalanceSettings::getSpiderBeamConnectivityCountFieldName() const
{
    return "beam_connectivity_count";
}

std::string BalanceSettings::getSpiderVolumeConnectivityCountFieldName() const
{
    return "volume_connectivity_count";
}

const stk::mesh::Field<int> * BalanceSettings::getSpiderBeamConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const
{
    return nullptr;
}

const stk::mesh::Field<int> * BalanceSettings::getSpiderVolumeConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const
{
    return nullptr;
}

bool BalanceSettings::useLocalIds() const
{
    return getGraphOption() == stk::balance::BalanceSettings::COLOR_MESH ||
           getGraphOption() == stk::balance::BalanceSettings::COLOR_MESH_BY_TOPOLOGY ||
           getGraphOption() == stk::balance::BalanceSettings::COLOR_MESH_AND_OUTPUT_COLOR_FIELDS;
}


bool BalanceSettings::useNodeBalancer() const
{
  return false;
}

double BalanceSettings::getNodeBalancerTargetLoadBalance() const
{
  return 1.0;
}

unsigned BalanceSettings::getNodeBalancerMaxIterations() const
{
  return 5;
}

void BalanceSettings::set_input_filename(const std::string& filename)
{
  m_inputFilename = filename;
}

std::string BalanceSettings::get_input_filename() const
{
  return m_inputFilename;
}

void BalanceSettings::set_output_filename(const std::string& filename)
{
  m_outputFilename = filename;
}

std::string BalanceSettings::get_output_filename() const
{
  return m_outputFilename;
}

//////////////////////////////////////

size_t GraphCreationSettings::getNumNodesRequiredForConnection(stk::topology element1Topology, stk::topology element2Topology) const
{
    const int noConnection = 1000;
    const int s = noConnection;
    const static int connectionTable[7][7] = {
        {1, 1, 1, 1, 1, 1, s}, // 0 dim
        {1, 1, 1, 1, 1, 1, s}, // 1 dim
        {1, 1, 2, 3, 2, 3, s}, // 2 dim linear
        {1, 1, 3, 3, 3, 3, s}, // 3 dim linear
        {1, 1, 2, 3, 3, 4, s}, // 2 dim higher-order
        {1, 1, 3, 3, 4, 4, s}, // 3 dim higher-order
        {s, s, s, s, s, s, s}  // super element
    };

    int element1Index = getConnectionTableIndex(element1Topology);
    int element2Index = getConnectionTableIndex(element2Topology);

    return connectionTable[element1Index][element2Index];
}

double GraphCreationSettings::getGraphEdgeWeightForSearch() const
{
    return edgeWeightForSearch;
}

double GraphCreationSettings::getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const
{
    const double defaultWeight = 1.0;
    const double lin3dlin3d = defaultWeight;
    const double noConnection = 0;
    const double s = noConnection;
    const double largeWeight = 5;
    const double L = largeWeight;
    const double twoDimWeight = 5;
    const double q = twoDimWeight;
    const double D = defaultWeight;
    const double heaviest = 9.;
    const double H = heaviest;
    const double N = lin3dlin3d;
    const static double weightTable[8][8] = {
        {L, L, L, L, L, L, H, s}, // 0 dim
        {L, L, L, L, L, L, H, s}, // 1 dim
        {L, L, q, q, q, q, H, s}, // 2 dim linear
        {L, L, q, N, q, D, H, s}, // 3 dim linear
        {L, L, q, q, q, q, H, s}, // 2 dim higher-order
        {L, L, q, D, q, D, H, s}, // 3 dim higher-order
        {H, H, H, H, H, H, H, s}, // misc heavy
        {s, s, s, s, s, s, s, s}  // super element        7
    };

    int element1Index = getEdgeWeightTableIndex(element1Topology);
    int element2Index = getEdgeWeightTableIndex(element2Topology);

    return weightTable[element1Index][element2Index];
}

double GraphCreationSettings::getGraphVertexWeight(stk::mesh::Entity entity, int criteria_index) const
{
    return 1.0;
}

int GraphCreationSettings::getGraphVertexWeight(stk::topology type) const
{
    switch(type)
    {
        case stk::topology::PARTICLE:
        case stk::topology::LINE_2:
        case stk::topology::BEAM_2:
        case stk::topology::BEAM_3:
        case stk::topology::SPRING_2:
        case stk::topology::SPRING_3:
            return 1;
        case stk::topology::SHELL_TRIANGLE_3:
            return 3;
        case stk::topology::SHELL_TRIANGLE_6:
            return 8;
        case stk::topology::SHELL_QUADRILATERAL_4:
            return 6;
        case stk::topology::SHELL_QUADRILATERAL_8:
            return 8;
        case stk::topology::TRI_3_2D:
            return 3;
        case stk::topology::QUAD_4_2D:
            return 6;
        case stk::topology::HEXAHEDRON_8:
            return 3;
        case stk::topology::HEXAHEDRON_20:
            return 8;
        case stk::topology::TETRAHEDRON_4:
            return 1;
        case stk::topology::TETRAHEDRON_10:
            return 6;
        case stk::topology::WEDGE_6:
            return 2;
        case stk::topology::WEDGE_12:
            // TODO
        case stk::topology::WEDGE_15:
            return 12;
        default:
            if ( type.is_superelement( ))
            {
                return 10;
            }
            throw("Invalid Element Type In WeightsOfElement");
    }
}

BalanceSettings::GraphOption GraphCreationSettings::getGraphOption() const
{
    return BalanceSettings::LOAD_BALANCE;
}

bool GraphCreationSettings::includeSearchResultsInGraph() const
{
    return m_includeSearchResultInGraph;
}

void GraphCreationSettings::setIncludeSearchResultsInGraph(bool doContactSearch)
{
    m_includeSearchResultInGraph = doContactSearch;
}

double GraphCreationSettings::getToleranceForParticleSearch() const
{
    return mToleranceForParticleSearch;
}

void GraphCreationSettings::setToleranceFunctionForFaceSearch(std::shared_ptr<stk::balance::FaceSearchTolerance> faceSearchTolerance)
{
    m_faceSearchToleranceFunction = faceSearchTolerance;
    m_UseConstantToleranceForFaceSearch = false;
}

bool GraphCreationSettings::isConstantFaceSearchTolerance() const
{
  return m_UseConstantToleranceForFaceSearch;
}

double GraphCreationSettings::getToleranceForFaceSearch(const stk::mesh::BulkData & mesh,
                                                        const stk::mesh::FieldBase & coordField,
                                                        const stk::mesh::Entity * faceNodes,
                                                        const unsigned numFaceNodes) const
{
    if (m_UseConstantToleranceForFaceSearch) {
        return mToleranceForFaceSearch;
    }
    else {
        return m_faceSearchToleranceFunction->compute(mesh, coordField, faceNodes, numFaceNodes);
    }
}

bool GraphCreationSettings::getEdgesForParticlesUsingSearch() const
{
    return false;
}

double GraphCreationSettings::getVertexWeightMultiplierForVertexInSearch() const
{
    return vertexWeightMultiplierForVertexInSearch;
}

std::string GraphCreationSettings::getDecompMethod() const
{
    return method;
}

void GraphCreationSettings::setDecompMethod(const std::string& input_method)
{
    method = input_method;
}

void GraphCreationSettings::setToleranceForFaceSearch(double tol)
{
    m_UseConstantToleranceForFaceSearch = true;
    mToleranceForFaceSearch = tol;
}
void GraphCreationSettings::setToleranceForParticleSearch(double tol)
{
    mToleranceForParticleSearch = tol;
}
void GraphCreationSettings::setEdgeWeightForSearch(double w)
{
    edgeWeightForSearch = w;
}
void GraphCreationSettings::setVertexWeightMultiplierForVertexInSearch(double w)
{
    vertexWeightMultiplierForVertexInSearch = w;
}
int GraphCreationSettings::getConnectionTableIndex(stk::topology elementTopology) const
{
    int tableIndex = -1;
    switch(elementTopology)
    {
        case stk::topology::PARTICLE:
            tableIndex = 0;
            break;
        case stk::topology::LINE_2:
        case stk::topology::LINE_2_1D:
        case stk::topology::LINE_3_1D:
        case stk::topology::BEAM_2:
        case stk::topology::BEAM_3:
        case stk::topology::SHELL_LINE_2:
        case stk::topology::SHELL_LINE_3:
        case stk::topology::SPRING_2:
        case stk::topology::SPRING_3:
            tableIndex = 1;
            break;
        case stk::topology::TRI_3_2D:
        case stk::topology::TRI_4_2D:
        case stk::topology::QUAD_4_2D:
        case stk::topology::SHELL_TRI_3:
        case stk::topology::SHELL_TRI_4:
        case stk::topology::SHELL_QUAD_4:
            tableIndex = 2;
            break;
        case stk::topology::TET_4:
        case stk::topology::PYRAMID_5:
        case stk::topology::WEDGE_6:
        case stk::topology::HEX_8:
            tableIndex = 3;
            break;
        case stk::topology::TRI_6_2D:
        case stk::topology::QUAD_8_2D:
        case stk::topology::QUAD_9_2D:
        case stk::topology::SHELL_TRI_6:
        case stk::topology::SHELL_QUAD_8:
        case stk::topology::SHELL_QUAD_9:
            tableIndex = 4;
            break;
        case stk::topology::TET_8:
        case stk::topology::TET_10:
        case stk::topology::TET_11:
        case stk::topology::PYRAMID_13:
        case stk::topology::PYRAMID_14:
        case stk::topology::WEDGE_15:
        case stk::topology::WEDGE_18:
        case stk::topology::HEX_20:
        case stk::topology::HEX_27:
            tableIndex = 5;
            break;
        default:
            if(elementTopology.is_superelement())
            {
                tableIndex = 6;
            }
            else
            {
                std::cerr << "Topology is " << elementTopology << std::endl;
                throw("Invalid Element Type in GetDimOfElement");
            }
            break;
    };
    return tableIndex;
}
int GraphCreationSettings::getEdgeWeightTableIndex(stk::topology elementTopology) const
{
    int tableIndex = -1;
    switch(elementTopology)
    {
        case stk::topology::PARTICLE:
            tableIndex = 0;
            break;
        case stk::topology::LINE_2:
        case stk::topology::LINE_2_1D:
        case stk::topology::LINE_3_1D:
        case stk::topology::BEAM_2:
        case stk::topology::BEAM_3:
        case stk::topology::SHELL_LINE_2:
        case stk::topology::SHELL_LINE_3:
        case stk::topology::SPRING_2:
        case stk::topology::SPRING_3:
            tableIndex = 1;
            break;
        case stk::topology::TRI_3_2D:
        case stk::topology::TRI_4_2D:
        case stk::topology::QUAD_4_2D:
        case stk::topology::SHELL_TRI_3:
        case stk::topology::SHELL_TRI_4:
        case stk::topology::SHELL_QUAD_4:
            tableIndex = 2;
            break;
        case stk::topology::TET_4:
        case stk::topology::PYRAMID_5:
        case stk::topology::WEDGE_6:
        case stk::topology::HEX_8:
            tableIndex = 3;
            break;
        case stk::topology::TRI_6_2D:
        case stk::topology::QUAD_8_2D:
        case stk::topology::QUAD_9_2D:
        case stk::topology::SHELL_TRI_6:
        case stk::topology::SHELL_QUAD_9:
            tableIndex = 4;
            break;
        case stk::topology::SHELL_QUAD_8:
        case stk::topology::TET_8:
        case stk::topology::TET_10:
        case stk::topology::TET_11:
        case stk::topology::PYRAMID_13:
        case stk::topology::PYRAMID_14:
        case stk::topology::WEDGE_12:
        case stk::topology::WEDGE_15:
        case stk::topology::WEDGE_18:
        case stk::topology::HEX_27:
            tableIndex = 5;
            break;
        case stk::topology::HEX_20:
            tableIndex = 6;
            break;
        default:
            if(elementTopology.is_superelement())
            {
                tableIndex = 7;
            }
            else
            {
                std::cerr << "Topology is " << elementTopology << std::endl;
                throw("Invalid Element Type in GetDimOfElement");
            }
            break;
    };
    return tableIndex;
}

void GraphCreationSettings::setShouldFixSpiders(bool fixSpiders)
{
    m_shouldFixSpiders = fixSpiders;
}

bool GraphCreationSettings::shouldFixMechanisms() const
{
    return true;
}

bool GraphCreationSettings::shouldFixSpiders() const
{
    return m_shouldFixSpiders;
}

const stk::mesh::Field<int> * GraphCreationSettings::getSpiderBeamConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const
{
    if (m_spiderBeamConnectivityCountField == nullptr) {
        m_spiderBeamConnectivityCountField =
            reinterpret_cast<stk::mesh::Field<int>*>(stkMeshBulkData.mesh_meta_data().get_field(stk::topology::NODE_RANK,
                                                     getSpiderBeamConnectivityCountFieldName()));
        ThrowRequireMsg(m_spiderBeamConnectivityCountField != nullptr,
                        "Must create nodal spider beam connectivity count field when stomping spiders.");
    }
    return m_spiderBeamConnectivityCountField;
}

const stk::mesh::Field<int> * GraphCreationSettings::getSpiderVolumeConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const
{
    if (m_spiderVolumeConnectivityCountField == nullptr) {
        m_spiderVolumeConnectivityCountField =
            reinterpret_cast<stk::mesh::Field<int>*>(stkMeshBulkData.mesh_meta_data().get_field(stk::topology::ELEM_RANK,
                                                     getSpiderVolumeConnectivityCountFieldName()));
        ThrowRequireMsg(m_spiderVolumeConnectivityCountField != nullptr,
                        "Must create element spider volume connectivity count field when stomping spiders.");
    }
    return m_spiderVolumeConnectivityCountField;
}

void GraphCreationSettings::setUseNodeBalancer(bool useBalancer)
{
  m_useNodeBalancer = useBalancer;
}

void GraphCreationSettings::setNodeBalancerTargetLoadBalance(double targetLoadBalance)
{
  m_nodeBalancerTargetLoadBalance = targetLoadBalance;
}

void GraphCreationSettings::setNodeBalancerMaxIterations(unsigned maxIterations)
{
  m_nodeBalancerMaxIterations = maxIterations;
}

bool GraphCreationSettings::useNodeBalancer() const
{
  return m_useNodeBalancer;
}

double GraphCreationSettings::getNodeBalancerTargetLoadBalance() const
{
  return m_nodeBalancerTargetLoadBalance;
}

unsigned GraphCreationSettings::getNodeBalancerMaxIterations() const
{
  return m_nodeBalancerMaxIterations;
}

const std::string& get_coloring_part_base_name()
{
    static std::string coloringPartBaseName = "STK_INTERNAL_COLORING_PART";
    return coloringPartBaseName;
}

stk::mesh::Part* get_coloring_part(const stk::mesh::BulkData& bulk, const stk::mesh::Entity& entity)
{
    const stk::mesh::Bucket& bucket = bulk.bucket(entity);
    const stk::mesh::PartVector& parts = bucket.supersets();
    stk::mesh::Part* colorPart = nullptr;
    unsigned numColors = 0;
    const std::string coloringBaseName = get_coloring_part_base_name();
    const unsigned length = coloringBaseName.length();
    for (stk::mesh::Part* part : parts)
    {
        std::string partSubName = part->name().substr(0, length);
        if (!sierra::case_strcmp(partSubName, coloringBaseName))
        {
            ++numColors;
            colorPart = part;
        }
    }
    ThrowRequireMsg(numColors <= 1, "Entity " << bulk.entity_key(entity) << " has " << numColors << " coloring parts.");
    return colorPart;
}

stk::mesh::PartVector get_root_topology_parts_for_rank(const stk::mesh::BulkData& bulk, stk::mesh::EntityRank rank)
{
  const stk::mesh::MetaData& meta = bulk.mesh_meta_data();
  std::set<stk::mesh::Part*> parts;
  stk::mesh::BucketVector buckets = bulk.buckets(rank);
  for (stk::mesh::Bucket* bucket : buckets)
  {
    stk::topology topo = bucket->topology();
    stk::mesh::Part& topoPart = meta.get_topology_root_part(topo);
    parts.insert(&topoPart);
  }

  return stk::mesh::PartVector(parts.begin(), parts.end());
}

}
}
