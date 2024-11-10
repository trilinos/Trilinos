#include "balanceUtils.hpp"
#include "mpi.h"
#include "search_tolerance/FaceSearchTolerance.hpp"
#include "stk_balance/search_tolerance_algs/SecondShortestEdgeFaceSearchTolerance.hpp"
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
  : m_numInputProcessors(0),
    m_numOutputProcessors(0),
    m_isRebalancing(false),
    m_shouldFixCoincidentElements(true),
    m_initialDecompMethod("RIB"),
    m_numCriteria(1),
    m_defaultFieldWeight(0.0),
    m_vertexWeightFieldNames{""},
    m_useNestedDecomp(false),
    m_shouldPrintDiagnostics(false),
    m_diagnosticElementWeightsField(nullptr),
    m_vertexConnectivityWeightField(nullptr),
    m_vertexWeightFields{nullptr},
    m_vertexWeightMethod(DefaultSettings::vertexWeightMethod),
    m_graphEdgeWeightMultiplier(DefaultSettings::graphEdgeWeightMultiplier)
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

double BalanceSettings::getFieldVertexWeight(const stk::mesh::BulkData &bulkData, stk::mesh::Entity entity, int criteria_index) const
{
    const stk::mesh::Field<double> &field = *getVertexWeightField(bulkData, criteria_index);
    const double *weight = stk::mesh::field_data(field, entity);
    if (weight != nullptr) {
      STK_ThrowRequireWithSierraHelpMsg(*weight >= 0);
      return *weight;
    }
    else {
      return m_defaultFieldWeight;
    }
}

BalanceSettings::GraphOption BalanceSettings::getGraphOption() const
{
  return BalanceSettings::LOAD_BALANCE;
}

double BalanceSettings::getGraphEdgeWeightMultiplier() const
{
  return m_graphEdgeWeightMultiplier;
}

void BalanceSettings::setGraphEdgeWeightMultiplier(double multiplier)
{
  m_graphEdgeWeightMultiplier = multiplier;
}

void BalanceSettings::setVertexWeightMethod(VertexWeightMethod method)
{
  m_vertexWeightMethod = method;
}

VertexWeightMethod BalanceSettings::getVertexWeightMethod() const
{
  return m_vertexWeightMethod;
}

bool
BalanceSettings::shouldFixCoincidentElements() const
{
  return m_shouldFixCoincidentElements;
}

void
BalanceSettings::setShouldFixCoincidentElements(bool fixCoincidentElements)
{
  m_shouldFixCoincidentElements = fixCoincidentElements;
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

void BalanceSettings::setVertexWeightBlockMultiplier(const std::string & blockName, double multiplier)
{
  m_vertexWeightBlockMultipliers[blockName] = multiplier;
}

const BlockWeightMultipliers &
BalanceSettings::getVertexWeightBlockMultipliers() const
{
  return m_vertexWeightBlockMultipliers;
}

bool BalanceSettings::isIncrementalRebalance() const
{
  return false;
}

bool BalanceSettings::isMultiCriteriaRebalance() const
{
  return false;
}

void BalanceSettings::setVertexWeightFieldName(std::string field_name, unsigned criteria_index)
{
  STK_ThrowRequireMsg(criteria_index < m_vertexWeightFieldNames.size(),
                  "The provided criteria index (" + std::to_string(criteria_index) + ") is too large for the " +
                  "supported number of criteria (" + std::to_string(m_vertexWeightFieldNames.size()) + ")");
  m_vertexWeightFieldNames[criteria_index] = field_name;
}

std::string BalanceSettings::getVertexWeightFieldName(unsigned criteria_index) const
{ 
  STK_ThrowRequireMsg(criteria_index < m_vertexWeightFieldNames.size(),
                  "The provided criteria index (" + std::to_string(criteria_index) + ") is too large for the " +
                  "supported number of criteria (" + std::to_string(m_vertexWeightFieldNames.size()) + ")");
  return m_vertexWeightFieldNames[criteria_index]; 
}

const stk::mesh::Field<double> * BalanceSettings::getVertexWeightField(const stk::mesh::BulkData & stkMeshBulkData, unsigned criteria_index) const
{
  STK_ThrowRequireMsg(criteria_index < m_vertexWeightFieldNames.size(),
                  "The provided criteria index (" + std::to_string(criteria_index) + ") is too large for the " +
                  "supported number of criteria (" + std::to_string(m_vertexWeightFieldNames.size()) + ")");
  if (m_vertexWeightFields[criteria_index] == nullptr) {
    m_vertexWeightFields[criteria_index] =
        stkMeshBulkData.mesh_meta_data().get_field<double>(stk::topology::ELEM_RANK,
                                                           getVertexWeightFieldName(criteria_index));
    STK_ThrowRequireMsg(m_vertexWeightFields[criteria_index] != nullptr,
                    "Must provide a field for criteria index (" + std::to_string(criteria_index) + ")");
  }
  return m_vertexWeightFields[criteria_index];
}

void BalanceSettings::setDefaultFieldWeight(double defaultFieldWeight)
{
  m_defaultFieldWeight = defaultFieldWeight;
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

void BalanceSettings::setShouldPrintDiagnostics(bool shouldPrint)
{
  m_shouldPrintDiagnostics = shouldPrint;
}

bool BalanceSettings::shouldPrintDiagnostics() const
{
  return m_shouldPrintDiagnostics;
}

bool BalanceSettings::shouldPrintMetrics() const
{
  return false;
}

void BalanceSettings::setNumCriteria(int num_criteria)
{ 
  m_numCriteria = num_criteria;
  m_vertexWeightFieldNames.resize(m_numCriteria); 
  m_vertexWeightFields.resize(m_numCriteria); 
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

// For graph based methods only
bool BalanceSettings::allowModificationOfVertexWeightsForSmallMeshes() const
{
  return true;
}

// For graph based methods only
bool BalanceSettings::shouldFixMechanisms() const
{
  return false;
}

bool BalanceSettings::shouldFixSpiders() const
{
  return false;
}

std::string BalanceSettings::getSpiderPartName() const
{
  return "stk_balance_spider_elements";
}

std::string BalanceSettings::getSpiderVolumeConnectivityCountFieldName() const
{
  return "stk_balance_volume_connectivity_count";
}

std::string BalanceSettings::getOutputSubdomainFieldName() const
{
  return "stk_balance_output_subdomain";
}

std::string BalanceSettings::getDiagnosticElementWeightFieldName() const
{

  return "stk_balance_diagnostic_element_weight" + std::to_string(getNumCriteria());
}

std::string BalanceSettings::getVertexConnectivityWeightFieldName() const
{
  return "stk_balance_vertex_connectivity_weight";
}

stk::mesh::Part * BalanceSettings::getSpiderPart(const stk::mesh::BulkData & stkMeshBulkData) const
{
  return nullptr;
}

const stk::mesh::Field<int> * BalanceSettings::getSpiderVolumeConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const
{
  return nullptr;
}

const stk::mesh::Field<int> * BalanceSettings::getOutputSubdomainField(const stk::mesh::BulkData & stkMeshBulkData) const
{
  return nullptr;
}

const stk::mesh::Field<double> * BalanceSettings::getDiagnosticElementWeightField(const stk::mesh::BulkData & stkMeshBulkData) const
{
  if (m_diagnosticElementWeightsField == nullptr ||
      (getDiagnosticElementWeightFieldName() != m_diagnosticElementWeightsField->name())) {
    m_diagnosticElementWeightsField =
        stkMeshBulkData.mesh_meta_data().get_field<double>(stk::topology::ELEM_RANK,
                                                           getDiagnosticElementWeightFieldName());
    STK_ThrowRequireMsg(m_diagnosticElementWeightsField != nullptr,
                    "Must create diagnostic element weight field when printing balance diagnostics.");
  }
  return m_diagnosticElementWeightsField;
}

const stk::mesh::Field<double> * BalanceSettings::getVertexConnectivityWeightField(const stk::mesh::BulkData & stkMeshBulkData) const
{
  if (m_vertexConnectivityWeightField == nullptr) {
    m_vertexConnectivityWeightField =
        stkMeshBulkData.mesh_meta_data().get_field<double>(stk::topology::ELEM_RANK,
                                                           getVertexConnectivityWeightFieldName());
    STK_ThrowRequireMsg(m_vertexConnectivityWeightField != nullptr,
                    "Must create vertex connectivity weight field when printing balance diagnostics.");
  }
  return m_vertexConnectivityWeightField;
}


bool BalanceSettings::usingColoring() const
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

void BalanceSettings::set_log_filename(const std::string& filename)
{
  m_logFilename = filename;
}

std::string BalanceSettings::get_log_filename() const
{
  return m_logFilename;
}

//////////////////////////////////////

GraphCreationSettings::GraphCreationSettings()
  : m_method(DefaultSettings::decompMethod),
    m_ToleranceForFaceSearch(DefaultSettings::faceSearchAbsTol),
    m_ToleranceForParticleSearch(DefaultSettings::particleSearchTol),
    m_vertexWeightMultiplierForVertexInSearch(DefaultSettings::faceSearchVertexMultiplier),
    m_edgeWeightForSearch(DefaultSettings::faceSearchEdgeWeight),
    m_UseConstantToleranceForFaceSearch(false),
    m_shouldFixSpiders(DefaultSettings::fixSpiders),
    m_shouldFixMechanisms(DefaultSettings::fixMechanisms),
    m_spiderPart(nullptr),
    m_spiderVolumeConnectivityCountField(nullptr),
    m_outputSubdomainField(nullptr),
    m_includeSearchResultInGraph(DefaultSettings::useContactSearch),
    m_useNodeBalancer(false),
    m_nodeBalancerTargetLoadBalance(1.0),
    m_nodeBalancerMaxIterations(5)
{
  setToleranceFunctionForFaceSearch(
      std::make_shared<stk::balance::SecondShortestEdgeFaceSearchTolerance>(DefaultSettings::faceSearchRelTol)
  );
}

GraphCreationSettings::GraphCreationSettings(double faceSearchTol, double particleSearchTol, double edgeWeightSearch,
                                             const std::string& decompMethod, double multiplierVWSearch)
  : GraphCreationSettings()
{
  m_method = decompMethod;
  m_ToleranceForFaceSearch = faceSearchTol;
  m_ToleranceForParticleSearch = particleSearchTol;
  m_vertexWeightMultiplierForVertexInSearch = multiplierVWSearch;
  m_edgeWeightForSearch = edgeWeightSearch;
}

size_t GraphCreationSettings::getNumNodesRequiredForConnection(stk::topology element1Topology, stk::topology element2Topology) const
{
  const int noConnection = 1000;
  const int s = noConnection;
  const static int connectionTable[7][7] = {
    {1, 1, 1, 1, 1, 1, s}, // 0 dim
    {1, 1, 1, 1, 1, 1, s}, // 1 dim
    {1, 1, 2, 2, 2, 2, s}, // 2 dim linear
    {1, 1, 2, 3, 3, 3, s}, // 3 dim linear
    {1, 1, 2, 3, 3, 3, s}, // 2 dim higher-order
    {1, 1, 2, 3, 3, 4, s}, // 3 dim higher-order
    {s, s, s, s, s, s, s}  // super element
  };

  int element1Index = getConnectionTableIndex(element1Topology);
  int element2Index = getConnectionTableIndex(element2Topology);

  return connectionTable[element1Index][element2Index];
}

double GraphCreationSettings::getGraphEdgeWeightForSearch() const
{
  return m_edgeWeightForSearch;
}

double GraphCreationSettings::getGraphEdgeWeight(stk::topology element1Topology, stk::topology element2Topology) const
{
  const double noConnection = 0;
  const double A = noConnection; // super element
  const double defaultWeight = 1.0;
  const double linTetlinTet = defaultWeight;
  const double linTetlinHex = defaultWeight;
  const double linHexlinHex = defaultWeight;
  const double B = linTetlinTet;
  const double C = linTetlinHex;
  const double D = linHexlinHex;
  const double E = defaultWeight;
  const double twoDimWeight = 5;
  const double F = twoDimWeight;
  const double largeWeight = 5;  
  const double G = largeWeight;
  const double H = 9.;

  const static double weightTable[9][9] = {
    {G, G, G, G, G, G, G, H, A}, // 0 dim
    {G, G, G, G, G, G, G, H, A}, // 1 dim
    {G, G, F, F, F, F, F, H, A}, // 2 dim linear
    {G, G, F, B, C, F, E, H, A}, // 3 dim linear tet
    {G, G, F, C, D, F, E, H, A}, // 3 dim linear hex
    {G, G, F, F, F, F, F, H, A}, // 2 dim higher-order
    {G, G, F, E, E, F, E, H, A}, // 3 dim higher-order
    {H, H, H, H, H, H, H, H, A}, // misc heavy
    {A, A, A, A, A, A, A, A, A}  // super element
  };

  int element1Index = getEdgeWeightTableIndex(element1Topology);
  int element2Index = getEdgeWeightTableIndex(element2Topology);

  return weightTable[element1Index][element2Index];
}

int GraphCreationSettings::getGraphVertexWeight(stk::topology type) const
{
  switch(type)
  {
  case stk::topology::PARTICLE:
  case stk::topology::LINE_2_1D:
  case stk::topology::LINE_3_1D:
  case stk::topology::BEAM_2:
  case stk::topology::BEAM_3:
  case stk::topology::SHELL_LINE_2:
  case stk::topology::SHELL_LINE_3:
  case stk::topology::SPRING_2:
  case stk::topology::SPRING_3:
    return 1;
  case stk::topology::TRI_3_2D:
    return 3;
  case stk::topology::TRI_4_2D:
    return 4;
  case stk::topology::TRI_6_2D:
    return 6;
  case stk::topology::QUAD_4_2D:
    return 6;
  case stk::topology::QUAD_8_2D:
    return 8;
  case stk::topology::QUAD_9_2D:
    return 9;
  case stk::topology::SHELL_TRI_3:
    return 3;
  case stk::topology::SHELL_TRI_4:
    return 4;
  case stk::topology::SHELL_TRI_6:
    return 8;
  case stk::topology::SHELL_TRI_3_ALL_FACE_SIDES:
    return 3;
  case stk::topology::SHELL_TRI_4_ALL_FACE_SIDES:
    return 4;
  case stk::topology::SHELL_TRI_6_ALL_FACE_SIDES:
    return 8;
  case stk::topology::SHELL_QUAD_4:
    return 6;
  case stk::topology::SHELL_QUAD_8:
    return 8;
  case stk::topology::SHELL_QUAD_9:
    return 9;
  case stk::topology::SHELL_QUAD_4_ALL_FACE_SIDES:
    return 6;
  case stk::topology::SHELL_QUAD_8_ALL_FACE_SIDES:
    return 8;
  case stk::topology::SHELL_QUAD_9_ALL_FACE_SIDES:
    return 9;
  case stk::topology::TET_4:
    return 1;
  case stk::topology::TET_8:
    return 5;
  case stk::topology::TET_10:
    return 6;
  case stk::topology::TET_11:
    return 7;
  case stk::topology::PYRAMID_5:
    return 1;
  case stk::topology::PYRAMID_13:
    return 5;
  case stk::topology::PYRAMID_14:
    return 6;
  case stk::topology::WEDGE_6:
    return 2;
  case stk::topology::WEDGE_12:
    return 10;
  case stk::topology::WEDGE_15:
    return 12;
  case stk::topology::WEDGE_18:
    return 14;
  case stk::topology::HEX_8:
    return 3;
  case stk::topology::HEX_20:
    return 8;
  case stk::topology::HEX_27:
    return 10;
  default:
    if ( type.is_superelement( ))
    {
      return 10;
    }
    STK_ThrowErrorMsg("Unrecognized element type (" << type << ") in getGraphVertexWeight()");
    return 0;  // Keep compiler happy about always having a return value
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
  return m_ToleranceForParticleSearch;
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
    return m_ToleranceForFaceSearch;
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
  return m_vertexWeightMultiplierForVertexInSearch;
}

std::string GraphCreationSettings::getDecompMethod() const
{
  return m_method;
}

void GraphCreationSettings::setDecompMethod(const std::string& input_method)
{
  m_method = input_method;
}

void GraphCreationSettings::setToleranceForFaceSearch(double tol)
{
  m_UseConstantToleranceForFaceSearch = true;
  m_ToleranceForFaceSearch = tol;
}
void GraphCreationSettings::setToleranceForParticleSearch(double tol)
{
  m_ToleranceForParticleSearch = tol;
}
void GraphCreationSettings::setEdgeWeightForSearch(double w)
{
  m_edgeWeightForSearch = w;
}
void GraphCreationSettings::setVertexWeightMultiplierForVertexInSearch(double w)
{
  m_vertexWeightMultiplierForVertexInSearch = w;
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
  case stk::topology::SHELL_TRI_3_ALL_FACE_SIDES:
  case stk::topology::SHELL_TRI_4:
  case stk::topology::SHELL_TRI_4_ALL_FACE_SIDES:
  case stk::topology::SHELL_QUAD_4:
  case stk::topology::SHELL_QUAD_4_ALL_FACE_SIDES:
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
  case stk::topology::SHELL_TRI_6_ALL_FACE_SIDES:
  case stk::topology::SHELL_QUAD_8:
  case stk::topology::SHELL_QUAD_8_ALL_FACE_SIDES:
  case stk::topology::SHELL_QUAD_9:
  case stk::topology::SHELL_QUAD_9_ALL_FACE_SIDES:
    tableIndex = 4;
    break;
  case stk::topology::TET_8:
  case stk::topology::TET_10:
  case stk::topology::TET_11:
  case stk::topology::PYRAMID_13:
  case stk::topology::PYRAMID_14:
  case stk::topology::WEDGE_12:
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
      STK_ThrowErrorMsg("Unrecognized element type (" << elementTopology << ") in getConnectionTableIndex()");
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
  case stk::topology::SHELL_TRI_3_ALL_FACE_SIDES:
  case stk::topology::SHELL_TRI_4:
  case stk::topology::SHELL_TRI_4_ALL_FACE_SIDES:
  case stk::topology::SHELL_QUAD_4:
  case stk::topology::SHELL_QUAD_4_ALL_FACE_SIDES:
    tableIndex = 2;
    break;
  case stk::topology::TET_4:
  case stk::topology::PYRAMID_5:
  case stk::topology::WEDGE_6:
    tableIndex = 3;
    break;
  case stk::topology::HEX_8:
    tableIndex = 4;
    break;
  case stk::topology::TRI_6_2D:
  case stk::topology::QUAD_8_2D:
  case stk::topology::QUAD_9_2D:
  case stk::topology::SHELL_TRI_6:
  case stk::topology::SHELL_TRI_6_ALL_FACE_SIDES:
  case stk::topology::SHELL_QUAD_9:
  case stk::topology::SHELL_QUAD_9_ALL_FACE_SIDES:
    tableIndex = 5;
    break;
  case stk::topology::SHELL_QUAD_8:
  case stk::topology::SHELL_QUAD_8_ALL_FACE_SIDES:
  case stk::topology::TET_8:
  case stk::topology::TET_10:
  case stk::topology::TET_11:
  case stk::topology::PYRAMID_13:
  case stk::topology::PYRAMID_14:
  case stk::topology::WEDGE_12:
  case stk::topology::WEDGE_15:
  case stk::topology::WEDGE_18:
  case stk::topology::HEX_27:
    tableIndex = 6;
    break;
  case stk::topology::HEX_20:
    tableIndex = 7;
    break;
  default:
    if(elementTopology.is_superelement())
    {
      tableIndex = 8;
    }
    else
    {
      STK_ThrowErrorMsg("Unrecognized element type (" << elementTopology << ") in getEdgeWeightTableIndex()");
    }
    break;
  };
  return tableIndex;
}

void GraphCreationSettings::setShouldFixSpiders(bool fixSpiders)
{
  m_shouldFixSpiders = fixSpiders;
}

void GraphCreationSettings::setShouldFixMechanisms(bool fixMechanisms)
{
  m_shouldFixMechanisms = fixMechanisms;
}

bool GraphCreationSettings::shouldFixMechanisms() const
{
  return m_shouldFixMechanisms;
}

bool GraphCreationSettings::shouldFixSpiders() const
{
  return m_shouldFixSpiders;
}

stk::mesh::Part * GraphCreationSettings::getSpiderPart(const stk::mesh::BulkData & stkMeshBulkData) const
{
  if (m_spiderPart == nullptr) {
    m_spiderPart = stkMeshBulkData.mesh_meta_data().get_part(getSpiderPartName());
    STK_ThrowRequireMsg(m_spiderPart != nullptr,
                    "Must create spider part when fixing spider elements.");
  }
  return m_spiderPart;
}

const stk::mesh::Field<int> * GraphCreationSettings::getSpiderVolumeConnectivityCountField(const stk::mesh::BulkData & stkMeshBulkData) const
{
  if (m_spiderVolumeConnectivityCountField == nullptr) {
    m_spiderVolumeConnectivityCountField =
        stkMeshBulkData.mesh_meta_data().get_field<int>(stk::topology::ELEM_RANK,
                                                        getSpiderVolumeConnectivityCountFieldName());
    STK_ThrowRequireMsg(m_spiderVolumeConnectivityCountField != nullptr,
                    "Must create element spider volume connectivity count field when fixing spider elements.");
  }
  return m_spiderVolumeConnectivityCountField;
}

const stk::mesh::Field<int> * GraphCreationSettings::getOutputSubdomainField(const stk::mesh::BulkData & stkMeshBulkData) const
{
  if (m_outputSubdomainField == nullptr) {
    m_outputSubdomainField =
        stkMeshBulkData.mesh_meta_data().get_field<int>(stk::topology::ELEM_RANK,
                                                        getOutputSubdomainFieldName());
    STK_ThrowRequireMsg(m_outputSubdomainField != nullptr,
                    "Must create output subdomain field when fixing spider elements.");
  }
  return m_outputSubdomainField;
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
  STK_ThrowRequireMsg(numColors <= 1, "Entity " << bulk.entity_key(entity) << " has " << numColors << " coloring parts.");
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
