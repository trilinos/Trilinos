#include "MasterElementProviderIntrepid2.hpp"
#include <Intrepid2_CellTools.hpp>
#include <Teuchos_RCPStdSharedPtrConversions.hpp>
#include "stk_search_util/ParametricDistance.hpp"
#include "stk_search_util/Intrepid2_HasParametricDistance.hpp"
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>


namespace stk::search {

using I2CellTools = Intrepid2::CellTools<typename BasisExecSpace::device_type>;
using I2RealSpaceTools = Intrepid2::RealSpaceTools<typename BasisExecSpace::device_type>;

void MasterElementProviderIntrepid2::evaluate_field(const Topology& topo,
    const std::vector<double>& paramCoords,     // (numParamCoords)
    const unsigned numFieldComponents,
    const std::vector<double>& fieldData,       // (numFieldComponents x numTopologyNodes)
    std::vector<double>& result) const     // (numFieldComponents)
{
  evaluate_field(topo, 1, paramCoords, numFieldComponents, fieldData, result);
}

void MasterElementProviderIntrepid2::evaluate_field(const Topology& topo,
    const unsigned numEvalPoints,
    const std::vector<double>& paramCoords,     // (numParamCoords x numEvalPoints)
    const unsigned numFieldComponents,
    const std::vector<double>& fieldData,       // (numFieldComponents x numTopologyNodes)
    std::vector<double>& result) const     // (numFieldComponents x numEvalPoints)
{
  check_consistent_topology(topo);
  stk::topology stkTopo = topo.get_topology();
  const unsigned numNodes = stkTopo.num_nodes();
  unsigned numParCoords = num_parametric_coordinates(topo);
  if (stkTopo.is_shell())
  {
    numParCoords -= 1;  // find_parametric_coords returns 3 parametric coords, the 3rd
                        // one from in the normal direction.  Evaluate field only
                        // needs the first two
  }

  if (stkTopo != m_lastTopology || numEvalPoints != m_numEvalPoints) {
    m_basis = Teuchos::rcp(lagrangeBasisFactory(stkTopo, m_useCompositeTet10));
    m_lastTopology = stkTopo;
    m_numEvalPoints = numEvalPoints;
    m_basisVals = BasisViewType("evaluate_field_basis_vals", stkTopo.num_nodes(), numEvalPoints);
    m_paramCoords = PointViewType("int_point", numEvalPoints, numParCoords);
  }

  STK_ThrowRequireMsg(paramCoords.size() >= numParCoords*numEvalPoints,
                      "Insufficient length for parametric coordinates: " << paramCoords.size() << " Expected: " << numParCoords*numEvalPoints);
  STK_ThrowRequireMsg(fieldData.size() >= numNodes*numFieldComponents,
                        "Insufficient length for fieldData: " << fieldData.size() << " Expected: " << numNodes*numFieldComponents);

  for (unsigned p=0; p < numEvalPoints; ++p){
    for (unsigned d=0; d < stkTopo.dimension(); ++d)
    {
      m_paramCoords(p, d) = paramCoords[p*numParCoords + d];
    }
  }

  m_basis->getValues(m_basisVals, m_paramCoords, Intrepid2::OPERATOR_VALUE);

  result.clear();
  result.assign(numFieldComponents*numEvalPoints, 0.0);

  for (unsigned node=0; node < numNodes; ++node)
  {
    for (unsigned comp=0; comp < numFieldComponents; ++comp)
    {
      unsigned idx = comp*numNodes + node;
      for (unsigned outputPt=0; outputPt < numEvalPoints; ++outputPt)
      {
        result[outputPt*numFieldComponents + comp] += m_basisVals(node, outputPt) * fieldData[idx];
      }
    }
  }
}

void MasterElementProviderIntrepid2::nodal_field_data(const EntityKey& key,
    const Field& meField,
    unsigned& numFieldComponents,
    unsigned& numNodes,
    std::vector<double>& fieldData) const  // (numFieldComponents x numTopologyNodes)
{
  const stk::mesh::Entity entity = key;

  // Extract transposed field data
  const stk::mesh::FieldBase& field = *meField.get_field();
  const stk::mesh::BulkData& bulk = field.get_mesh();

  // load nodal coordinates from entity
  stk::mesh::Entity const* nodes = bulk.begin_nodes(entity);

  numNodes = bulk.num_nodes(entity);

  // Implicit assumption that number of field components is same for all nodes
  numFieldComponents = numNodes > 0 ? stk::mesh::field_extent0_per_entity(field, nodes[0]) : 0;

  fieldData.resize(static_cast<size_t>(numFieldComponents) * numNodes, 0.0);

  stk::search::CachedEntityFieldData<double> data;

  for(unsigned ni = 0u; ni < numNodes; ++ni) {
    stk::mesh::Entity node = nodes[ni];

    if(field.defined_on(node)) {
      meField.populate_entity_data(node, data);

      for(unsigned j(0); j < numFieldComponents; ++j) {
        const auto offSet = ni + j * numNodes;
        fieldData[offSet] = meField.transform(data.constPointer[j * data.componentStride]);
      }
    } else {
      for(unsigned j = 0; j < numFieldComponents; ++j) {
        const auto offSet = ni + j * numNodes;
        fieldData[offSet] = meField.default_value();
      }
    }
  }
}

void MasterElementProviderIntrepid2::nodal_field_data(const std::vector<EntityKey>& nodeKeys,
    const Field& meField,
    unsigned& numFieldComponents,
    std::vector<double>& fieldData) const  // (numFieldComponents x numTopologyNodes)
{
  // Extract transposed field data
  const stk::mesh::FieldBase& field = *meField.get_field();
  unsigned numNodes = nodeKeys.size();

  // Implicit assumption that number of field components is same for all nodes
  numFieldComponents = numNodes > 0 ? stk::mesh::field_extent0_per_entity(field, nodeKeys[0]) : 0;
  fieldData.resize(static_cast<size_t>(numFieldComponents) * numNodes, 0.0);

  stk::search::CachedEntityFieldData<double> data;

  for (unsigned ni = 0u; ni < numNodes; ++ni) {
    stk::mesh::Entity node = nodeKeys[ni];

    unsigned numNodeFieldComponents = stk::mesh::field_extent0_per_entity(field, node);

    if(field.defined_on(node)) {
      meField.populate_entity_data(node, data);

      for(int j(0); j < static_cast<int>(numNodeFieldComponents); ++j) {
        const auto offSet = ni + j * numNodes;
        fieldData[offSet] = meField.transform(data.constPointer[j * data.componentStride]);
      }
    } else {
      for(unsigned j = 0; j < numNodeFieldComponents; ++j) {
        const auto offSet = ni + j * numNodes;
        fieldData[offSet] = meField.default_value();
      }
    }
  }
}

//This has to be a templated function
template<bool Intrepid2HasNewFunction, typename ViewType, typename BasisType>
void compute_parametric_coords(unsigned numParamCoords, std::vector<double>& paramCoords,
                               ViewType& refCellCenter, ViewType& physCoordsView,
                               ViewType& elemNodesView, BasisType& basis)
{
  constexpr unsigned numCells=1, numPoints=1;
  ViewType paramCoordsView("param_coords_view", numCells, numPoints, numParamCoords);

  if constexpr (Intrepid2HasNewFunction) {
    I2RealSpaceTools::clone(paramCoordsView, refCellCenter);
    I2CellTools::mapToReferenceFrame(paramCoordsView, physCoordsView, elemNodesView, basis);
  }
  else {
    ViewType initGuess("init_guess", numCells, numPoints, numParamCoords);
    I2RealSpaceTools::clone(initGuess, refCellCenter);
    I2CellTools::mapToReferenceFrameInitGuess(paramCoordsView, initGuess, physCoordsView, elemNodesView, basis);
  }

  paramCoords.resize(numParamCoords);
  for (unsigned d=0; d < numParamCoords; ++d) {
    paramCoords[d] = paramCoordsView(0, 0, d);
  }
}

void MasterElementProviderIntrepid2::find_parametric_coordinates(const Topology& topo,
    const unsigned numCoordComponents,
    const std::vector<double>& elemNodeCoords,  // (numCoordComponents x numTopologyNodes)
    const std::vector<double>& inputCoords,     // (numCoordComponents)
    std::vector<double>& paramCoords,
    double& paramDistance) const
{
  check_consistent_topology(topo);
  using ViewType = Kokkos::DynRankView<double, typename BasisExecSpace::device_type>;

  stk::topology stkTopo = topo.get_topology();
  constexpr unsigned numCells=1, numPoints=1;
  unsigned numNodes = stkTopo.num_nodes();
  unsigned spatialDim = numCoordComponents;

  STK_ThrowRequireMsg(inputCoords.size() >= numCoordComponents,
                      "Insufficient length for inputCoords: " << inputCoords.size() << " Expected: " << numCoordComponents);

  STK_ThrowRequireMsg(elemNodeCoords.size() >= stkTopo.num_nodes()*numCoordComponents,
                      "Insufficient length for elemNodeCoords: " << elemNodeCoords.size() << " Expected: " <<
                      stkTopo.num_nodes()*numCoordComponents);

  ViewType physCoordsView("phys_coords_view", numCells, numPoints, spatialDim);
  ViewType elemNodesView("elem_nodes_view", numCells, numNodes, spatialDim);

  for (unsigned d=0; d < spatialDim; ++d) {
    physCoordsView(0, 0, d) = inputCoords[d];
  }

  for (unsigned node=0; node < stkTopo.num_nodes(); ++node) {
    for (unsigned d=0; d < spatialDim; ++d) {
      elemNodesView(0, node, d) = elemNodeCoords[d*stkTopo.num_nodes() + node];
    }
  }

  if (stkTopo != m_lastTopology) {
    m_basis = Teuchos::rcp(lagrangeBasisFactory(stkTopo, m_useCompositeTet10));
    m_lastTopology = stkTopo;
  }

  auto shardsCellTopo = stk::mesh::get_cell_topology(stkTopo);

  ViewType refCellCenter("ref_cell_center", spatialDim);
  I2CellTools::getReferenceCellCenter(refCellCenter, shardsCellTopo);

  const unsigned numParCoords = num_parametric_coordinates(topo);
  compute_parametric_coords<stk::search::intrepid2_has_param_dist>(numParCoords, paramCoords, refCellCenter, physCoordsView, elemNodesView, m_basis);

  paramDistance = compute_parametric_distance(stkTopo, paramCoords);
}

void MasterElementProviderIntrepid2::coordinate_center(const Topology& topo, std::vector<double>& coords) const
{
  check_consistent_topology(topo);

  using ViewType = Kokkos::DynRankView<double, typename BasisExecSpace::memory_space>;

  stk::topology stkTopo = topo.get_topology();
  coords.resize(stkTopo.dimension(), 0.0);

  ViewType refCellCenterView("coordinate_center_view", stkTopo.dimension());
  I2CellTools::getReferenceCellCenter(refCellCenterView, stk::mesh::get_cell_topology(topo));

  for (unsigned d=0; d < stkTopo.dimension(); ++d)
  {
    coords[d] = refCellCenterView(0, 0, d);
  }
}

unsigned MasterElementProviderIntrepid2::num_parametric_coordinates(const Topology& topo) const
{
  check_consistent_topology(topo);
  return topo.get_topology().dimension();
}

unsigned MasterElementProviderIntrepid2::num_integration_points(const Topology& topo) const
{
  check_consistent_topology(topo);
  return get_default_cubature(topo)->getNumPoints();
}

constexpr unsigned MAX_NUM_INTG_PTS = 27;
std::array<int,MAX_NUM_INTG_PTS> get_permutation(stk::topology topo)
{
  std::array<int,MAX_NUM_INTG_PTS> perm;
  switch(topo) {
  case stk::topology::HEX_8:
    perm = {7,3,5,1,6,2,4,0};
    break;
  case stk::topology::HEX_20:
    perm = {26,17,8,23,14,5,20,11,2,25,16,7,22,13,4,19,10,1,24,15,6,21,12,3,18,9,0};
    break;
  case stk::topology::WEDGE_6:
    perm = {5,4,3,2,1,0};
    break;
  case stk::topology::PYRAMID_5:
    perm = {7,5,6,4,3,1,2,0};
    break;
  default:
    perm = {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26};
    break;
  }
  return perm;
}

void MasterElementProviderIntrepid2::integration_points(const Topology& topo, std::vector<double>& gaussPoints) const
{
  check_consistent_topology(topo);
  stk::topology stkTopo = topo.get_topology();
  Teuchos::RCP<Cubature> cubature = get_default_cubature(topo);

  const unsigned dim = stkTopo.dimension();
  PointViewType intPoints("int_point", cubature->getNumPoints(), dim);
  WeightViewType intWeights("int_weights", cubature->getNumPoints());
  cubature->getCubature(intPoints, intWeights);

  const auto perm = get_permutation(stkTopo);

  gaussPoints.resize(cubature->getNumPoints()*dim);
  for (unsigned i=0; i < static_cast<unsigned>(cubature->getNumPoints()); ++i)
  {
    for (unsigned d=0; d < dim; ++d)
    {
      gaussPoints[i*dim + d] = intPoints(perm[i], d);
    }
  }
}

using PDViewType = Kokkos::DynRankView<double, typename BasisExecSpace::memory_space>;

template<unsigned CellTopoKey>
double intrepid2_param_dist(const PDViewType& paramPoint)
{
  double paramDist = 999999.9;
  if constexpr ((stk::search::intrepid2_has_param_dist)) {
    paramDist = Intrepid2::ParametricDistance<CellTopoKey>::distance(paramPoint);
  }
  else {
    STK_ThrowErrorMsg("Intrepid2 doesn't support shells yet");
  }

  return paramDist;
}

double MasterElementProviderIntrepid2::compute_parametric_distance(const stk::topology& topo, const std::vector<double>& paramCoords) const
{
  check_consistent_topology(topo);
  // compute a distance d in the range [0, inft] where
  // 0 <= d <= 1 means the point is inside the element
  // otherwise the point is outside the element
  double parametricDistance = 0.0;

  PDViewType paramPoint("paramPoint", paramCoords.size());
  for(unsigned i=0; i<paramCoords.size(); ++i) {
    paramPoint(i) = paramCoords[i];
  }

  switch (topo) {
    case stk::topology::HEX_8:
    case stk::topology::HEX_20:
      parametricDistance = parametric_distance_hex(paramCoords);
      break;
    case stk::topology::TRI_3_2D:
      parametricDistance = parametric_distance_tri_2d(paramCoords);
      break;
    case stk::topology::TET_4:
    case stk::topology::TET_10:
      parametricDistance = parametric_distance_tet(paramCoords);
      break;
    case stk::topology::PYRAMID_5:
      parametricDistance = parametric_distance_pyramid(paramCoords);
      break;
    case stk::topology::WEDGE_6:
      parametricDistance = parametric_distance_wedge(paramCoords);
      break;
    case stk::topology::QUAD_4_2D:
      parametricDistance = parametric_distance_quad_2d(paramCoords);
      break;
    case stk::topology::SHELL_QUAD_4:
      parametricDistance = intrepid2_param_dist<shards::ShellQuadrilateral<>::key>(paramPoint);
      break;
    case stk::topology::SHELL_TRI_3:
      parametricDistance = intrepid2_param_dist<shards::ShellTriangle<>::key>(paramPoint);
      break;

    case stk::topology::LINE_2:
    default:
      STK_ThrowErrorMsg("Topology " + topo.name() + " not supported for compute_parametric_distance.");
  }

  return parametricDistance;
}

Teuchos::RCP<MasterElementProviderIntrepid2::Cubature>
MasterElementProviderIntrepid2::get_default_cubature(const Topology& topo) const
{
  stk::topology stkTopo = topo.get_topology();
  shards::CellTopology cellTopo = stk::mesh::get_cell_topology(stkTopo);
  STK_ThrowRequireMsg(cellTopo.isValid(),"stkTopo " << stkTopo << " doesn't map to valid Shards Cell Topology");

  Intrepid2::DefaultCubatureFactory cubatureFactory;
  const bool linearElement = stkTopo.num_vertices() == stkTopo.num_nodes();
  unsigned degree = linearElement ? 2 : 4;
  if (stkTopo == stk::topology::TET_10) {
      degree = 2;
  }

  const auto cellCubature = cubatureFactory.create<BasisExecSpace,double,double>(cellTopo, degree);

  return cellCubature;
}

void MasterElementProviderIntrepid2::check_consistent_topology(const Topology& topo) const
{
  stk::topology stkTopo = topo.get_topology();
  STK_ThrowRequireMsg(stkTopo == stk::topology::QUAD_4_2D ||
                      stkTopo == stk::topology::TRI_3_2D ||
                      stkTopo == stk::topology::HEX_8 ||
                      stkTopo == stk::topology::TET_4 ||
                      stkTopo == stk::topology::TET_10 ||
                      stkTopo == stk::topology::PYRAMID_5 ||
                      stkTopo == stk::topology::WEDGE_6 ||
                      stkTopo == stk::topology::HEX_20 ||
                      stkTopo == stk::topology::SHELL_QUAD_4 ||
                      stkTopo == stk::topology::SHELL_TRI_3,
                     "Invalid topology " << topo.get_topology());
}

}
