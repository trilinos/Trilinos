#include <Teuchos_GlobalMPISession.hpp>
#include <Teuchos_TimeMonitor.hpp>

#include "Teuchos_UnitTestRepository.hpp"

#include "Kokkos_Core.hpp"
#include "VerificationDriverHelpers.hpp"

#include <fstream>


#include <Intrepid2_DefaultCubatureFactory.hpp>
#include <Intrepid2_HGRAD_TET_Cn_FEM.hpp>
#include <Intrepid2_HGRAD_HEX_Cn_FEM.hpp>

#include <Intrepid2_HGRAD_LINE_Cn_FEM.hpp>
#include <Intrepid2_HVOL_LINE_Cn_FEM.hpp>

#include <Intrepid2_HGRAD_QUAD_Cn_FEM.hpp>
#include <Intrepid2_HCURL_QUAD_In_FEM.hpp>
#include <Intrepid2_HDIV_QUAD_In_FEM.hpp>
#include <Intrepid2_HVOL_QUAD_Cn_FEM.hpp>

#include <Intrepid2_HGRAD_HEX_Cn_FEM.hpp>
#include <Intrepid2_HCURL_HEX_In_FEM.hpp>
#include <Intrepid2_HDIV_HEX_In_FEM.hpp>
#include <Intrepid2_HVOL_HEX_Cn_FEM.hpp>

#include <Intrepid2_CellTools.hpp>
#include <Intrepid2_FunctionSpaceTools.hpp>

#include "Intrepid2_DerivedBasisFamily.hpp"
#include "Intrepid2_HierarchicalBasisFamily.hpp"
#include "Intrepid2_NodalBasisFamily.hpp"

#include "Intrepid2_IntegratedLegendreBasis_HGRAD_LINE.hpp"
#include "Intrepid2_LegendreBasis_HVOL_LINE.hpp"
#include "Intrepid2_Polynomials.hpp"
#include "Intrepid2_TensorBasis.hpp"

#include "Intrepid2_ESEAS_Interface.hpp"
#include "Intrepid2_TestUtils.hpp"

#include "Intrepid2_IntegratedLegendreBasis_HGRAD_TET.hpp"

using namespace Teuchos;
using namespace Kokkos;
using namespace Intrepid2;

using DeviceType = Kokkos::DefaultExecutionSpace::device_type;

bool testSegment(Intrepid2::EFunctionSpace fs, bool defineVertexFunctions)
{
  using namespace Intrepid2;
  using BasisFamilyCG = HierarchicalBasisFamily<DeviceType>;
  using BasisFamilyDG = DGHierarchicalBasisFamily<DeviceType>;
  using BasisPtr = typename BasisFamilyCG::BasisPtr;
  
  bool success = true;
  int numPoints = 8;
  
  BasisPtr line_basis;
  
  if (defineVertexFunctions)
  {
    line_basis = getLineBasis<BasisFamilyCG>(fs, N);
  }
  else
  {
    line_basis = getLineBasis<BasisFamilyDG>(fs, N);
  }
  
  auto lineTopo    = line_basis->getBaseCellTopology();
  auto inputPoints = getInputPointsView<double,DeviceType>(lineTopo, numPoints);
  ViewType<double,DeviceType> outputValues("compatible vertices H^1 value output",N+1,numPoints);
  ViewType<double,DeviceType> outputValues_dx("compatible vertices H^1 derivative output",N+1,numPoints);
  
  line_basis->getValues(outputValues,    inputPoints, OPERATOR_VALUE);
  if (fs == FUNCTION_SPACE_HGRAD)
  {
    line_basis->getValues(outputValues_dx, inputPoints, OPERATOR_GRAD);
  }
  
  for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
  {
    const double x_intrepid2 = inputPoints(pointOrdinal,0);
    const double xi = toZeroOne(x_intrepid2);
    const int polyOrder = N;
    const int size[2] = {N,VALUE_LENGTH_1D};
    int dofCount = -1;
    double values[VALUE_LENGTH_1D];
    double values_dx[VALUE_LENGTH_1D];
    
    if (fs == FUNCTION_SPACE_HGRAD)
    {
      shape1dhseg_(xi, polyOrder, size, dofCount, values, values_dx);
    }
    else if (fs == FUNCTION_SPACE_HVOL)
    {
      const int polyOrder_H1 = polyOrder+1; // ESEAS uses a the same polynomial order across the exact sequence, so that L^2 functions are of order one less than the specified order.  Intrepid2 always uses the maximum polynomial degree of the basis functions to define its polynomial order, so we add 1 here.
      shape1dqseg_(xi, polyOrder_H1, size, dofCount, values);
    }
    
    if (!defineVertexFunctions && (fs == FUNCTION_SPACE_HGRAD))
    {
      // overwrite the field ordinal 0 case
      values[0] = 1.0;
      values_dx[0] = 0.0;
    }
    
    for (int fieldOrdinal=0; fieldOrdinal<VALUE_LENGTH_1D; fieldOrdinal++)
    {
      auto actual_value   = outputValues(fieldOrdinal,pointOrdinal);
      auto expected_value = values[fieldOrdinal];
      
      bool values_match = valuesMatchToTolerance(actual_value, expected_value);
      
      if (!values_match)
      {
        std::cout << "values for x = "  << x_intrepid2 << " differ for field ordinal " << fieldOrdinal;
        std::cout << ": expected " << expected_value << "; actual " << actual_value;
        std::cout << " (diff: " << expected_value-actual_value << ")" << std::endl;
        success = false;
      }
      
      if (fs == FUNCTION_SPACE_HGRAD)
      {
        auto actual_value_dx   = outputValues_dx(fieldOrdinal,pointOrdinal);
        auto expected_value_dx = fromZeroOne_dx(values_dx[fieldOrdinal]);

        bool dxes_match = valuesMatchToTolerance(actual_value_dx, expected_value_dx);

        if (!dxes_match)
        {
          std::cout << "derivatives for x = "  << x_intrepid2 << " differ for field ordinal " << fieldOrdinal;
          std::cout << ": expected " << expected_value_dx << "; actual " << actual_value_dx;
          std::cout << " (diff: " << expected_value_dx-actual_value_dx << ")" << std::endl;
          success = false;
        }
      }
    }
  }
  return success;
}

bool testTensorTopologyMapLineLine()
{
  bool success = true;
  std::vector< std::pair<unsigned,unsigned> > expectedNodePairs;
  expectedNodePairs.push_back({0,0}); // quad vertex 0
  expectedNodePairs.push_back({1,0}); // quad vertex 1
  expectedNodePairs.push_back({1,1}); // quad vertex 2
  expectedNodePairs.push_back({0,1}); // quad vertex 3
  
  shards::CellTopology line = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
  shards::CellTopology quad = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
  
  std::vector< std::pair<unsigned,unsigned> > nodePairs = Intrepid2::TensorTopologyMap::defaultNodePairs(line, line, quad);
  if (nodePairs.size() != expectedNodePairs.size())
  {
    success = false;
    std::cout << "FAILURE: Expected node pairs size does not match default node pairs size from TensorTopologyMap.\n";
  }
  else
  {
    for (int i=0; i<int(expectedNodePairs.size()); i++)
    {
      if (expectedNodePairs[i] != nodePairs[i])
      {
        success = false;
        std::cout << "FAILURE: expected node pair " << i << " to be (" << expectedNodePairs[i].first;
        std::cout << "," << expectedNodePairs[i].second << "), but default node pair is (";
        std::cout << nodePairs[i].first << "," << nodePairs[i].second << ".\n";
      }
    }
  }
  
  Intrepid2::TensorTopologyMap topoMap(line,line);
  // test vertices -- these should just match the node pairs specified above
  unsigned expectedNodeOrdinal = 0;
  for (auto nodePair : expectedNodePairs)
  {
    unsigned vertexDim = 0;
    unsigned actualNodeOrdinal = topoMap.getCompositeSubcellOrdinal(vertexDim, nodePair.first, vertexDim, nodePair.second);
    if (actualNodeOrdinal != expectedNodeOrdinal)
    {
      std::cout << "FAILURE: Expected node ordinal of " << expectedNodeOrdinal << " but got " << actualNodeOrdinal << std::endl;
      success = false;
    }
    expectedNodeOrdinal++;
  }
  using Subcell = std::pair<unsigned,unsigned>;
  using SubcellPair = std::pair<Subcell,Subcell>;
  std::map<SubcellPair,unsigned> expectedCompositeSubcells;
  Subcell subcell1, subcell2;
  // edges
  // edge 0 times vertex 0 is edge 0
  subcell1 = {1,0};
  subcell2 = {0,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 0;
  // edge 0 times vertex 1 is edge 2
  subcell1 = {1,0};
  subcell2 = {0,1};
  expectedCompositeSubcells[{subcell1,subcell2}] = 2;
  // vertex 0 times edge 0 is edge 3
  subcell1 = {0,0};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 3;
  // vertex 1 times edge 0 is edge 1
  subcell1 = {0,1};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 1;
  
  // face
  // edge 0 times edge 0 is face 0
  subcell1 = {1,0};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 0;
  
  for (auto entry : expectedCompositeSubcells)
  {
    auto subcell1 = entry.first.first;
    auto subcell2 = entry.first.second;
    auto expectedCompositeOrdinal = entry.second;
    unsigned actualOrdinal = topoMap.getCompositeSubcellOrdinal(subcell1.first, subcell1.second,
                                                                subcell2.first, subcell2.second);
    if (actualOrdinal != expectedCompositeOrdinal)
    {
      std::cout << "FAILURE: Expected ordinal of " << expectedCompositeOrdinal << " but got " << actualOrdinal << std::endl;
      success = false;
    }
  }
  
  return success;
}

bool testTensorTopologyMapLineQuad()
{
  bool success = true;
  std::vector< std::pair<unsigned,unsigned> > expectedNodePairs;
  expectedNodePairs.push_back({0,0}); // hex vertex 0
  expectedNodePairs.push_back({1,0}); // hex vertex 1
  expectedNodePairs.push_back({1,1}); // hex vertex 2
  expectedNodePairs.push_back({0,1}); // hex vertex 3
  expectedNodePairs.push_back({0,3}); // hex vertex 4
  expectedNodePairs.push_back({1,3}); // hex vertex 5
  expectedNodePairs.push_back({1,2}); // hex vertex 6
  expectedNodePairs.push_back({0,2}); // hex vertex 7
  
  shards::CellTopology line = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
  shards::CellTopology quad = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
  shards::CellTopology hex  = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
  
  std::vector< std::pair<unsigned,unsigned> > nodePairs = Intrepid2::TensorTopologyMap::defaultNodePairs(line, quad, hex);
  if (nodePairs.size() != expectedNodePairs.size())
  {
    success = false;
    std::cout << "FAILURE: Expected node pairs size does not match default node pairs size from TensorTopologyMap.\n";
  }
  else
  {
    for (int i=0; i<int(expectedNodePairs.size()); i++)
    {
      if (expectedNodePairs[i] != nodePairs[i])
      {
        success = false;
        std::cout << "FAILURE: expected node pair " << i << " to be (" << expectedNodePairs[i].first;
        std::cout << "," << expectedNodePairs[i].second << "), but default node pair is (";
        std::cout << nodePairs[i].first << "," << nodePairs[i].second << ").\n";
      }
    }
  }
  
  Intrepid2::TensorTopologyMap topoMap(line,quad);
  // test vertices -- these should just match the node pairs specified above
  unsigned expectedNodeOrdinal = 0;
  for (auto nodePair : expectedNodePairs)
  {
    unsigned vertexDim = 0;
    unsigned actualNodeOrdinal = topoMap.getCompositeSubcellOrdinal(vertexDim, nodePair.first, vertexDim, nodePair.second);
    if (actualNodeOrdinal != expectedNodeOrdinal)
    {
      std::cout << "FAILURE: Expected node ordinal of " << expectedNodeOrdinal << " but got " << actualNodeOrdinal << std::endl;
      success = false;
    }
    expectedNodeOrdinal++;
  }
  // TODO: fix the expected values / comments below
  using Subcell = std::pair<unsigned,unsigned>;
  using SubcellPair = std::pair<Subcell,Subcell>;
  std::map<SubcellPair,unsigned> expectedCompositeSubcells;
  Subcell subcell1, subcell2;
  // edges
  // vertex 0 times edge 0 is edge 3
  subcell1 = {0,0};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 3;
  // vertex 0 times edge 1 is edge 11
  subcell1 = {0,0};
  subcell2 = {1,1};
  expectedCompositeSubcells[{subcell1,subcell2}] = 11;
  // vertex 0 times edge 2 is edge 7
  subcell1 = {0,0};
  subcell2 = {1,2};
  expectedCompositeSubcells[{subcell1,subcell2}] = 7;
  // vertex 0 times edge 3 is edge 8
  subcell1 = {0,0};
  subcell2 = {1,3};
  expectedCompositeSubcells[{subcell1,subcell2}] = 8;
  // vertex 1 times edge 0 is edge 1
  subcell1 = {0,1};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 1;
  // vertex 1 times edge 1 is edge 10
  subcell1 = {0,1};
  subcell2 = {1,1};
  expectedCompositeSubcells[{subcell1,subcell2}] = 10;
  // vertex 1 times edge 2 is edge 5
  subcell1 = {0,1};
  subcell2 = {1,2};
  expectedCompositeSubcells[{subcell1,subcell2}] = 5;
  // vertex 1 times edge 3 is edge 8
  subcell1 = {0,1};
  subcell2 = {1,3};
  expectedCompositeSubcells[{subcell1,subcell2}] = 9;
  
  // hex faces go in this order: bottom, right, top, left, front, back
  // vertex 0 times face 0 is left -- face 3
  subcell1 = {0,0};
  subcell2 = {2,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 3;
  // vertex 1 times face 0 is right -- face 1
  subcell1 = {0,1};
  subcell2 = {2,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 1;
  
  // edge 0 times edge 0 is front -- face 4
  subcell1 = {1,0};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 4;
  // edge 0 times edge 1 is top -- face 2
  subcell1 = {1,0};
  subcell2 = {1,1};
  expectedCompositeSubcells[{subcell1,subcell2}] = 2;
  // edge 0 times edge 2 is back -- face 5
  subcell1 = {1,0};
  subcell2 = {1,2};
  expectedCompositeSubcells[{subcell1,subcell2}] = 5;
  // edge 0 times edge 3 is bottom -- face 0
  subcell1 = {1,0};
  subcell2 = {1,3};
  expectedCompositeSubcells[{subcell1,subcell2}] = 0;
  
  // volume
  // edge 0 times face 0 is interior -- subcell 0
  subcell1 = {1,0};
  subcell2 = {2,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 0;
  
  for (auto entry : expectedCompositeSubcells)
  {
    auto subcell1 = entry.first.first;
    auto subcell2 = entry.first.second;
    auto expectedCompositeOrdinal = entry.second;
    unsigned actualOrdinal = topoMap.getCompositeSubcellOrdinal(subcell1.first, subcell1.second,
                                                                subcell2.first, subcell2.second);
    if (actualOrdinal != expectedCompositeOrdinal)
    {
      std::cout << "FAILURE: Expected ordinal of " << expectedCompositeOrdinal << " but got " << actualOrdinal << std::endl;
      success = false;
    }
  }
  
  return success;
}

bool testTensorTopologyMapQuadLine()
{
  bool success = true;
  std::vector< std::pair<unsigned,unsigned> > expectedNodePairs;
  expectedNodePairs.push_back({0,0}); // hex vertex 0
  expectedNodePairs.push_back({1,0}); // hex vertex 1
  expectedNodePairs.push_back({2,0}); // hex vertex 2
  expectedNodePairs.push_back({3,0}); // hex vertex 3
  expectedNodePairs.push_back({0,1}); // hex vertex 4
  expectedNodePairs.push_back({1,1}); // hex vertex 5
  expectedNodePairs.push_back({2,1}); // hex vertex 6
  expectedNodePairs.push_back({3,1}); // hex vertex 7
  
  shards::CellTopology line = shards::CellTopology(shards::getCellTopologyData<shards::Line<2> >() );
  shards::CellTopology quad = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<4> >() );
  shards::CellTopology hex  = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<8> >() );
  
  std::vector< std::pair<unsigned,unsigned> > nodePairs = Intrepid2::TensorTopologyMap::defaultNodePairs(quad, line, hex);
  if (nodePairs.size() != expectedNodePairs.size())
  {
    success = false;
    std::cout << "FAILURE: Expected node pairs size does not match default node pairs size from TensorTopologyMap.\n";
  }
  else
  {
    for (int i=0; i<int(expectedNodePairs.size()); i++)
    {
      if (expectedNodePairs[i] != nodePairs[i])
      {
        success = false;
        std::cout << "FAILURE: expected node pair " << i << " to be (" << expectedNodePairs[i].first;
        std::cout << "," << expectedNodePairs[i].second << "), but default node pair is (";
        std::cout << nodePairs[i].first << "," << nodePairs[i].second << ").\n";
      }
    }
  }
  
  Intrepid2::TensorTopologyMap topoMap(quad,line);
  // test vertices -- these should just match the node pairs specified above
  unsigned expectedNodeOrdinal = 0;
  for (auto nodePair : expectedNodePairs)
  {
    unsigned vertexDim = 0;
    unsigned actualNodeOrdinal = topoMap.getCompositeSubcellOrdinal(vertexDim, nodePair.first, vertexDim, nodePair.second);
    if (actualNodeOrdinal != expectedNodeOrdinal)
    {
      std::cout << "FAILURE: Expected node ordinal of " << expectedNodeOrdinal << " but got " << actualNodeOrdinal << std::endl;
      success = false;
    }
    expectedNodeOrdinal++;
  }
  using Subcell = std::pair<unsigned,unsigned>;
  using SubcellPair = std::pair<Subcell,Subcell>;
  std::map<SubcellPair,unsigned> expectedCompositeSubcells;
  Subcell subcell1, subcell2;
  // edges
  // edge 0 times vertex 0 is edge 0
  subcell1 = {1,0};
  subcell2 = {0,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 0;
  // edge 1 times vertex 0 is edge 1
  subcell1 = {1,1};
  subcell2 = {0,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 1;
  // edge 2 times vertex 0 is edge 2
  subcell1 = {1,2};
  subcell2 = {0,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 2;
  // edge 3 times vertex 0 is edge 3
  subcell1 = {1,3};
  subcell2 = {0,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 3;
  
  // edge 0 times vertex 1 is edge 4
  subcell1 = {1,0};
  subcell2 = {0,1};
  expectedCompositeSubcells[{subcell1,subcell2}] = 4;
  // edge 1 times vertex 1 is edge 5
  subcell1 = {1,1};
  subcell2 = {0,1};
  expectedCompositeSubcells[{subcell1,subcell2}] = 5;
  // edge 2 times vertex 1 is edge 6
  subcell1 = {1,2};
  subcell2 = {0,1};
  expectedCompositeSubcells[{subcell1,subcell2}] = 6;
  // edge 3 times vertex 1 is edge 7
  subcell1 = {1,3};
  subcell2 = {0,1};
  expectedCompositeSubcells[{subcell1,subcell2}] = 7;
  
  // vertex 0 times edge 0 is edge 4
  subcell1 = {0,0};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 8;
  // vertex 1 times edge 0 is edge 5
  subcell1 = {0,1};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 9;
  // vertex 2 times edge 0 is edge 6
  subcell1 = {0,2};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 10;
  // vertex 3 times edge 0 is edge 7
  subcell1 = {0,3};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 11;

  // hex faces go in this order: bottom, right, top, left, front, back
  // face 0 times vertex 0 is front -- face 4
  subcell1 = {2,0};
  subcell2 = {0,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 4;
  // face 0 times vertex 1 is back -- face 5
  subcell1 = {2,0};
  subcell2 = {0,1};
  expectedCompositeSubcells[{subcell1,subcell2}] = 5;

  // edge 0 times edge 0 is bottom -- face 0
  subcell1 = {1,0};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 0;
  // edge 1 times edge 0 is right -- face 1
  subcell1 = {1,1};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 1;
  // edge 2 times edge 0 is top -- face 2
  subcell1 = {1,2};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 2;
  // edge 3 times edge 0 is top -- face 3
  subcell1 = {1,3};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 3;
  
  // volume
  // face 0 times edge 0 is interior -- subcell 0
  subcell1 = {2,0};
  subcell2 = {1,0};
  expectedCompositeSubcells[{subcell1,subcell2}] = 0;
  
  for (auto entry : expectedCompositeSubcells)
  {
    auto subcell1 = entry.first.first;
    auto subcell2 = entry.first.second;
    auto expectedCompositeOrdinal = entry.second;
    unsigned actualOrdinal = topoMap.getCompositeSubcellOrdinal(subcell1.first, subcell1.second,
                                                                subcell2.first, subcell2.second);
    if (actualOrdinal != expectedCompositeOrdinal)
    {
      std::cout << "FAILURE: Expected ordinal of " << expectedCompositeOrdinal << " but got " << actualOrdinal << std::endl;
      success = false;
    }
  }
  
  return success;
}

bool testTensorTopologyMapCompositeAssignment()
{
  bool success = true;
  // test that the default assignment of shards::CellTopology for each supported pairing does what we expect
  shards::CellTopology node  = shards::CellTopology(shards::getCellTopologyData<shards::Node >() );
  shards::CellTopology line  = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
  shards::CellTopology quad  = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >() );
  shards::CellTopology hex   = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >() );
  shards::CellTopology tri   = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<> >() );
  shards::CellTopology wedge = shards::CellTopology(shards::getCellTopologyData<shards::Wedge<> >() );

  std::vector<shards::CellTopology> allTopos = {node,line,quad,hex,tri,wedge};
  // test that node x topo, topo x node, return the topo
  for (auto topo : allTopos)
  {
    auto composite = Intrepid2::TensorTopologyMap::compositeCellTopology(node, topo);
    if (composite.getBaseKey() != topo.getBaseKey())
    {
      success = false;
      std::cout << "FAILURE: Node x " << topo.getBaseName() << " is not " << topo.getBaseName() << std::endl;
    }
    composite = Intrepid2::TensorTopologyMap::compositeCellTopology(topo, node);
    if (composite.getBaseKey() != topo.getBaseKey())
    {
      success = false;
      std::cout << "FAILURE: " << topo.getBaseName() << " x Node is not " << topo.getBaseName() << std::endl;
    }
  }
  // line x blah
  std::vector<shards::CellTopology> secondTopos    = {line,quad};
  std::vector<shards::CellTopology> compositeTopos = {quad,hex};
  for (int i=0; i<int(secondTopos.size()); i++)
  {
    auto second = secondTopos[i];
    auto expectedComposite = compositeTopos[i];
    auto composite = Intrepid2::TensorTopologyMap::compositeCellTopology(line, second);
    
    if (composite.getBaseKey() != expectedComposite.getBaseKey())
    {
      success = false;
      std::cout << "FAILURE: Expected Line x " << second.getBaseName() << " to be " << expectedComposite.getBaseName();
      std::cout << "; got " << composite.getBaseName() << " instead.\n";
    }
  }
  // blah x line
  std::vector<shards::CellTopology> firstTopos = {line,tri,quad};
  compositeTopos = {quad,wedge,hex};
  for (int i=0; i<int(firstTopos.size()); i++)
  {
    auto first = firstTopos[i];
    auto expectedComposite = compositeTopos[i];
    auto composite = Intrepid2::TensorTopologyMap::compositeCellTopology(first, line);
    
    if (composite.getBaseKey() != expectedComposite.getBaseKey())
    {
      success = false;
      std::cout << "FAILURE: Expected " << first.getBaseName() << " x Line to be " << expectedComposite.getBaseName();
      std::cout << "; got " << composite.getBaseName() << " instead.\n";
    }
  }
  
  return success;
}

bool testSubBasis(shards::CellTopology cellTopo, Intrepid2::EFunctionSpace fs,
                  int polyOrder_x, int polyOrder_y=-1, int polyOrder_z = -1, bool defineVertexFunctions = true)
{
  bool success = true;
  
  using ExecSpace = Kokkos::DefaultExecutionSpace;
  using Scalar = double;
  using namespace Intrepid2;
  using Basis = Basis<DeviceType,Scalar,double>;
  
  auto vectorsMatch = [](std::vector<int> lhs,std::vector<int> rhs)
  {
    if (lhs.size() != rhs.size()) return false;
    for (int i=0; i<int(lhs.size()); i++)
    {
      if (lhs[i] != rhs[i]) return false;
    }
    return true;
  };
  
  Teuchos::RCP<Basis> basis;
  if (defineVertexFunctions)
  {
    basis = getHierarchicalBasis<true>(cellTopo, fs, polyOrder_x, polyOrder_y, polyOrder_z);
  }
  else
  {
    basis = getHierarchicalBasis<false>(cellTopo, fs, polyOrder_x, polyOrder_y, polyOrder_z);
  }
  int spaceDim = cellTopo.getDimension();
  int minDegree = (fs == Intrepid2::FUNCTION_SPACE_HVOL) ? 0 : 1;
  
  auto subBasisDegreeTestCases = getTestCasesUpToDegree(spaceDim, minDegree, polyOrder_x, polyOrder_y, polyOrder_z);
  
  std::vector<int> degrees(spaceDim);
  degrees[0] = polyOrder_x;
  if (spaceDim > 1) degrees[1] = polyOrder_y;
  if (spaceDim > 2) degrees[2] = polyOrder_z;
//  // TODO: consider what to do when non-hypercubes are being tested
  
  int numPoints_1D = 5;
  auto inputPoints = getInputPointsView<double,DeviceType>(cellTopo, numPoints_1D);
  int numPoints = inputPoints.extent_int(0);
  
  auto op = Intrepid2::OPERATOR_VALUE;
  auto outputValues = getOutputView<double,DeviceType>(fs, op, basis->getCardinality(), numPoints, spaceDim);
  
  basis->getValues(outputValues, inputPoints, op);
  
  for (auto testCase : subBasisDegreeTestCases)
  {
    // pad test case with -1s to get to 3D (so we can make one call to getBasis below)
    auto paddedTestCase = testCase;
    for (int d=testCase.size(); d<3; d++)
    {
      paddedTestCase.push_back(-1);
    }
    Teuchos::RCP<Basis> subBasis;
    if (defineVertexFunctions)
    {
      subBasis = getHierarchicalBasis<true>(cellTopo, fs, paddedTestCase[0], paddedTestCase[1], paddedTestCase[2]);
    }
    else
    {
      subBasis = getHierarchicalBasis<false>(cellTopo, fs, paddedTestCase[0], paddedTestCase[1], paddedTestCase[2]);
    }
    auto allLowerOrderFieldOrdinals = basis->getFieldOrdinalsForDegree(testCase);
    // for H(curl) and H(div), allLowerOrderFieldOrdinals may have some extra entries beyond what the subBasis has
    // (this has to do with the fact that the different families have different degrees in each dimension)
    // we filter according to polynomial degree of the fields, and place the result in lowerOrderFieldOrdinals
    std::vector<int> lowerOrderFieldOrdinals;
    if ((fs == Intrepid2::FUNCTION_SPACE_HCURL) || (fs == Intrepid2::FUNCTION_SPACE_HDIV))
    {
      int allFieldOrdinalIndex = 0;
      for (int subBasisFieldOrdinal=0; subBasisFieldOrdinal < subBasis->getCardinality(); subBasisFieldOrdinal++)
      {
        while (   (allFieldOrdinalIndex<int(allLowerOrderFieldOrdinals.size()))
               && (
                   !vectorsMatch(
                      basis   ->getPolynomialDegreeOfFieldAsVector(allLowerOrderFieldOrdinals[allFieldOrdinalIndex]),
                      subBasis->getPolynomialDegreeOfFieldAsVector(subBasisFieldOrdinal) // same poly order in each dimension
                    )
                  )
              )
        {
          allFieldOrdinalIndex++;
        }
        if (allFieldOrdinalIndex < int(allLowerOrderFieldOrdinals.size()))
        {
          // then the subcell dim and subcell ordinal must match -- add this into our filtered list
          lowerOrderFieldOrdinals.push_back(allLowerOrderFieldOrdinals[allFieldOrdinalIndex]);
        }
        // if all has gone well, we've established a correspondence between subBasisFieldOrdinal and allLowerOrderFieldOrdinals[allFieldOrdinalIndex]
        // we therefore should increment allFieldOrdinalIndex
        allFieldOrdinalIndex++;
      }
    }
    else
    {
      lowerOrderFieldOrdinals = allLowerOrderFieldOrdinals;
    }
    if (int(lowerOrderFieldOrdinals.size()) != subBasis->getCardinality())
    {
      success = false;
      std::cout << "FAILURE: for test case {";
      for (int d=0; d<int(testCase.size()); d++)
      {
        std::cout << testCase[d];
        if (d<int(testCase.size()-1)) std::cout << ",";
      }
      std::cout << "}, expected fieldOrdinals for degree to have " << subBasis->getCardinality() << " entries, but had " << lowerOrderFieldOrdinals.size() << std::endl;
      std::cout.flush();
      continue; // next test case
    }
//    std::cout << "Sub-basis ordinals of degree {";
//    for (int d=0; d<testCase.size(); d++)
//    {
//      std::cout << testCase[d];
//      if (d<testCase.size()-1) std::cout << ",";
//    }
//    std::cout << "} in basis of degree {";
//    for (int d=0; d<testCase.size(); d++)
//    {
//      std::cout << degrees[d];
//      if (d<testCase.size()-1) std::cout << ",";
//    }
//    std::cout << "}: {";
//    for (int fieldOrdinal=0; fieldOrdinal<lowerOrderFieldOrdinals.size(); fieldOrdinal++)
//    {
//      std::cout << lowerOrderFieldOrdinals[fieldOrdinal];
//      if (fieldOrdinal<lowerOrderFieldOrdinals.size()-1) std::cout << ",";
//    }
//    std::cout << "}\n";
    auto subBasisOutputValues = getOutputView<double,DeviceType>(fs, op, subBasis->getCardinality(), numPoints, spaceDim);
    subBasis->getValues(subBasisOutputValues, inputPoints, op);
    
    bool vectorValued = (outputValues.rank() == 3); // F,P,D -- if scalar-valued, F,P
    
    for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
    {
      int pointPassed = true;
      // by construction, the sub-basis should have fields in the same order as the original basis
      // (the original basis just has extra fields interspersed)
      int subBasisFieldOrdinal = 0;
      for (int fieldOrdinal : lowerOrderFieldOrdinals)
      {
        if (!vectorValued)
        {
          double originalValue = outputValues(fieldOrdinal,pointOrdinal);
          double subBasisValue = subBasisOutputValues(subBasisFieldOrdinal,pointOrdinal);
          
          bool valuesMatch = valuesMatchToTolerance(originalValue, subBasisValue);
          
          if (!valuesMatch)
          {
            pointPassed = false;
            if (fs == Intrepid2::FUNCTION_SPACE_HCURL)
            {
              // scalar values are the curls
              std::cout << "curl ";
            }
            else if (fs == Intrepid2::FUNCTION_SPACE_HDIV)
            {
              // scalar values are the div values
              std::cout << "div ";
            }
            double x = inputPoints(pointOrdinal,0);
            double y = (spaceDim > 1) ? inputPoints(pointOrdinal,1) : -2.0;
            double z = (spaceDim > 2) ? inputPoints(pointOrdinal,2) : -2.0;
            
            if (spaceDim == 1)
              std::cout << "values for "  << x  << " differ for field ordinal " << fieldOrdinal;
            else if (spaceDim == 2)
              std::cout << "values for ("  << x << "," << y << ") differ for field ordinal " << fieldOrdinal;
            else
              std::cout << "values for ("  << x << "," << y << "," << z << ") differ for field ordinal " << fieldOrdinal;
            std::cout << ": expected " << subBasisValue << "; actual " << originalValue;
            std::cout << " (diff: " << subBasisValue-originalValue << ")" << std::endl;
            success = false;
          }
        }
        else // vector-valued
        {
          bool valuesMatch = true;
          for (int d=0; d<spaceDim; d++)
          {
            double originalValue = outputValues(fieldOrdinal,pointOrdinal,d);
            double subBasisValue = subBasisOutputValues(subBasisFieldOrdinal,pointOrdinal,d);
            
            if (!valuesMatchToTolerance(originalValue,subBasisValue))
            {
              valuesMatch = false;
            }
          }
          
          if (!valuesMatch)
          {
            pointPassed = false;
            double x = inputPoints(pointOrdinal,0);
            double y = (spaceDim > 1) ? inputPoints(pointOrdinal,1) : -2.0;
            double z = (spaceDim > 2) ? inputPoints(pointOrdinal,2) : -2.0;
            
            if (spaceDim == 1)
              std::cout << "values for "  << x  << " differ for field ordinal " << fieldOrdinal;
            else if (spaceDim == 2)
              std::cout << "values for ("  << x << "," << y << ") differ for field ordinal " << fieldOrdinal;
            else
              std::cout << "values for ("  << x << "," << y << "," << z << ") differ for field ordinal " << fieldOrdinal;
            std::cout << ": expected value (lower-order basis fieldOrdinal " << subBasisFieldOrdinal << "): (";
            for (int d=0; d<spaceDim; d++)
            {
              std::cout << subBasisOutputValues(subBasisFieldOrdinal,pointOrdinal,d);
              if (d<spaceDim-1) std::cout << ",";
            }
            std::cout << "); actual (larger basis fieldOrdinal " << fieldOrdinal << ") was (";
            for (int d=0; d<spaceDim; d++)
            {
              std::cout << outputValues(fieldOrdinal,pointOrdinal,d);
              if (d<spaceDim-1) std::cout << ",";
            }
            std::cout << ")" << std::endl;
            success = false;
          }
        }
        subBasisFieldOrdinal++;
      }
    }
  }
  
  return success;
}

// compare basis values between standard Intrepid2 nodal basis and our derived nodal basis family members
template<class DerivedNodalBasisFamily, class StandardNodalBasisFamily>
bool testNodalBasisMatches(shards::CellTopology &cellTopo, Intrepid2::EFunctionSpace fs, Intrepid2::EOperator op, int polyOrder)
{
  using namespace Intrepid2;
  using BasisPtr = typename DerivedNodalBasisFamily::BasisPtr;
  using Teuchos::rcp;
 
  bool success = true;
  
  BasisPtr derivedBasis, standardBasis;
  
  if (cellTopo.getKey() == shards::Line<>::key)
  {
    derivedBasis  = getLineBasis<DerivedNodalBasisFamily> (fs, polyOrder);
    standardBasis = getLineBasis<StandardNodalBasisFamily>(fs, polyOrder);
  }
  else if (cellTopo.getKey() == shards::Quadrilateral<>::key)
  {
    derivedBasis  = getQuadrilateralBasis<DerivedNodalBasisFamily> (fs, polyOrder); // derived basis supports both isotropic and anisotropic polyOrder
    standardBasis = getQuadrilateralBasis<StandardNodalBasisFamily>(fs, polyOrder); // standard basis is isotropic
  }
  else if (cellTopo.getKey() == shards::Hexahedron<>::key)
  {
    derivedBasis  = getHexahedronBasis<DerivedNodalBasisFamily> (fs, polyOrder); // derived basis supports both isotropic and anisotropic polyOrder
    standardBasis = getHexahedronBasis<StandardNodalBasisFamily>(fs, polyOrder); // isotropic
  }
  
  int standardCardinality = standardBasis->getCardinality();
  int derivedCardinality = derivedBasis->getCardinality();
  
  if (standardCardinality != derivedCardinality)
  {
    success = false;
    std::cout << "FAILURE: standard basis cardinality (" << standardCardinality;
    std::cout << ") does not match derived basis cardinality (" << derivedCardinality << ")\n";
  }
  
  int spaceDim = cellTopo.getDimension();
  // we do allow ordering to be different.  dofMapToDerived maps from standard field ordinal to the derived.
  // we use getDofCoords() to perform the mapping.
  using OutputValueType        = typename DerivedNodalBasisFamily::Basis::OutputValueType;
  using ExecutionSpace         = typename DerivedNodalBasisFamily::Basis::ExecutionSpace;
  using ViewType               = Kokkos::DynRankView<OutputValueType,ExecutionSpace>;
  using OrdinalTypeArray1DHost = typename DerivedNodalBasisFamily::Basis::OrdinalTypeArray1DHost;
  using OrdinalTypeArray1D     = typename DerivedNodalBasisFamily::Basis::OrdinalTypeArray1D;
  OrdinalTypeArray1D     dofMapToDerived     = OrdinalTypeArray1D("dofMapToDerived",standardCardinality);
  OrdinalTypeArray1DHost dofMapToDerivedHost = Kokkos::create_mirror_view(dofMapToDerived);
  
  ViewType dofCoordsStandard("dofCoordsStandard", standardCardinality, spaceDim);
  ViewType dofCoordsDerived ("dofCoordsDerived",  standardCardinality, spaceDim);
  standardBasis->getDofCoords(dofCoordsStandard);
  derivedBasis-> getDofCoords(dofCoordsDerived );
  
  Kokkos::deep_copy(dofMapToDerived, -1);
  
  Kokkos::parallel_for(standardCardinality, KOKKOS_LAMBDA (const int fieldOrdinalStandard)
  {
    // search for a fieldOrdinalDerived that matches
    bool matchFound = false;
    for (int fieldOrdinalDerived=0; fieldOrdinalDerived<standardCardinality; fieldOrdinalDerived++)
    {
      bool matches = true;
      for (int d=0; d<spaceDim; d++)
      {
        if (dofCoordsStandard(fieldOrdinalStandard,d) != dofCoordsDerived(fieldOrdinalDerived,d))
        {
          matches = false;
          break;
        }
      }
      
      if (matches)
      {
        dofMapToDerived(fieldOrdinalStandard) = fieldOrdinalDerived;
        matchFound = true;
        break;
      }
    }
    if (!matchFound)
    {
      // DEBUGGING -- note that this code is invalid on device!
      std::cout << "match not found.\n";
    }
  });
  // copy dofMapToDerived to host view
  Kokkos::deep_copy(dofMapToDerivedHost, dofMapToDerived);
  
  int numPoints_1D = 5;
  auto inputPointsView = getInputPointsView<double,DeviceType>(cellTopo, numPoints_1D);
  int numPoints = inputPointsView.extent_int(0);
  auto standardOutputView = getOutputView<double,DeviceType>(fs, op, standardCardinality, numPoints, spaceDim);
  auto derivedOutputView  = getOutputView<double,DeviceType>(fs, op, standardCardinality, numPoints, spaceDim);
  
  standardBasis->getValues(standardOutputView, inputPointsView, op);
  derivedBasis->getValues(derivedOutputView, inputPointsView, op);
  
  bool scalarValued = (standardOutputView.rank() == 2); // F,P -- if vector-valued, F,P,D or F,P,Dk
  
  for (int pointOrdinal=0; pointOrdinal<numPoints; pointOrdinal++)
  {
    int pointPassed = true;
    for (int fieldOrdinalStandard=0; fieldOrdinalStandard<derivedCardinality; fieldOrdinalStandard++)
    {
      int fieldOrdinalDerived = dofMapToDerivedHost(fieldOrdinalStandard);
      if (scalarValued)
      {
        double standardValue = standardOutputView(fieldOrdinalStandard,pointOrdinal);
        double derivedValue  =  derivedOutputView(fieldOrdinalDerived, pointOrdinal);
        
        bool valuesMatch = valuesMatchToTolerance(standardValue, derivedValue);
        
        if (!valuesMatch)
        {
          pointPassed = false;
          if (fs == Intrepid2::FUNCTION_SPACE_HCURL)
          {
            // scalar values are the curls
            std::cout << "curl ";
          }
          else if (fs == Intrepid2::FUNCTION_SPACE_HDIV)
          {
            // scalar values are the div values
            std::cout << "div ";
          }
          double x = inputPointsView(pointOrdinal,0);
          double y = (spaceDim > 1) ? inputPointsView(pointOrdinal,1) : -2.0;
          double z = (spaceDim > 2) ? inputPointsView(pointOrdinal,2) : -2.0;
          
          if (spaceDim == 1)
            std::cout << "values for "  << x  << " differ for field ordinal " << fieldOrdinalStandard;
          else if (spaceDim == 2)
            std::cout << "values for ("  << x << "," << y << ") differ for field ordinal " << fieldOrdinalStandard;
          else
            std::cout << "values for ("  << x << "," << y << "," << z << ") differ for field ordinal " << fieldOrdinalStandard;
          std::cout << ": expected " << standardValue << "; actual " << derivedValue;
          std::cout << " (diff: " << standardValue-derivedValue << ")" << std::endl;
          success = false;
        }
      }
      else // vector-valued
      {
        bool valuesMatch = true;
        int dkcard = standardOutputView.extent_int(2);
        for (int d=0; d<dkcard; d++)
        {
          double standardValue = standardOutputView(fieldOrdinalStandard,pointOrdinal,d);
          double derivedValue  =  derivedOutputView(fieldOrdinalDerived, pointOrdinal,d);
          
          if (! valuesMatchToTolerance(standardValue, derivedValue))
          {
            valuesMatch = false;
          }
        }
        
        if (!valuesMatch)
        {
          pointPassed = false;
          double x = inputPointsView(pointOrdinal,0);
          double y = (spaceDim > 1) ? inputPointsView(pointOrdinal,1) : -2.0;
          double z = (spaceDim > 2) ? inputPointsView(pointOrdinal,2) : -2.0;
          
          if (spaceDim == 1)
            std::cout << "values for "  << x  << " differ for field ordinal " << fieldOrdinalStandard;
          else if (spaceDim == 2)
            std::cout << "values for ("  << x << "," << y << ") differ for field ordinal " << fieldOrdinalStandard;
          else
            std::cout << "values for ("  << x << "," << y << "," << z << ") differ for field ordinal " << fieldOrdinalStandard;
          std::cout << ": expected (";
          for (int d=0; d<dkcard; d++)
          {
            std::cout << standardOutputView(fieldOrdinalStandard,pointOrdinal,d);
            if (d<dkcard-1) std::cout << ",";
          }
          std::cout << "); actual was (";
          for (int d=0; d<dkcard; d++)
          {
            std::cout << derivedOutputView(fieldOrdinalStandard,pointOrdinal,d);
            if (d<dkcard-1) std::cout << ",";
          }
          std::cout << ")" << std::endl;
          success = false;
        }
      }
    }
  }
  
  return success;
}

bool testHierarchicalHGRAD_LINEMatchesAnalyticValues(Intrepid2::EOperator op)
{
  bool success = true;
  using namespace Intrepid2;
  using BasisFamily = HierarchicalBasisFamily<DeviceType>;
  
  const int polyOrder = 4;
  auto hgradBasis = getLineBasis<BasisFamily>(FUNCTION_SPACE_HGRAD, polyOrder);
  
  int numPoints_1D = 5;
  shards::CellTopology lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
  auto inputPointsView = getInputPointsView<double,DeviceType>(lineTopo, numPoints_1D);
  
  auto hgradOutputView = getOutputView<double,DeviceType>(FUNCTION_SPACE_HGRAD, op, hgradBasis->getCardinality(), numPoints_1D, 1);
  
  hgradBasis->getValues(hgradOutputView, inputPointsView, op);
  
  for (int pointOrdinal=0; pointOrdinal<numPoints_1D; pointOrdinal++)
  {
    int pointPassed = true;
    double x = inputPointsView(pointOrdinal,0);
    
    std::vector<double> expectedValues(hgradBasis->getCardinality());
    switch (op)
    {
      case Intrepid2::OPERATOR_VALUE:
        expectedValues[0] = (1.0-x)/2.0;                   // left vertex function (node at -1)
        expectedValues[1] = (1.0+x)/2.0;                   // right vertex function (node at 1)
        expectedValues[2] = (x*x-1.0)/4.0;                 // L_2 : (x^2 - 1) / 4
        expectedValues[3] = (x*x*x-x)/4.0;                 // L_3 : (x^3 - x) / 4
        expectedValues[4] = (5.*x*x*x*x - 6.*x*x+1.)/16.0; // L_4 : (5x^4-6x^2+1) / 16
        break;
      case Intrepid2::OPERATOR_GRAD:
      case Intrepid2::OPERATOR_D1:
        // first derivatives of the above:
        expectedValues[0] = -1.0/2.0;              // left vertex function (node at -1)
        expectedValues[1] = 1.0/2.0;               // right vertex function (node at 1)
        expectedValues[2] = x/2.0;                 // L_2 : (x^2 - 1) / 4
        expectedValues[3] = (3.0*x*x-1.0)/4.0;     // L_3 : (x^3 - x) / 4
        expectedValues[4] = (5.*x*x*x - 3.*x)/4.0; // L_4 : (5x^4-6x^2+1) / 16
        break;
      case Intrepid2::OPERATOR_D2:
        // second derivatives:
        expectedValues[0] = 0.0;               // left vertex function (node at -1)
        expectedValues[1] = 0.0;               // right vertex function (node at 1)
        expectedValues[2] = 1.0/2.0;           // L_2 : (x^2 - 1) / 4
        expectedValues[3] = 3.0*x/2.0;         // L_3 : (x^3 - x) / 4
        expectedValues[4] = (15.*x*x - 3.)/4.; // L_4 : (5x^4-6x^2+1) / 16
        break;
      case Intrepid2::OPERATOR_D3:
        // third derivatives:
        expectedValues[0] = 0.0;               // left vertex function (node at -1)
        expectedValues[1] = 0.0;               // right vertex function (node at 1)
        expectedValues[2] = 0.0;               // L_2 : (x^2 - 1) / 4
        expectedValues[3] = 3.0/2.0;           // L_3 : (x^3 - x) / 4
        expectedValues[4] = (15.*x)/2.;        // L_4 : (5x^4-6x^2+1) / 16
        break;
      case Intrepid2::OPERATOR_D4:
        // fourth derivatives:
        expectedValues[0] = 0.0;           // left vertex function (node at -1)
        expectedValues[1] = 0.0;           // right vertex function (node at 1)
        expectedValues[2] = 0.0;           // L_2 : (x^2 - 1) / 4
        expectedValues[3] = 0.0;           // L_3 : (x^3 - x) / 4
        expectedValues[4] = 15./2.;        // L_4 : (5x^4-6x^2+1) / 16
        break;
      case Intrepid2::OPERATOR_D5:
      case Intrepid2::OPERATOR_D6:
      case Intrepid2::OPERATOR_D7:
      case Intrepid2::OPERATOR_D8:
      case Intrepid2::OPERATOR_D9:
      case Intrepid2::OPERATOR_D10:
        // fourth derivatives:
        expectedValues[0] = 0.0;           // left vertex function (node at -1)
        expectedValues[1] = 0.0;           // right vertex function (node at 1)
        expectedValues[2] = 0.0;           // L_2 : (x^2 - 1) / 4
        expectedValues[3] = 0.0;           // L_3 : (x^3 - x) / 4
        expectedValues[4] = 0.0;           // L_4 : (5x^4-6x^2+1) / 16
        break;
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator");
    }
    
    
    for (int fieldOrdinal=0; fieldOrdinal<hgradBasis->getCardinality(); fieldOrdinal++)
    {
      double actual    = hgradOutputView.access(fieldOrdinal,pointOrdinal,0);
      double expected  = expectedValues[fieldOrdinal];
      
      bool valuesMatch = valuesMatchToTolerance(actual, expected);
      
      if (!valuesMatch)
      {
        pointPassed = false;
        double x = inputPointsView(pointOrdinal,0);
        if (op == OPERATOR_VALUE) std::cout << "values";
        else
        {
          int derivativeOrder = getOperatorOrder(op);
          if (derivativeOrder == 1)
          {
            std::cout << "first ";
          }
          else if (derivativeOrder == 2)
          {
            std::cout << "second ";
          }
          else if (derivativeOrder == 3)
          {
            std::cout << "third ";
          }
          else
          {
            std::cout << derivativeOrder << "th ";
          }
          std::cout << "derivatives";
        }
        std::cout << " for "  << x  << " differ for field ordinal " << fieldOrdinal;
        std::cout << ": expected " << expected << "; actual " << actual;
        std::cout << " (diff: " << expected-actual << ")" << std::endl;
        success = false;
      }
    }
  }
  
  return success;
}

bool testHierarchicalHVOL_LINEMatchesAnalyticValues(Intrepid2::EOperator op)
{
  bool success = true;
  using namespace Intrepid2;
  using BasisFamily = HierarchicalBasisFamily<DeviceType>;
  
  const int polyOrder = 4;
  auto hgradBasis = getLineBasis<BasisFamily>(FUNCTION_SPACE_HVOL, polyOrder);
  
  int numPoints_1D = 5;
  shards::CellTopology lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
  auto inputPointsView = getInputPointsView<double,DeviceType>(lineTopo, numPoints_1D);
  
  auto hgradOutputView = getOutputView<double,DeviceType>(FUNCTION_SPACE_HGRAD, op, hgradBasis->getCardinality(), numPoints_1D, 1);
  
  hgradBasis->getValues(hgradOutputView, inputPointsView, op);
  
  for (int pointOrdinal=0; pointOrdinal<numPoints_1D; pointOrdinal++)
  {
    int pointPassed = true;
    double x = inputPointsView(pointOrdinal,0);
    
    std::vector<double> expectedValues(hgradBasis->getCardinality());
    switch (op)
    {
      case Intrepid2::OPERATOR_VALUE:
        expectedValues[0] = 1;                              // P_0 : 1
        expectedValues[1] = x;                              // P_1 : x
        expectedValues[2] = (3.*x*x-1.0)/2.0;               // P_2 : (3x^2 - 1) / 2
        expectedValues[3] = (5.*x*x*x-3.*x)/2.0;            // P_3 : (5x^3 - 3x) / 2
        expectedValues[4] = (35.*x*x*x*x - 30.*x*x+3.)/8.0; // P_4 : (35x^4-30x^2+3) / 8
        break;
      case Intrepid2::OPERATOR_GRAD:
      case Intrepid2::OPERATOR_D1:
        // first derivatives of the above:
        expectedValues[0] = 0.0;                      // P_0 : 1
        expectedValues[1] = 1.0;                      // P_1 : x
        expectedValues[2] = 3.*x;                     // P_2 : (3x^2 - 1) / 2
        expectedValues[3] = (15.0*x*x-3.0)/2.0;       // P_3 : (5x^3 - 3x) / 2
        expectedValues[4] = (35.*x*x*x - 15.*x)/2.0; // P_4 : (35x^4-30x^2+3) / 8
        break;
      case Intrepid2::OPERATOR_D2:
        // second derivatives:
        expectedValues[0] = 0.0;                 // P_0 : 1
        expectedValues[1] = 0.0;                 // P_1 : x
        expectedValues[2] = 3.0;                 // P_2 : (3x^2 - 1) / 2
        expectedValues[3] = 15.0*x;              // P_3 : (5x^3 - 3x) / 2
        expectedValues[4] = (105.*x*x - 15.)/2.; // P_4 : (35x^4-30x^2+3) / 8
        break;
      case Intrepid2::OPERATOR_D3:
        // third derivatives:
        expectedValues[0] = 0.0;           // P_0 : 1
        expectedValues[1] = 0.0;           // P_1 : x
        expectedValues[2] = 0.0;           // P_2 : (3x^2 - 1) / 2
        expectedValues[3] = 15.0;          // P_3 : (5x^3 - 3x) / 2
        expectedValues[4] = 105.*x;        // P_4 : (35x^4-30x^2+3) / 8
        break;
      case Intrepid2::OPERATOR_D4:
        // fourth derivatives:
        expectedValues[0] = 0.0;           // P_0 : 1
        expectedValues[1] = 0.0;           // P_1 : x
        expectedValues[2] = 0.0;           // P_2 : (3x^2 - 1) / 2
        expectedValues[3] = 0.0;           // P_3 : (5x^3 - 3x) / 2
        expectedValues[4] = 105.;          // P_4 : (35x^4-30x^2+3) / 8
        break;
      case Intrepid2::OPERATOR_D5:
      case Intrepid2::OPERATOR_D6:
      case Intrepid2::OPERATOR_D7:
      case Intrepid2::OPERATOR_D8:
      case Intrepid2::OPERATOR_D9:
      case Intrepid2::OPERATOR_D10:
        // fourth derivatives:
        expectedValues[0] = 0.0;           // P_0 : 1
        expectedValues[1] = 0.0;           // P_1 : x
        expectedValues[2] = 0.0;           // P_2 : (3x^2 - 1) / 2
        expectedValues[3] = 0.0;           // P_3 : (5x^3 - 3x) / 2
        expectedValues[4] = 0.0;           // P_4 : (35x^4-30x^2+3) / 8
        break;
      default:
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported operator");
    }
    
    for (int fieldOrdinal=0; fieldOrdinal<hgradBasis->getCardinality(); fieldOrdinal++)
    {
      double actual    = hgradOutputView.access(fieldOrdinal,pointOrdinal,0);
      double expected  = expectedValues[fieldOrdinal];
      
      bool valuesMatch = valuesMatchToTolerance(actual, expected);
      
      if (!valuesMatch)
      {
        pointPassed = false;
        double x = inputPointsView(pointOrdinal,0);
        if (op == OPERATOR_VALUE) std::cout << "values";
        else
        {
          int derivativeOrder = getOperatorOrder(op);
          if (derivativeOrder == 1)
          {
            std::cout << "first ";
          }
          else if (derivativeOrder == 2)
          {
            std::cout << "second ";
          }
          else if (derivativeOrder == 3)
          {
            std::cout << "third ";
          }
          else
          {
            std::cout << derivativeOrder << "th ";
          }
          std::cout << "derivatives";
        }
        std::cout << " for "  << x  << " differ for field ordinal " << fieldOrdinal;
        std::cout << ": expected " << expected << "; actual " << actual;
        std::cout << " (diff: " << expected-actual << ")" << std::endl;
        success = false;
      }
    }
  }
  
  return success;
}

// compare derivatives of order derivativeOrder in H(grad) with derivatives of order (derivativeOrder-1) in H(vol)
template<class LineBasisFamily>
bool testDerivativesMatch(int polyOrder, int derivativeOrder)
{
  using namespace Intrepid2;
  using Teuchos::rcp;
  
  bool success = true;
  
  auto hgradBasis = getLineBasis<LineBasisFamily>(FUNCTION_SPACE_HGRAD, polyOrder);
  auto hvolBasis  = getLineBasis<LineBasisFamily>(FUNCTION_SPACE_HVOL, polyOrder);
  
  auto opGrad = Intrepid2::EOperator(OPERATOR_D1 - 1 + derivativeOrder);
  auto opVol  = (derivativeOrder > 1) ? Intrepid2::EOperator(OPERATOR_D1 - 2 + derivativeOrder) : OPERATOR_VALUE;
  
  int numPoints_1D = 5;
  shards::CellTopology lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
  auto inputPointsView = getInputPointsView<double,DeviceType>(lineTopo, numPoints_1D);
  const int spaceDim = 1;
  
  auto hgradOutputView = getOutputView<double,DeviceType>(FUNCTION_SPACE_HGRAD, opGrad, hgradBasis->getCardinality(), numPoints_1D, spaceDim);
  auto hvolOutputView  = getOutputView<double,DeviceType>(FUNCTION_SPACE_HVOL,   opVol, hvolBasis->getCardinality(),  numPoints_1D, spaceDim);
  
  hgradBasis->getValues(hgradOutputView, inputPointsView, opGrad);
  hvolBasis->getValues(  hvolOutputView, inputPointsView, opVol);
  
  for (int pointOrdinal=0; pointOrdinal<numPoints_1D; pointOrdinal++)
  {
    int pointPassed = true;
    for (int fieldOrdinal=2; fieldOrdinal<hgradBasis->getCardinality(); fieldOrdinal++)
    {
      // relationship between our H^1 basis and our L^2 (see the "Analytic" tests above for explicit polynomial expressions)
      // - nth derivative of H^1 fieldOrdinal i is one half of the (n-1)th derivative of L^2 fieldOrdinal i-1 (for i>=2)
      double hgradValue = hgradOutputView(fieldOrdinal,  pointOrdinal,0);
      double hvolValue  =  hvolOutputView(fieldOrdinal-1,pointOrdinal,0);
      
      bool valuesMatch = valuesMatchToTolerance(hgradValue, 0.5 * hvolValue);
      
      if (!valuesMatch)
      {
        pointPassed = false;
        double x = inputPointsView(pointOrdinal,0);
        
        std::string gradDerivative, volDerivative;
        if (derivativeOrder == 1)
        {
          gradDerivative = "st";
          volDerivative  = "th";
        }
        else if (derivativeOrder == 2)
        {
          gradDerivative = "nd";
          volDerivative  = "st";
        }
        else if (derivativeOrder == 3)
        {
          gradDerivative = "rd";
          volDerivative  = "nd";
        }
        else if (derivativeOrder == 4)
        {
          gradDerivative = "th";
          volDerivative  = "rd";
        }
        else
        {
          gradDerivative = "th";
          volDerivative  = "th";
        }
        
        std::cout << "values for "  << x  << " differ for field ordinal " << fieldOrdinal;
        std::cout << ": h(grad) " << derivativeOrder << gradDerivative << " derivative is " << hgradValue << "; h(vol) " << derivativeOrder-1 << volDerivative << " derivative is " << hvolValue;
        std::cout << " (diff: " << hgradValue-hvolValue << ")" << std::endl;
        success = false;
      }
      
    }
  }
  
  return success;
}

template<class DerivedNodalBasisFamily, class StandardNodalBasisFamily>
void runNodalBasisComparisonTests(int &successCount, int &testCount,
                                  shards::CellTopology &cellTopo, std::vector<Intrepid2::EFunctionSpace> functionSpaces, std::vector<Intrepid2::EOperator> ops)
{
  using namespace Intrepid2;
  const int polyOrderForNodalBasisTest = 3;
  
  for (auto fs : functionSpaces)
  {
    for (auto op : ops)
    {
      if (testNodalBasisMatches<DerivedNodalBasisFamily,StandardNodalBasisFamily>(cellTopo, fs, op, polyOrderForNodalBasisTest)) successCount++;
      testCount++;
    }
  }
}

int main(int argc, char *argv[])
{
  {
    // Note that the dtor for GlobalMPISession will call Kokkos::finalize_all() but does not call Kokkos::initialize()...
    Teuchos::GlobalMPISession mpiSession(&argc, &argv);
    Kokkos::initialize(argc,argv);
    Teuchos::UnitTestRepository::setGloballyReduceTestResult(true);
    
    // Teuchos tests:
    int result = Teuchos::UnitTestRepository::runUnitTestsFromMain(argc, argv);
    
    {
      // home-grown tests -- these run after the Teuchos tests (and probably should be converted to Teuchos tests):
      int testCount = 0;
      int successCount = 0;
            
      std::vector<Intrepid2::EFunctionSpace> functionSpaces_2D = {FUNCTION_SPACE_HGRAD,FUNCTION_SPACE_HCURL,FUNCTION_SPACE_HDIV,FUNCTION_SPACE_HVOL};
      std::vector<Intrepid2::EFunctionSpace> functionSpaces_3D = {FUNCTION_SPACE_HGRAD,FUNCTION_SPACE_HCURL,FUNCTION_SPACE_HDIV,FUNCTION_SPACE_HVOL};
      
      shards::CellTopology lineTopo = shards::CellTopology(shards::getCellTopologyData<shards::Line<> >() );
      shards::CellTopology quadTopo = shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >() );
      shards::CellTopology triTopo = shards::CellTopology(shards::getCellTopologyData<shards::Triangle<> >() );
      shards::CellTopology hexTopo = shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >() );
      shards::CellTopology tetTopo = shards::CellTopology(shards::getCellTopologyData<shards::Tetrahedron<> >() );
      
      std::cout << "Testing TensorTopologyMap…\n";
      if (testTensorTopologyMapLineLine()) successCount++;
      testCount++;
      if (testTensorTopologyMapLineQuad()) successCount++;
      testCount++;
      if (testTensorTopologyMapQuadLine()) successCount++;
      testCount++;
      if (testTensorTopologyMapCompositeAssignment()) successCount++;
      testCount++;
      
      std::vector<shards::CellTopology> cellTopos = {quadTopo, hexTopo};
      
      std::vector<Intrepid2::EFunctionSpace> functionSpaces_1D = {FUNCTION_SPACE_HGRAD,FUNCTION_SPACE_HVOL};
      std::vector<std::vector<Intrepid2::EFunctionSpace> > functionSpacesForDimension = {functionSpaces_1D,functionSpaces_2D,functionSpaces_3D};
      
      std::cout << "Testing sub-basis inclusion…\n";
      for (auto cellTopo : cellTopos)
      {
        int spaceDim = cellTopo.getDimension();
        auto functionSpaces = functionSpacesForDimension[spaceDim-1];
        
        int maxDegree = 4;
        for (auto fs : functionSpaces)
        {
          std::vector<bool> continuousBasisValues;
          if (fs != FUNCTION_SPACE_HVOL)
          {
            continuousBasisValues = {true,false};
          }
          else
          {
            continuousBasisValues = {true}; // false case not supported by the dof tag stuff that testSubBasis() does
          }
          for (auto continuousBasis : continuousBasisValues) // corresponds to "defineVertexFunctions" in line basis definitions
          {
            for (int degree=1; degree<=maxDegree; degree++)
            {
              if (testSubBasis(cellTopo, fs, degree,degree,degree,continuousBasis)) successCount++;
              testCount++;
            }
          }
        }
      }
      
      std::vector<EOperator> operators_dk = {OPERATOR_D1,OPERATOR_D2,OPERATOR_D3,OPERATOR_D4,OPERATOR_D5,OPERATOR_D6,OPERATOR_D7,OPERATOR_D8,OPERATOR_D9,OPERATOR_D10};
      
      std::cout << "Running tests comparing H(grad) basis on line against analytic expected functions.\n";
      if (testHierarchicalHGRAD_LINEMatchesAnalyticValues(OPERATOR_VALUE)) successCount++;
      testCount++;
      
      for (auto op : operators_dk)
      {
        if (testHierarchicalHGRAD_LINEMatchesAnalyticValues(op)) successCount++;
        testCount++;
      }
      
      std::cout << "Running tests comparing H(vol) basis on line against analytic expected functions.\n";
      if (testHierarchicalHVOL_LINEMatchesAnalyticValues(OPERATOR_VALUE)) successCount++;
      testCount++;
      
      for (auto op : operators_dk)
      {
        if (testHierarchicalHVOL_LINEMatchesAnalyticValues(op)) successCount++;
        testCount++;
      }
      
      const int polyOrderForTests = 10;
      std::cout << "Running tests comparing high-order derivatives of H(grad) and H(vol) bases\n";
      using HierarchicalBasisFamily = HierarchicalBasisFamily<DeviceType>;
      for (int derivativeOrder=1; derivativeOrder<=5; derivativeOrder++)
      {
        if (testDerivativesMatch<HierarchicalBasisFamily>(polyOrderForTests, derivativeOrder)) successCount++;
        testCount++;
      }
      
      // the following should match in H(grad) and H(vol), but are not expected to match in H(curl) and H(div)
      // (the latter differ in that Standard uses another H(grad) instance to represent L^2, while Derived uses H(vol))
      using DerivedNodalBasisFamily  = Intrepid2::DerivedNodalBasisFamily<DeviceType>;
      using StandardNodalBasisFamily = Intrepid2::NodalBasisFamily<DeviceType>;
      
      std::cout << "Running 1D nodal basis comparison tests…\n";
      runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(successCount,testCount, lineTopo, {FUNCTION_SPACE_HVOL}, {OPERATOR_VALUE});
      runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(successCount,testCount, lineTopo, {FUNCTION_SPACE_HGRAD}, {OPERATOR_VALUE,OPERATOR_GRAD});
      
      runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(successCount,testCount, lineTopo, {FUNCTION_SPACE_HGRAD}, operators_dk);
      
      std::cout << "Running 2D nodal quadrilateral basis comparison tests…\n";
      runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(successCount,testCount, quadTopo, {FUNCTION_SPACE_HVOL}, {OPERATOR_VALUE});
      
      runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(successCount,testCount, quadTopo, {FUNCTION_SPACE_HGRAD}, {OPERATOR_VALUE,OPERATOR_GRAD});
      
      runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(successCount,testCount, quadTopo, {FUNCTION_SPACE_HGRAD}, operators_dk);
      
      std::cout << "Running 3D nodal hexahedral basis comparison tests…\n";
      runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(successCount,testCount, hexTopo, {FUNCTION_SPACE_HVOL}, {OPERATOR_VALUE});
      
      runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(successCount,testCount, hexTopo, {FUNCTION_SPACE_HGRAD}, {OPERATOR_VALUE,OPERATOR_GRAD});
      
      runNodalBasisComparisonTests<DerivedNodalBasisFamily, StandardNodalBasisFamily>(successCount,testCount, hexTopo, {FUNCTION_SPACE_HGRAD}, operators_dk);

      std::cout << "Passed " << successCount << " of " << testCount << " tests.\n";
      
      if (result != 0)
      {
        if (successCount != testCount)
        {
          result = -1;
        }
      }
    }
    
    return result;
  }
}
