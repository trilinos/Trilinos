// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   TensorTopologyMapTests.cpp
    \brief  Tests to verify TensorTopologyMap.
    \author Created by N.V. Roberts.
 */

#include "Teuchos_UnitTestHarness.hpp"

#include "Intrepid2_TensorTopologyMap.hpp"
#include "Intrepid2_Types.hpp"

#include "Intrepid2_TestUtils.hpp"

#include "Kokkos_Core.hpp"

namespace
{
  using namespace Intrepid2;
  
  void testTensorTopologyMapLineLine(Teuchos::FancyOStream &out, bool &success)
  {
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
      out << "FAILURE: Expected node pairs size does not match default node pairs size from TensorTopologyMap.\n";
    }
    else
    {
      for (unsigned i=0; i<expectedNodePairs.size(); i++)
      {
        if (expectedNodePairs[i] != nodePairs[i])
        {
          success = false;
          out << "FAILURE: expected node pair " << i << " to be (" << expectedNodePairs[i].first;
          out << "," << expectedNodePairs[i].second << "), but default node pair is (";
          out << nodePairs[i].first << "," << nodePairs[i].second << ".\n";
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
        out << "FAILURE: Expected node ordinal of " << expectedNodeOrdinal << " but got " << actualNodeOrdinal << std::endl;
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
        out << "FAILURE: Expected ordinal of " << expectedCompositeOrdinal << " but got " << actualOrdinal << std::endl;
        success = false;
      }
    }
  }
  
  void testTensorTopologyMapLineQuad(Teuchos::FancyOStream &out, bool &success)
  {
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
      out << "FAILURE: Expected node pairs size does not match default node pairs size from TensorTopologyMap.\n";
    }
    else
    {
      for (unsigned i=0; i<expectedNodePairs.size(); i++)
      {
        if (expectedNodePairs[i] != nodePairs[i])
        {
          success = false;
          out << "FAILURE: expected node pair " << i << " to be (" << expectedNodePairs[i].first;
          out << "," << expectedNodePairs[i].second << "), but default node pair is (";
          out << nodePairs[i].first << "," << nodePairs[i].second << ").\n";
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
        out << "FAILURE: Expected node ordinal of " << expectedNodeOrdinal << " but got " << actualNodeOrdinal << std::endl;
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
        out << "FAILURE: Expected ordinal of " << expectedCompositeOrdinal << " but got " << actualOrdinal << std::endl;
        success = false;
      }
    }
  }
  
  void testTensorTopologyMapQuadLine(Teuchos::FancyOStream &out, bool &success)
  {
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
      out << "FAILURE: Expected node pairs size does not match default node pairs size from TensorTopologyMap.\n";
    }
    else
    {
      for (unsigned i=0; i<expectedNodePairs.size(); i++)
      {
        if (expectedNodePairs[i] != nodePairs[i])
        {
          success = false;
          out << "FAILURE: expected node pair " << i << " to be (" << expectedNodePairs[i].first;
          out << "," << expectedNodePairs[i].second << "), but default node pair is (";
          out << nodePairs[i].first << "," << nodePairs[i].second << ").\n";
        }
      }
    }
    
    //  for (int i=0; i<6; i++)
    //  {
    //    // DEBUGGING
    //    out << "Hexahedron face " << i << " has nodes: ";
    //    out << hex.getNodeMap(2, i, 0) << ", " << hex.getNodeMap(2, i, 1) << ", " << hex.getNodeMap(2, i, 2) << ", "<< hex.getNodeMap(2, i, 3) << std::endl;
    //  }
    
    Intrepid2::TensorTopologyMap topoMap(quad,line);
    // test vertices -- these should just match the node pairs specified above
    unsigned expectedNodeOrdinal = 0;
    for (auto nodePair : expectedNodePairs)
    {
      unsigned vertexDim = 0;
      unsigned actualNodeOrdinal = topoMap.getCompositeSubcellOrdinal(vertexDim, nodePair.first, vertexDim, nodePair.second);
      if (actualNodeOrdinal != expectedNodeOrdinal)
      {
        out << "FAILURE: Expected node ordinal of " << expectedNodeOrdinal << " but got " << actualNodeOrdinal << std::endl;
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
        out << "FAILURE: Expected ordinal of " << expectedCompositeOrdinal << " but got " << actualOrdinal << std::endl;
        success = false;
      }
    }
  }
  
  void testTensorTopologyMapCompositeAssignment(Teuchos::FancyOStream &out, bool &success)
  {
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
        out << "FAILURE: Node x " << topo.getBaseName() << " is not " << topo.getBaseName() << std::endl;
      }
      composite = Intrepid2::TensorTopologyMap::compositeCellTopology(topo, node);
      if (composite.getBaseKey() != topo.getBaseKey())
      {
        success = false;
        out << "FAILURE: " << topo.getBaseName() << " x Node is not " << topo.getBaseName() << std::endl;
      }
    }
    // line x blah
    std::vector<shards::CellTopology> secondTopos    = {line,quad};
    std::vector<shards::CellTopology> compositeTopos = {quad,hex};
    for (unsigned i=0; i<secondTopos.size(); i++)
    {
      auto second = secondTopos[i];
      auto expectedComposite = compositeTopos[i];
      auto composite = Intrepid2::TensorTopologyMap::compositeCellTopology(line, second);
      
      if (composite.getBaseKey() != expectedComposite.getBaseKey())
      {
        success = false;
        out << "FAILURE: Expected Line x " << second.getBaseName() << " to be " << expectedComposite.getBaseName();
        out << "; got " << composite.getBaseName() << " instead.\n";
      }
    }
    // blah x line
    std::vector<shards::CellTopology> firstTopos = {line,tri,quad};
    compositeTopos = {quad,wedge,hex};
    for (unsigned i=0; i<firstTopos.size(); i++)
    {
      auto first = firstTopos[i];
      auto expectedComposite = compositeTopos[i];
      auto composite = Intrepid2::TensorTopologyMap::compositeCellTopology(first, line);
      
      if (composite.getBaseKey() != expectedComposite.getBaseKey())
      {
        success = false;
        out << "FAILURE: Expected " << first.getBaseName() << " x Line to be " << expectedComposite.getBaseName();
        out << "; got " << composite.getBaseName() << " instead.\n";
      }
    }
  }

  TEUCHOS_UNIT_TEST( TensorTopology, LineLine )
  {
    testTensorTopologyMapLineLine(out, success);
  }
  
  TEUCHOS_UNIT_TEST( TensorTopology, LineQuad )
  {
    testTensorTopologyMapLineQuad(out, success);
  }
  
  TEUCHOS_UNIT_TEST( TensorTopology, QuadLine )
  {
    testTensorTopologyMapQuadLine(out, success);
  }
  
  TEUCHOS_UNIT_TEST( TensorTopology, CompositeAssignment )
  {
    testTensorTopologyMapCompositeAssignment(out, success);
  }
  
} // namespace
