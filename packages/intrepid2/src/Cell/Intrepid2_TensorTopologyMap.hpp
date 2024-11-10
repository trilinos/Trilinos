// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_TensorTopologyMap.hpp
    \brief  Class that defines mappings from component cell topologies to their tensor product topologies.
    \author Created by N.V. Roberts.
 */

#ifndef Intrepid2_TensorTopologyMap_h
#define Intrepid2_TensorTopologyMap_h

#include "Intrepid2_Utils.hpp"

#include "Shards_CellTopology.hpp"

#include <map>
#include <set>
#include <vector>

namespace Intrepid2
{
  /** \class  Intrepid2::TensorTopologyMap
      \brief  For two cell topologies whose tensor product is a third, this class establishes a mapping from subcell pairs in the component topologies to the tensor product topology.

     This is used in the construction of tensor-product bases; \see Intrepid2::TensorBasis.
   
     This establishes the mapping from component topologies to an existing shards CellTopology, which might not have a "tensor numbering" of subcells
     (the shards quadrilateral is one example: vertices are numbered counter-clockwise rather than as the tensor product of two line topologies).
     */
  class TensorTopologyMap
  {
    shards::CellTopology cellTopo1_;
    shards::CellTopology cellTopo2_;
    shards::CellTopology compositeCellTopo_;
    
    using Subcell     = std::pair<unsigned,unsigned>; // first: dimension, second: subcell ordinal
    using SubcellPair = std::pair<Subcell,Subcell>;   // first: subcell in cellTopo1_, second: subcell in cellTopo2_
    std::map<SubcellPair,Subcell> subcellMap_;        // values are the subcell in compositeCellTopo (the dimension of this subcell is the sum of the dimensions in the SubcellPair)
  public:
    /** \brief  Constructor.
        \param [in] cellTopo1         - the first component cell topology
        \param [in] cellTopo2         - the second component cell topology
        \param [in] compositeCellTopo - the composite cell topology
        \param [in] nodePairs         - ordered list of node pairs that correspond to the nodes in compositeCellTopo
     
     If nodePairs is empty, then TensorTopologyMap will assign a default mapping.  Default mappings are provided for the
     relevant lowest-order (linear) shards topologies in up to three dimensions.
     
     If nodePairs is not empty, then the pair (n1,n2) in position tn indicates a correspondence between the tnth node in
     compositeCellTopo and node n1 in cellTopo1 and node n2 in cellTopo2.
     
     The higher-dimensional subcell mappings are inferred from the node map.
     */
    TensorTopologyMap(shards::CellTopology cellTopo1, shards::CellTopology cellTopo2, shards::CellTopology compositeCellTopo,
                      std::vector< std::pair<unsigned,unsigned> > nodePairs = std::vector< std::pair<unsigned,unsigned> >())
    :
    cellTopo1_(cellTopo1),
    cellTopo2_(cellTopo2),
    compositeCellTopo_(compositeCellTopo)
    {
      if (nodePairs.size() == 0)
      {
        nodePairs = defaultNodePairs(cellTopo1, cellTopo2, compositeCellTopo);
      }
      
      auto spaceDim1 = cellTopo1.getDimension();
      auto spaceDim2 = cellTopo2.getDimension();
      auto spaceDim  = compositeCellTopo.getDimension();
      INTREPID2_TEST_FOR_EXCEPTION(spaceDim1 + spaceDim2 != spaceDim, std::invalid_argument, "incompatible spatial dimensions");
      std::map<std::pair<unsigned,unsigned>,unsigned> compositeNodeOrdinalMap;
      for (unsigned compositeNodeOrdinal=0; compositeNodeOrdinal<nodePairs.size(); compositeNodeOrdinal++)
      {
        compositeNodeOrdinalMap[nodePairs[compositeNodeOrdinal]] = compositeNodeOrdinal;
      }
      for (unsigned d1=0; d1<=spaceDim1; d1++)
      {
        unsigned subcellCount1 = cellTopo1.getSubcellCount(d1);
        for (unsigned subcellOrdinal1=0; subcellOrdinal1<subcellCount1; subcellOrdinal1++)
        {
          Subcell subcell1 = {d1, subcellOrdinal1};
          std::set<unsigned> subcell1Nodes; // set because we don't care about ordering
          unsigned nodeCount1 = cellTopo1_.getNodeCount(d1, subcellOrdinal1);
          // unfortunately, the node count for vertices is given as 0 by shards::CellTopology.  This seems like a bug.
          if (d1 == 0)
          {
            subcell1Nodes.insert(subcellOrdinal1);
          }
          else
          {
            for (unsigned nodeOrdinal1=0; nodeOrdinal1<nodeCount1; nodeOrdinal1++)
            {
              subcell1Nodes.insert(cellTopo1_.getNodeMap(d1, subcellOrdinal1, nodeOrdinal1));
            }
          }
          for (unsigned d2=0; d2<=spaceDim2; d2++)
          {
            unsigned subcellCount2 = cellTopo2.getSubcellCount(d2);
            for (unsigned subcellOrdinal2=0; subcellOrdinal2<subcellCount2; subcellOrdinal2++)
            {
              Subcell subcell2 = {d2, subcellOrdinal2};
              std::set<unsigned> subcell2Nodes; // set because we don't care about ordering
              unsigned nodeCount2 = cellTopo2_.getNodeCount(d2, subcellOrdinal2);
              // unfortunately, the node count for vertices is given as 0 by shards::CellTopology.  This seems like a bug.
              if (d2 == 0)
              {
                subcell2Nodes.insert(subcellOrdinal2);
              }
              else
              {
                for (unsigned nodeOrdinal2=0; nodeOrdinal2<nodeCount2; nodeOrdinal2++)
                {
                  subcell2Nodes.insert(cellTopo2_.getNodeMap(d2, subcellOrdinal2, nodeOrdinal2));
                }
              }
              
              std::set<unsigned> compositeNodes; // all the nodes from subcell1 times nodes from subcell2
              for (auto subcellNode1 : subcell1Nodes)
              {
                for (auto subcellNode2 : subcell2Nodes)
                {
                  INTREPID2_TEST_FOR_EXCEPTION(compositeNodeOrdinalMap.find({subcellNode1,subcellNode2}) == compositeNodeOrdinalMap.end(),
                                             std::invalid_argument, "Node combination not found in map");
                  compositeNodes.insert(compositeNodeOrdinalMap[{subcellNode1,subcellNode2}]);
                }
              }
              // now, search the composite topology for the unique subcell that involves all those nodes
              // we do a brute-force search; we'll never have very big topologies, so the cost is not an issue
              unsigned compositeSubcellDim = d1 + d2;
              unsigned compositeSubcellCount = compositeCellTopo.getSubcellCount(compositeSubcellDim);
              bool compositeSubcellFound = false;
              for (unsigned compositeSubcellOrdinal=0; compositeSubcellOrdinal<compositeSubcellCount; compositeSubcellOrdinal++)
              {
                unsigned compositeSubcellNodeCount = (compositeSubcellDim > 0) ? compositeCellTopo.getNodeCount(compositeSubcellDim, compositeSubcellOrdinal)
                                                                               : 1; // again, dealing with the fact that the node count for vertices is defined as 0, not 1
                if (compositeSubcellNodeCount != compositeNodes.size()) continue; // node counts don't match, so this is not the subcell we're looking for
                bool matches = true; // if we don't find a node that isn't contained in compositeNodes, then this is the subcell we're looking for
                for (unsigned compositeSubcellNode=0; compositeSubcellNode<compositeSubcellNodeCount; compositeSubcellNode++)
                {
                  unsigned nodeInCell = compositeCellTopo.getNodeMap(compositeSubcellDim, compositeSubcellOrdinal, compositeSubcellNode);
                  if (compositeNodes.find(nodeInCell) == compositeNodes.end())
                  {
                    matches = false;
                    break;
                  }
                }
                if (matches)
                {
                  compositeSubcellFound = true;
                  subcellMap_[{subcell1,subcell2}] = {compositeSubcellDim, compositeSubcellOrdinal};
                  break;
                }
              }
              INTREPID2_TEST_FOR_EXCEPTION(!compositeSubcellFound, std::invalid_argument, "Composite subcell not found");
            }
          }
        }
      }
    }
    TensorTopologyMap(shards::CellTopology cellTopo1, shards::CellTopology cellTopo2)
    :
    TensorTopologyMap(cellTopo1, cellTopo2, compositeCellTopology(cellTopo1,cellTopo2)) {}
    
    /** \brief  Map from component subcell ordinals to the corresponding composite subcell ordinal.

        \param  subcell1Dim     [in] - spatial dimension of the subcell in cellTopo1
        \param  subcell1Ordinal [in] - ordinal of the subcell in cellTopo1
        \param  subcell2Dim     [in] - spatial dimension of the subcell in cellTopo2
        \param  subcell2Ordinal [in] - ordinal of the subcell in cellTopo2
     
        \return the subcell ordinal of the corresponding subcell in the composite cell topology.

        The dimension of the composite subcell is subcell1Dim + subcell2Dim.
    */
    unsigned getCompositeSubcellOrdinal(unsigned subcell1Dim, unsigned subcell1Ordinal, unsigned subcell2Dim, unsigned subcell2Ordinal)
    {
      Subcell subcell1 = {subcell1Dim, subcell1Ordinal};
      Subcell subcell2 = {subcell2Dim, subcell2Ordinal};
      auto entry = subcellMap_.find({subcell1,subcell2});
      INTREPID2_TEST_FOR_EXCEPTION(entry == subcellMap_.end(), std::invalid_argument, "entry not found");
      auto subcell = entry->second;
      return subcell.second; // subcell ordinal
    }
    
    /** \brief  A topological representation of the vertex geometries corresponding to a cell topology.

        \param  cellTopo [in] - the cell topology for which the geometry is requested.
     
        \return a vector, whose entries correspond to vertices; each entry is an integer vector of length corresponding to the spatial dimension.

        This method is used internally to TensorTopologyMap to determine default node mappings from component cell topologies to the composite cell topology.
     
       TensorTopologyMap is agnostic about translations and scaling of the reference cell.
       This method gives the cell coordinates as integers, to emphasize that these are not coordinates in reference space,
       but rather a representation of the geometry concerned with x/y/z orientation but little else.
    */
    static std::vector< std::vector<int> > referenceCellGeometry(shards::CellTopology cellTopo)
    {
      std::vector< std::vector<int> > nodes(cellTopo.getVertexCount());
      auto key = cellTopo.getKey();
      switch (key)
      {
        case shards::Node::key:
          nodes.push_back({}); // node.getVertexCount() returns 0; we do want a single (empty) entry, though, even though there's no spatial variation
          break;
        case shards::Line<2>::key:
          nodes[0] = {0};
          nodes[1] = {1};
          break;
        case shards::Triangle<3>::key:
          nodes[0] = {0,0};
          nodes[1] = {1,0};
          nodes[2] = {0,1};
          break;
        case shards::Quadrilateral<4>::key:
          nodes[0] = {0,0};
          nodes[1] = {1,0};
          nodes[2] = {1,1};
          nodes[3] = {0,1};
          break;
        case shards::Hexahedron<8>::key:
          nodes[0] = {0,0,0};
          nodes[1] = {1,0,0};
          nodes[2] = {1,1,0};
          nodes[3] = {0,1,0};
          nodes[4] = {0,0,1};
          nodes[5] = {1,0,1};
          nodes[6] = {1,1,1};
          nodes[7] = {0,1,1};
          break;
        case shards::Wedge<6>::key: // wedge is triangle in (x,y) times line in z
          nodes[0] = {0,0,0};
          nodes[1] = {1,0,0};
          nodes[2] = {0,1,0};
          nodes[3] = {0,0,1};
          nodes[4] = {1,0,1};
          nodes[5] = {0,1,1};
          break;
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "Unsupported CellTopology");
      }
      return nodes;
    }
    
    /** \brief  Default node pairs corresponding to the provided component and composite cell topologies.

        \param  cellTopo1         [in] - the first component cell topology
        \param  cellTopo2         [in] - the second component cell topology
        \param  compositeCellTopo [in] - the composite cell topology
     
        \return ordered list of node pairs that correspond to the nodes in compositeCellTopo

        This method determines a suitable mapping from the component cell topologies to the composite cell topology;
        supported composite topologies are the lowest-order (linear) shards topologies in up to three dimensions, and
        any pair of component cell topologies that generate those composite topologies.  For example, the hexahedron could
        be generated from line x quad, quad x line, node x hexahedron, hexahedron x node.
    */
    static std::vector< std::pair<unsigned,unsigned> > defaultNodePairs(shards::CellTopology cellTopo1, shards::CellTopology cellTopo2, shards::CellTopology compositeCellTopo)
    {
      unsigned compositeNodeCount = compositeCellTopo.getVertexCount();
      unsigned nodeCount1 = cellTopo1.getVertexCount();
      unsigned nodeCount2 = cellTopo2.getVertexCount();
      INTREPID2_TEST_FOR_EXCEPTION(compositeNodeCount != nodeCount1 * nodeCount2, std::invalid_argument, "Incompatible topologies");
      std::vector< std::pair<unsigned,unsigned> > nodePairs(compositeNodeCount);
      auto nodes1 = referenceCellGeometry(cellTopo1);
      auto nodes2 = referenceCellGeometry(cellTopo2);
      auto compositeNodes = referenceCellGeometry(compositeCellTopo);
      std::map< std::vector<int>, unsigned > compositeNodeMap;
      for (unsigned compositeNodeOrdinal=0; compositeNodeOrdinal<compositeNodeCount; compositeNodeOrdinal++)
      {
        compositeNodeMap[compositeNodes[compositeNodeOrdinal]] = compositeNodeOrdinal;
      }
      for (unsigned nodeOrdinal1=0; nodeOrdinal1<nodeCount1; nodeOrdinal1++)
      {
        std::vector<int> node1 = nodes1[nodeOrdinal1];
        for (unsigned nodeOrdinal2=0; nodeOrdinal2<nodeCount2; nodeOrdinal2++)
        {
          std::vector<int> node2 = nodes2[nodeOrdinal2];
          std::vector<int> compositeNode(node1); // copy node1 into the first slots of compositeNode
          compositeNode.insert(compositeNode.end(), node2.begin(), node2.end());
          auto compositeNodeMapEntry = compositeNodeMap.find(compositeNode);
          INTREPID2_TEST_FOR_EXCEPTION(compositeNodeMapEntry == compositeNodeMap.end(), std::invalid_argument, "composite node not found in map");
          nodePairs[compositeNodeMapEntry->second] = {nodeOrdinal1,nodeOrdinal2};
        }
      }
      return nodePairs;
    }
    
    /** \brief  The composite cell topology generated as the tensor product of the specified component cell topologies.

        \param  cellTopo1         [in] - The first component cell topology.
        \param  cellTopo2         [in] - The second component cell topology.
     
        \return The composite cell topology.
    */
    static shards::CellTopology compositeCellTopology(shards::CellTopology cellTopo1, shards::CellTopology cellTopo2)
    {
      if (cellTopo1.getBaseKey() == shards::Node::key)
      {
        return cellTopo2;
      }
      else if (cellTopo2.getBaseKey() == shards::Node::key)
      {
        return cellTopo1;
      }
      
      // if we get here, neither is a Node
      auto key1 = cellTopo1.getBaseKey();
      auto key2 = cellTopo2.getBaseKey();
      switch (key1)
      {
        case shards::Line<2>::key:
          switch (key2)
        {
          case shards::Line<2>::key:
            return shards::CellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >() );
          case shards::Triangle<3>::key:
            // unsupported:
            INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"line x triangle is not supported at present.  Wedge is triangle x line.");
          case shards::Quadrilateral<4>::key:
            return shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >() );
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Unsupported cell topology pair");
        }
        case shards::Triangle<3>::key:
          switch (key2)
        {
          case shards::Line<2>::key:
            return shards::CellTopology(shards::getCellTopologyData<shards::Wedge<> >() );
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Unsupported cell topology pair");
        }
        case shards::Quadrilateral<4>::key:
          switch (key2)
        {
          case shards::Line<2>::key:
            return shards::CellTopology(shards::getCellTopologyData<shards::Hexahedron<> >() );
          default:
            INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Unsupported cell topology pair");
        }
        default:
          INTREPID2_TEST_FOR_EXCEPTION(true,std::invalid_argument,"Unsupported cell topology pair");
      }
    }
  };
} // end namespace Intrepid2

#endif /* Intrepid2_TensorTopology_h */
