// @HEADER
// *****************************************************************************
//                           Intrepid2 Package
//
// Copyright 2007 NTESS and the Intrepid2 contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

/** \file   Intrepid2_CellTopology.hpp
    \brief  Implements arbitrary-dimensional extrusion of a base shards::CellTopology.
 
  This implementation borrows heavily from the Camellia::CellTopology class, which likewise wraps a shards::CellTopology.  Much of the implementation is in fact a simple copy and paste from the Camellia implementation.  The primary distinction between this and what is provided by Camellia is that Camellia additionally includes support for nodal permutations; here, we do not, though we may revisit this in the future.

    \author Nathan V. Roberts
*/

#ifndef Intrepid2_CellTopology_h
#define Intrepid2_CellTopology_h

#include <Shards_CellTopology.hpp>

namespace Intrepid2
{
  /** \class Intrepid2::CellTopology
      \brief Implements arbitrary-dimensional extrusion of a base shards::CellTopology.
  */
  class CellTopology {
  public:
    using CellTopoPtr = Teuchos::RCP<CellTopology>;
    using CellTopologyKey = std::pair<ordinal_type,ordinal_type>;
  protected:
    shards::CellTopology shardsBaseTopology_;
    ordinal_type tensorialDegree_; // number of tensorial extrusions of the base topology
    
    std::string name_;
    
    std::vector< std::vector<CellTopoPtr> > subcells_; // ordered by dimension, then ordinal
  public:
    /**
     \brief Constructor.
     \param [in] baseTopo - the base shards CellTopology.
     \param [in] tensorialDegree - the number of dimensions to extrude into.
     
     The dimension of the CellTopology will be <var>baseTopo.getDimension() + tensorialDegree</var>.
    */
    CellTopology(const shards::CellTopology &baseTopo, ordinal_type tensorialDegree)
    :
    shardsBaseTopology_(baseTopo),
    tensorialDegree_(tensorialDegree)
    {
      using std::vector;
      if (tensorialDegree_ == 0)
      {
        name_ = baseTopo.getName();
      }
      else
      {
        std::ostringstream nameStream;
        nameStream << baseTopo.getName();
        for (int tensorialOrdinal = 0; tensorialOrdinal < tensorialDegree; tensorialOrdinal++)
        {
          nameStream << " x Line_2";
        }
        name_ = nameStream.str();
      }

      int baseDim = baseTopo.getDimension();
      vector<ordinal_type> subcellCounts = vector<ordinal_type>(baseDim + tensorialDegree_ + 1);
      subcells_ = vector< vector< CellTopoPtr > >(baseDim + tensorialDegree_ + 1);

      if (tensorialDegree_==0)
      {
        for (int d=0; d<=baseDim; d++)
        {
          subcellCounts[d] = static_cast<ordinal_type>(baseTopo.getSubcellCount(d));
        }
      }
      else
      {
        CellTopoPtr tensorComponentTopo = getTensorialComponent();
        subcellCounts[0] = 2 * tensorComponentTopo->getSubcellCount(0);
        for (int d=1; d < baseDim+tensorialDegree_; d++)
        {
          subcellCounts[d] = 2 * tensorComponentTopo->getSubcellCount(d) + tensorComponentTopo->getSubcellCount(d-1);
        }
        subcellCounts[baseDim + tensorialDegree_] = 1; // the volume topology
      }
      for (int d=0; d<baseDim+tensorialDegree_; d++)
      {
        subcells_[d] = vector< CellTopoPtr >(subcellCounts[d]);
        int subcellCount = subcells_[d].size();
        for (int scord=0; scord<subcellCount; scord++)
        {
          subcells_[d][scord] = getSubcell(d, scord);
        }
      }
      subcells_[baseDim+tensorialDegree_] = vector<CellTopoPtr>(1);
      subcells_[baseDim+tensorialDegree_][0] = Teuchos::rcp(this, false); // false: does not own memory (self-reference)
    }
    
    /** \brief  Returns the underlying shards CellTopology */
    const shards::CellTopology & getBaseTopology() const
    {
      return shardsBaseTopology_;
    }
    
    /** \brief  The number of times we have taken a tensor product between a line topology and the shards topology to form this cell topology */
    ordinal_type getTensorialDegree() const
    {
      return tensorialDegree_;
    }
    
    /** \brief  Dimension of this tensor topology */
    ordinal_type getDimension() const
    {
      return shardsBaseTopology_.getDimension() + tensorialDegree_;
    }
    
    static ordinal_type getNodeCount(const shards::CellTopology &shardsTopo)
    {
      if (shardsTopo.getDimension()==0) return 1; // Node topology; by my lights shards returns the wrong thing (0) here
      return shardsTopo.getNodeCount();
    }

    /** \brief  Node count of this cell topology */
    ordinal_type getNodeCount() const
    {
      ordinal_type two_pow = 1 << tensorialDegree_;
      return getNodeCount(shardsBaseTopology_) * two_pow;
    }
    
    /** \brief  Vertex count of this cell topology */
    ordinal_type getVertexCount() const
    {
      return getNodeCount();
    }

    /** \brief  Edge (dimension 1) subcell count of this cell topology */
    ordinal_type getEdgeCount() const
    {
      return getSubcellCount(1);
    }

    /** \brief  Face (dimension 2) subcell count of this cell topology */
    ordinal_type getFaceCount() const
    {
      return getSubcellCount(2);
    }

    /** \brief  Side (dimension N-1) subcell count of this cell topology */
    ordinal_type getSideCount() const
    {
      ordinal_type spaceDim = getDimension();
      if (spaceDim == 0)
      {
        return 0;
      }
      else
      {
        int sideDim = spaceDim - 1;
        return getSubcellCount(sideDim);
      }
    }

    /** \brief  Key that's unique for standard shards topologies and any tensorial degree.
     */
    std::pair<ordinal_type,ordinal_type> getKey() const
    {
      return std::make_pair(static_cast<ordinal_type>(shardsBaseTopology_.getKey()), tensorialDegree_);
    }

    /** \brief  Node count of a subcell of the given dimension and ordinal.
     *  \param  subcell_dim    [in]  - spatial dimension of the subcell
     *  \param  subcell_ord    [in]  - subcell ordinal
     */
    ordinal_type getNodeCount( const ordinal_type subcell_dim ,
                               const ordinal_type subcell_ord ) const
    {
      return subcells_[subcell_dim][subcell_ord]->getNodeCount();
    }

    /** \brief  Human-readable name of the CellTopology.
     */
    std::string getName() const
    {
      return name_;
    }

    /** \brief  Vertex count of a subcell of the given dimension and ordinal.
     *  \param  subcell_dim    [in]  - spatial dimension of the subcell
     *  \param  subcell_ord    [in]  - subcell ordinal
     */
    ordinal_type getVertexCount( const ordinal_type subcell_dim ,
                                 const ordinal_type subcell_ord ) const
    {
      return subcells_[subcell_dim][subcell_ord]->getVertexCount();
    }

    /** \brief  Edge count of a subcell of the given dimension and ordinal.
     *  \param  subcell_dim    [in]  - spatial dimension of the subcell
     *  \param  subcell_ord    [in]  - subcell ordinal
     */
    ordinal_type getEdgeCount( const ordinal_type subcell_dim ,
                               const ordinal_type subcell_ord ) const
    {
      return subcells_[subcell_dim][subcell_ord]->getEdgeCount();
    }

    /** \brief  Mapping from the tensorial component CellTopology's subcell ordinal to the corresponding
     *          subcell ordinal of the extruded subcell in the tensor product topology; that is if
     *              this = (shardsTopo x Line_2 x Line_2 ...) x Line_2,
     *          the mapping takes the subcell of dimension subcell_dim_in_component_topo and ordinal subcell_ord_in_component_topo in
     *              (shardsTopo x Line_2 x Line_2 ...)
     *          and returns the ordinal of that subcell extruded in the final Line_2 dimension.
     */
    ordinal_type getExtrudedSubcellOrdinal( const ordinal_type subcell_dim_in_component_topo ,
                                            const ordinal_type subcell_ord_in_component_topo ) const
    {
      // The rule is that the two copies (unextruded) of subcells of dimension
      // (subcell_dim_in_component_topo + 1) come first, and then the ones built from the extrusion of
      // subcells of dimension subcell_dim_in_component_topo in the component topology.
      if (tensorialDegree_==0)
      {
        INTREPID2_TEST_FOR_EXCEPTION(true, std::invalid_argument, "getExtrudedSubcellOrdinal() is not valid for un-extruded topologies");
      }
      else
      {
        ordinal_type componentSubcellCount = getTensorialComponent()->getSubcellCount(subcell_dim_in_component_topo + 1);
        return subcell_ord_in_component_topo + componentSubcellCount * 2;
      }
    }

    /** \brief  Side count of a subcell of the given dimension and ordinal.
     *  \param  subcell_dim    [in]  - spatial dimension of the subcell
     *  \param  subcell_ord    [in]  - subcell ordinal
     */
    ordinal_type getSideCount( const ordinal_type subcell_dim ,
                              const ordinal_type subcell_ord ) const
    {
      return subcells_[subcell_dim][subcell_ord]->getSideCount();
    }


    /** \brief  Subcell count of subcells of the given dimension.
     *  \param  subcell_dim    [in]  - spatial dimension of the subcell
     */
    ordinal_type getSubcellCount( const ordinal_type subcell_dim ) const
    {
      if (subcell_dim >= ordinal_type(subcells_.size())) return 0;
      else return subcells_[subcell_dim].size();
    }

    /** \brief  Mapping from the tensorial component node ordinals to the
     *          node ordinal of this tensor cell topology.
     *  \param  tensorComponentNodes      [in]  - node ordinals in the tensorial components.
     */
    ordinal_type getNodeFromTensorialComponentNodes(const std::vector<ordinal_type> &tensorComponentNodes) const
    {
      if (ordinal_type(tensorComponentNodes.size()) != tensorialDegree_ + 1)
      {
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "tensorComponentNodes.size() != _tensorialDegree + 1");
      }
      /*
       Example: we have a base topology of 4 nodes x line x line.  Read addresses from right to left.

       address (1,0,0) --> 0 * (2 * 4) + 0 * 4 + 1 =  1
       address (0,1,0) --> 0 * (2 * 4) + 1 * 4 + 0 =  4
       address (0,0,1) --> 1 * (2 * 4) + 0 * 4 + 0 =  8
       address (0,1,1) --> 1 * (2 * 4) + 1 * 4 + 0 = 12

       */

      ordinal_type node = 0;
      CellTopoPtr line = CellTopology::line();
      std::vector<CellTopoPtr> componentTopos(tensorialDegree_ + 1, line);
      componentTopos[0] = cellTopology(shardsBaseTopology_);
      for (int i=tensorComponentNodes.size()-1; i >= 0; i--)
      {
        ordinal_type componentNode = tensorComponentNodes[i];
        node *= componentTopos[i]->getNodeCount();
        node += componentNode;
      }
      return node;
    }

    /** \brief  Mapping from a subcell's node ordinal to a
     *          node ordinal of this parent cell topology.
     *  \param  subcell_dim      [in]  - spatial dimension of the subcell
     *  \param  subcell_ord      [in]  - subcell ordinal
     *  \param  subcell_node_ord [in]  - node ordinal relative to subcell
     */
    ordinal_type getNodeMap( const ordinal_type subcell_dim ,
                             const ordinal_type subcell_ord ,
                             const ordinal_type subcell_node_ord ) const
    {
      if (subcell_dim == getDimension())
      {
        // map from topology to itself
        if (subcell_ord != 0)
        {
          TEUCHOS_TEST_FOR_EXCEPTION(subcell_ord != 0, std::invalid_argument, "subcell ordinal out of bounds");
        }
        return subcell_node_ord;
      }
      else if (subcell_dim==0)
      {
        // mapping a node--the subcell_node_ord must be 0, then, and we should just return the subcell_ord (which is the node ordinal)
        if (subcell_node_ord != 0)
        {
          TEUCHOS_TEST_FOR_EXCEPTION(subcell_node_ord != 0, std::invalid_argument, "subcell node ordinal out of bounds");
        }
        return subcell_ord;
      }
      if (tensorialDegree_==0)
      {
        return shardsBaseTopology_.getNodeMap(subcell_dim, subcell_ord, subcell_node_ord);
      }
      else
      {
        CellTopoPtr tensorComponentTopo = CellTopology::cellTopology(shardsBaseTopology_, tensorialDegree_ - 1);
        ordinal_type componentSubcellCount = tensorComponentTopo->getSubcellCount(subcell_dim);
        if (subcell_ord < componentSubcellCount * 2)   // subcell belongs to one of the two component topologies
        {
          ordinal_type subcell_ord_comp = subcell_ord % componentSubcellCount;  // subcell ordinal in the component topology
          ordinal_type compOrdinal = subcell_ord / componentSubcellCount; // which component topology? 0 or 1.
          ordinal_type mappedNodeInsideComponentTopology = tensorComponentTopo->getNodeMap(subcell_dim, subcell_ord_comp, subcell_node_ord);
          return mappedNodeInsideComponentTopology + compOrdinal * tensorComponentTopo->getNodeCount();
        }
        else
        {
          // otherwise, the subcell is a tensor product of a component's (subcell_dim-1)-dimensional subcell with the line topology.
          ordinal_type subcell_ord_comp = subcell_ord - componentSubcellCount * 2;
          ordinal_type subcell_dim_comp = subcell_dim - 1;
          CellTopoPtr subcellTensorComponent = tensorComponentTopo->getSubcell(subcell_dim_comp, subcell_ord_comp);
          // which of the two copies of the subcell tensor component owns the node subcell_node_ord?
          ordinal_type scCompOrdinal = subcell_node_ord / subcellTensorComponent->getNodeCount(); // 0 or 1
          // what's the node ordinal inside the subcell component?
          ordinal_type scCompNodeOrdinal = subcell_node_ord % subcellTensorComponent->getNodeCount();
          ordinal_type mappedNodeInsideComponentTopology = tensorComponentTopo->getNodeMap(subcell_dim_comp, subcell_ord_comp, scCompNodeOrdinal);
          return mappedNodeInsideComponentTopology + scCompOrdinal * tensorComponentTopo->getNodeCount();
        }
      }
    }

    /** \brief  Number of node permutations defined for this cell */
    ordinal_type getNodePermutationCount() const;

    /** \brief  Permutation of a cell's node ordinals.
     *  \param  permutation_ordinal [in]
     *  \param  node_ordinal        [in]
     */
    ordinal_type getNodePermutation( const ordinal_type permutation_ord ,
                                 const ordinal_type node_ord ) const;

    /** \brief  Inverse permutation of a cell's node ordinals.
     *  \param  permutation_ordinal [in]
     *  \param  node_ordinal        [in]
     */
    ordinal_type getNodePermutationInverse( const ordinal_type permutation_ord ,
                                        const ordinal_type node_ord ) const;

    /** \brief  Returns a CellTopoPtr for the specified side.
     */
    CellTopoPtr getSide( ordinal_type sideOrdinal ) const;

    /** \brief  Get the subcell of dimension scdim with ordinal scord.
     *  \param  scdim        [in]
     *  \param  scord        [in]
     *  For tensor-product topologies T x L (L being the line topology), there are two "copies" of T, T0 and T1,
     *  and the enumeration of subcells of dimension d goes as follows:
     - d-dimensional subcells from T0
     - d-dimensional subcells from T1
     - ((d-1)-dimensional subcells of T) x L.
     */
    CellTopoPtr getSubcell( ordinal_type scdim, ordinal_type scord ) const
    {
      if (tensorialDegree_==0)
      {
        return cellTopology(shardsBaseTopology_.getCellTopologyData(scdim, scord), 0);
      }
      else
      {
        CellTopoPtr tensorComponentTopo = getTensorialComponent();
        ordinal_type componentSubcellCount = tensorComponentTopo->getSubcellCount(scdim);
        if (scord < componentSubcellCount * 2)
        {
          scord = scord % componentSubcellCount;
          return tensorComponentTopo->getSubcell(scdim, scord);
        }
        // otherwise, the subcell is a tensor product of one of the components (scdim-1)-dimensional subcells with the line topology.
        scord = scord - componentSubcellCount * 2;
        scdim = scdim - 1;
        CellTopoPtr subcellTensorComponent = tensorComponentTopo->getSubcell(scdim, scord);
        return cellTopology(subcellTensorComponent->getBaseTopology(), subcellTensorComponent->getTensorialDegree() + 1);
      }
    }
    
    /** \brief  Maps the from a subcell within a subcell of the present CellTopology to the subcell in the present CellTopology; returns the ordinal of the subcell's subcell within the cell.
     
     Adapted from CamelliaCellTools::getSubcellOrdinalMap().
     */
    static ordinal_type getSubcellOrdinalMap(CellTopoPtr cellTopo, ordinal_type subcdim, ordinal_type subcord, ordinal_type subsubcdim, ordinal_type subsubcord)
    {
      // static map<> …
      // TODO: implement this method
      
      // maps from a subcell's ordering of its subcells (the sub-subcells) to the cell topology's ordering of those subcells.
      typedef ordinal_type SubcellOrdinal;
      typedef ordinal_type SubcellDimension;
      typedef ordinal_type SubSubcellOrdinal;
      typedef ordinal_type SubSubcellDimension;
      typedef ordinal_type SubSubcellOrdinalInCellTopo;
      typedef std::pair< SubcellDimension, SubcellOrdinal > SubcellIdentifier; // dim, ord in cellTopo
      typedef std::pair< SubSubcellDimension, SubSubcellOrdinal > SubSubcellIdentifier; // dim, ord in subcell
      typedef std::map< SubcellIdentifier, std::map< SubSubcellIdentifier, SubSubcellOrdinalInCellTopo > > OrdinalMap;
      static std::map< CellTopologyKey, OrdinalMap > ordinalMaps;

      if (subsubcdim==subcdim)
      {
        if (subsubcord==0)   // i.e. the "subsubcell" is really just the subcell
        {
          return subcord;
        }
        else
        {
          std::cout << "request for subsubcell of the same dimension as subcell, but with subsubcord > 0.\n";
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "request for subsubcell of the same dimension as subcell, but with subsubcord > 0.");
        }
      }

      if (subcdim==cellTopo->getDimension())
      {
        if (subcord==0)   // i.e. the subcell is the cell itself
        {
          return subsubcord;
        }
        else
        {
          std::cout << "request for subcell of the same dimension as cell, but with subsubcord > 0.\n";
          TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "request for subcell of the same dimension as cell, but with subsubcord > 0.");
        }
      }

      CellTopologyKey key = cellTopo->getKey();
      if (ordinalMaps.find(key) == ordinalMaps.end())
      {
        // then we construct the map for this cellTopo
        OrdinalMap ordinalMap;
        ordinal_type sideDim = cellTopo->getDimension() - 1;
        typedef ordinal_type NodeOrdinal;
        std::map< std::set<NodeOrdinal>, SubcellIdentifier > subcellMap; // given set of nodes in cellTopo, what subcell is it?)

        for (ordinal_type d=1; d<=sideDim; d++)   // only things of dimension >= 1 will have subcells
        {
          ordinal_type subcellCount = cellTopo->getSubcellCount(d);
          for (ordinal_type subcellOrdinal=0; subcellOrdinal<subcellCount; subcellOrdinal++)
          {
            std::set<NodeOrdinal> nodes;
            ordinal_type nodeCount = cellTopo->getNodeCount(d, subcellOrdinal);
            for (NodeOrdinal subcNode=0; subcNode<nodeCount; subcNode++)
            {
              nodes.insert(cellTopo->getNodeMap(d, subcellOrdinal, subcNode));
            }
            SubcellIdentifier subcell = std::make_pair(d, subcellOrdinal);
            subcellMap[nodes] = subcell;

            CellTopoPtr subcellTopo = cellTopo->getSubcell(d, subcellOrdinal);
            // now, go over all the subsubcells, and look them up...
            for (ordinal_type subsubcellDim=0; subsubcellDim<d; subsubcellDim++)
            {
              ordinal_type subsubcellCount = subcellTopo->getSubcellCount(subsubcellDim);
              for (ordinal_type subsubcellOrdinal=0; subsubcellOrdinal<subsubcellCount; subsubcellOrdinal++)
              {
                SubSubcellIdentifier subsubcell = std::make_pair(subsubcellDim,subsubcellOrdinal);
                if (subsubcellDim==0)   // treat vertices separately
                {
                  ordinalMap[subcell][subsubcell] = cellTopo->getNodeMap(subcell.first, subcell.second, subsubcellOrdinal);
                  continue;
                }
                ordinal_type nodeCount_inner = subcellTopo->getNodeCount(subsubcellDim, subsubcellOrdinal);
                std::set<NodeOrdinal> subcellNodes; // NodeOrdinals index into cellTopo, though!
                for (NodeOrdinal subsubcNode=0; subsubcNode<nodeCount_inner; subsubcNode++)
                {
                  NodeOrdinal subcNode = subcellTopo->getNodeMap(subsubcellDim, subsubcellOrdinal, subsubcNode);
                  NodeOrdinal node = cellTopo->getNodeMap(d, subcellOrdinal, subcNode);
                  subcellNodes.insert(node);
                }

                SubcellIdentifier subsubcellInCellTopo = subcellMap[subcellNodes];
                ordinalMap[ subcell ][ subsubcell ] = subsubcellInCellTopo.second;
                //              cout << "ordinalMap( (" << subcell.first << "," << subcell.second << "), (" << subsubcell.first << "," << subsubcell.second << ") ) ";
                //              cout << " ---> " << subsubcellInCellTopo.second << endl;
              }
            }
          }
        }
        ordinalMaps[key] = ordinalMap;
      }
      SubcellIdentifier subcell = std::make_pair(subcdim, subcord);
      SubSubcellIdentifier subsubcell = std::make_pair(subsubcdim, subsubcord);
      if (ordinalMaps[key][subcell].find(subsubcell) != ordinalMaps[key][subcell].end())
      {
        return ordinalMaps[key][subcell][subsubcell];
      }
      else
      {
        std::cout << "For topology " << cellTopo->getName() << " and subcell " << subcord << " of dim " << subcdim;
        std::cout << ", subsubcell " << subsubcord << " of dim " << subsubcdim << " not found.\n";
        TEUCHOS_TEST_FOR_EXCEPTION(true, std::invalid_argument, "subsubcell not found");
        return -1; // NOT FOUND
      }
    }

    /** \brief  For cell topologies of positive tensorial degree, returns the cell topology of tensorial degree one less.
                For cell topologies of tensorial degree zero, returns Teuchos::null.
     */
    CellTopoPtr getTensorialComponent() const
    {
      if (tensorialDegree_ > 0)
      {
        return cellTopology(shardsBaseTopology_, tensorialDegree_ - 1);
      }
      else
      {
        return Teuchos::null;
      }
    }
    
    /** \brief  Returns the side corresponding to the provided node in the final extrusion dimension.
     *  \param  extrusionNodeOrdinal      [in]  - 0 or 1, the node number for the vertex in the extrusion line topology.
     */
    ordinal_type getTensorialComponentSideOrdinal(ordinal_type extrusionNodeOrdinal)
    {
      // our ordering places the "copies" of the tensorial component first, so these have side ordinal 0 or 1.
      return extrusionNodeOrdinal;
    }

    /** \brief  Returns true if the specified side has extension in the final tensorial dimension.  For topologies with zero tensorialDegree_, returns false.
     *  \param  sideOrdinal [in] Ordinal of the side.
     */
    bool sideIsExtrudedInFinalDimension( ordinal_type sideOrdinal ) const
    {
      int sideCount = getSideCount();
      if (tensorialDegree_ == 0) return false;
      
      return (sideOrdinal > 1) && (sideOrdinal < sideCount);
    }
    
    /** \brief  static accessor that returns a CellTopoPtr; these are lazily constructed and cached.
     *  \param [in] baseTopo - the base shards CellTopology.
     *  \param [in] tensorialDegree - the number of dimensions to extrude into.
     */
    static CellTopoPtr cellTopology(const shards::CellTopology &shardsCellTopo, ordinal_type tensorialDegree = 0)
    {
      ordinal_type shardsKey = static_cast<ordinal_type>(shardsCellTopo.getBaseKey());
      std::pair<ordinal_type,ordinal_type> key = std::make_pair(shardsKey, tensorialDegree);
      
      static std::map< CellTopologyKey, CellTopoPtr > tensorizedShardsTopologies; // (shards key, n) --> our CellTopoPtr for that cellTopo's nth-order tensor product with a line topology.  I.e. (shard::CellTopology::Line<2>::key, 2) --> a tensor-product hexahedron.  (This differs from the Shards hexahedron, because the enumeration of the sides of the quad in Shards goes counter-clockwise.)

      if (tensorizedShardsTopologies.find(key) == tensorizedShardsTopologies.end())
      {
        tensorizedShardsTopologies[key] = Teuchos::rcp( new CellTopology(shardsCellTopo, tensorialDegree));
      }
      return tensorizedShardsTopologies[key];
    }
    
    static CellTopoPtr point()
    {
      return cellTopology(shards::getCellTopologyData<shards::Node >());
    }
    
    static CellTopoPtr line()
    {
      return cellTopology(shards::getCellTopologyData<shards::Line<> >());
    }
    
    static CellTopoPtr quad()
    {
      return cellTopology(shards::getCellTopologyData<shards::Quadrilateral<> >());
    }
    
    static CellTopoPtr hexahedron()
    {
      return cellTopology(shards::getCellTopologyData<shards::Hexahedron<> >());
    }

    static CellTopoPtr triangle()
    {
      return cellTopology(shards::getCellTopologyData<shards::Triangle<> >());
    }
    
    static CellTopoPtr tetrahedron()
    {
      return cellTopology(shards::getCellTopologyData<shards::Tetrahedron<> >());
    }
    
    static CellTopoPtr wedge()
    {
      return cellTopology(shards::getCellTopologyData<shards::Wedge<> >());
    }
  };
}

#endif /* Intrepid2_CellTopology_h */
