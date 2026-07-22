// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include "Panzer_IntrepidFieldPattern.hpp"

#include "Teuchos_Assert.hpp"
#include "Intrepid2_CellTools.hpp"
#include "Shards_CellTopology.hpp"

namespace panzer {

  Intrepid2FieldPattern::
  Intrepid2FieldPattern(const Teuchos::RCP< Intrepid2::Basis<PHX::Device,double,double> > &intrepidBasis)
    : intrepidBasis_(intrepidBasis) {
    const auto dofOrd = intrepidBasis_->getAllDofOrdinal(); // rank 3 view
    const auto dofTag = intrepidBasis_->getAllDofTags(); // rank 2 view
    
    const int 
      iend = dofOrd.extent(0),
      jend = dofOrd.extent(1);
    
    subcellIndicies_.resize(iend);
    for (int i=0;i<iend;++i) {
      subcellIndicies_[i].resize(jend);
      for (int j=0;j<jend;++j) {
        const int ord = dofOrd(i, j, 0);
        if (ord >= 0) {
          const int ndofs = dofTag(ord, 3);
          subcellIndicies_[i][j].resize(ndofs); 
          for (int k=0;k<ndofs;++k) 
            subcellIndicies_[i][j][k] = dofOrd(i, j, k);
        } else {
          // if ordinal does not exist empty container.
          subcellIndicies_[i][j].clear();
        }
      }
    }
  }
  
  int 
  Intrepid2FieldPattern::
  getSubcellCount(int dim) const
  {
    const shards::CellTopology ct = intrepidBasis_->getBaseCellTopology();
    return ct.getSubcellCount(dim);
  }

  const std::vector<int> &
  Intrepid2FieldPattern::getSubcellIndices(int dim, int cellIndex) const
  {
    if ((dim       < static_cast<int>(subcellIndicies_.size()     ))  and 
        (cellIndex < static_cast<int>(subcellIndicies_[dim].size())))
      return subcellIndicies_[dim][cellIndex];

    return empty_;
  }

  void 
  Intrepid2FieldPattern::
  getSubcellClosureIndices(int dim,int cellIndex,std::vector<int> & indices) const
  {
    // wipe out previous values
    indices.clear(); 
    
    if (dim == 0) {
      indices = getSubcellIndices(dim,cellIndex);
    } else {
      // construct full topology
      const shards::CellTopology ct = intrepidBasis_->getBaseCellTopology();

      std::set<std::pair<unsigned,unsigned> > closure;
      Intrepid2FieldPattern::buildSubcellClosure(ct,dim,cellIndex,closure);
      
      // grab basis indices on the closure of the sub cell
      std::set<std::pair<unsigned,unsigned> >::const_iterator itr;
      for (itr=closure.begin();itr!=closure.end();++itr) {
        // grab indices for this sub cell
        const std::vector<int> & subcellIndices = getSubcellIndices(itr->first,itr->second);

        // concatenate the two vectors
        indices.insert(indices.end(),subcellIndices.begin(),subcellIndices.end());
      }
    }
  }
  
  int 
  Intrepid2FieldPattern::
  getDimension() const
  {
    return intrepidBasis_->getBaseCellTopology().getDimension();
  }

  shards::CellTopology 
  Intrepid2FieldPattern::
  getCellTopology() const
  {
    // TODO BWR Probably should change the name of this call...
    // TODO BWR However this is a virtual function
    return intrepidBasis_->getBaseCellTopology();
  }

  void 
  Intrepid2FieldPattern::
  getSubcellNodes(const shards::CellTopology & cellTopo,unsigned dim,unsigned subCell,
                  std::vector<unsigned> & nodes)
  {
    if(dim==0) {
      nodes.push_back(subCell);
      return;
    }

    // get all nodes on requested sub cell
    unsigned subCellNodeCount = cellTopo.getNodeCount(dim,subCell);
    for(unsigned node=0;node<subCellNodeCount;++node)
      nodes.push_back(cellTopo.getNodeMap(dim,subCell,node));

    // sort them so they are ordered correctly for "includes" call
    std::sort(nodes.begin(),nodes.end());
  }

  void 
  Intrepid2FieldPattern::
  findContainedSubcells(const shards::CellTopology & cellTopo,unsigned dim,
                        const std::vector<unsigned> & nodes,
                        std::set<std::pair<unsigned,unsigned> > & subCells)
  {
    unsigned subCellCount = cellTopo.getSubcellCount(dim); 
    for(unsigned subCellOrd=0;subCellOrd<subCellCount;++subCellOrd) {
      // get all nodes in sub cell
      std::vector<unsigned> subCellNodes;
      getSubcellNodes(cellTopo,dim,subCellOrd,subCellNodes);

      // if subCellNodes \subset nodes => add (dim,subCellOrd) to subCells
      bool isSubset = std::includes(       nodes.begin(),        nodes.end(),
                                           subCellNodes.begin(), subCellNodes.end());
      if(isSubset)
        subCells.insert(std::make_pair(dim,subCellOrd));
       
    }

    // stop recursion base case
    if (dim==0) return;

    // find subcells in next sub dimension
    findContainedSubcells(cellTopo,dim-1,nodes,subCells);
  }

  void 
  Intrepid2FieldPattern::
  buildSubcellClosure(const shards::CellTopology & cellTopo,unsigned dim,unsigned subCell,
                      std::set<std::pair<unsigned,unsigned> > & closure)
  {
    switch(dim) {
    case 0:
      closure.insert(std::make_pair(0,subCell));
      break;
    case 1:
      closure.insert(std::make_pair(0,cellTopo.getNodeMap(dim,subCell,0)));
      closure.insert(std::make_pair(0,cellTopo.getNodeMap(dim,subCell,1)));
      closure.insert(std::make_pair(1,subCell));
      break;
    case 2:
      {
        unsigned cnt = (shards::CellTopology(cellTopo.getCellTopologyData(dim,subCell))).getSubcellCount(dim-1);
        for(unsigned i=0;i<cnt;i++) {
          int edge = mapCellFaceEdge(cellTopo.getCellTopologyData(),subCell,i);
          buildSubcellClosure(cellTopo,dim-1,edge,closure);
        }
        closure.insert(std::make_pair(2,subCell));
      }
      break;
    default:
      // beyond a two dimension surface this thing crashes!
      TEUCHOS_ASSERT(false);
    };
  }

  bool 
  Intrepid2FieldPattern::
  supportsInterpolatoryCoordinates() const
  {
    // we no longer use CoordsInterface
    return true;
  }
  
  /** Get the local coordinates for this field. This is independent of element
   * locations.
   *
   * \param[in,out] coords   Coordinates associated with this field type.
   */
  void 
  Intrepid2FieldPattern::
  getInterpolatoryCoordinates(Kokkos::DynRankView<double,PHX::Device> & coords) const
  {
    // this may not be efficient if coords is allocated every time this function is called
    coords = Kokkos::DynRankView<double,PHX::Device>("coords",intrepidBasis_->getCardinality(),getDimension());
    intrepidBasis_->getDofCoords(coords);
  }

  /** Get the local coordinates for this field. This is independent of element
   * locations.
   *
   * \param[in,out] coords   Coordinates associated with this field type.
   */
  void 
  Intrepid2FieldPattern::
  getInterpolatoryCoordinates(const Kokkos::DynRankView<double,PHX::Device> & cellNodes,
                              Kokkos::DynRankView<double,PHX::Device> & coords,
                              Teuchos::RCP<const shards::CellTopology> meshCellTopology) const
  {
    TEUCHOS_ASSERT(cellNodes.rank()==3);

    int numCells = cellNodes.extent(0);

    // grab the local coordinates
    Kokkos::DynRankView<double,PHX::Device> localCoords;
    getInterpolatoryCoordinates(localCoords);

    // resize the coordinates field container
    coords = Kokkos::DynRankView<double,PHX::Device>("coords",numCells,localCoords.extent(0),getDimension());

    if(numCells>0) {
      Intrepid2::CellTools<PHX::Device> cellTools;
      // For backwards compatability, allow the FE basis to supply the mesh cell topology (via the FEM base cell topo)
      // This will occur if no meshCellTopology is provided (defaults to Teuchos::null)
      // If provided, use the mesh topology directly
      if (meshCellTopology==Teuchos::null) {
        cellTools.mapToPhysicalFrame(coords,localCoords,cellNodes,intrepidBasis_->getBaseCellTopology());
      } else {
        cellTools.mapToPhysicalFrame(coords,localCoords,cellNodes,*meshCellTopology);
      }
    }
  }

  Teuchos::RCP< Intrepid2::Basis<PHX::Device,double,double> >
  Intrepid2FieldPattern::getIntrepidBasis() const
  { return intrepidBasis_; }
  
}
