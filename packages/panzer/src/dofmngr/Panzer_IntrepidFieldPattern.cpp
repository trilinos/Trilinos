// @HEADER
// ***********************************************************************
//
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//                 Copyright (2011) Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roger P. Pawlowski (rppawlo@sandia.gov) and
// Eric C. Cyr (eccyr@sandia.gov)
// ***********************************************************************
// @HEADER

#include "Panzer_IntrepidFieldPattern.hpp"

#include "Intrepid_CellTools.hpp"
#include "Shards_CellTopology.hpp"

namespace panzer {

int IntrepidFieldPattern::getSubcellCount(int dim) const
{
   const shards::CellTopology ct = intrepidBasis_->getBaseCellTopology();
   return ct.getSubcellCount(dim);
}

const std::vector<int> & IntrepidFieldPattern::getSubcellIndices(int dim,int cellIndex) const
{
   const std::vector<std::vector<std::vector<int> > > & ordData = intrepidBasis_->getDofOrdinalData();

   // enough dimemnsions from intrepid?
   if((std::size_t) dim<ordData.size()) {
      // enough cells from intrepid?
      if((std::size_t) cellIndex<ordData[dim].size()) {
         if(ordData[dim][cellIndex].size()!=1)
            return ordData[dim][cellIndex];
         else if(ordData[dim][cellIndex][0]<0)
            return empty_; 
         else
            return ordData[dim][cellIndex];
      }
   }
   
   // not enough information!
   return empty_;
}

void IntrepidFieldPattern::getSubcellClosureIndices(int dim,int cellIndex,std::vector<int> & indices) const
{
   // recursive base case
   if(dim==0) {
      indices = getSubcellIndices(dim,cellIndex);
      return;
   }

   indices.clear(); // wipe out previous values...we are using insert!

   // use topology to build closure of sub cell
   const shards::CellTopology ct = intrepidBasis_->getBaseCellTopology();
   std::set<std::pair<unsigned,unsigned> > closure;
   IntrepidFieldPattern::buildSubcellClosure(ct,dim,cellIndex,closure);

   // grab basis indices on the closure of the sub cell
   std::set<std::pair<unsigned,unsigned> >::const_iterator itr;
   for(itr=closure.begin();itr!=closure.end();++itr) {
      // grab indices for this sub cell
      const std::vector<int> & subcellIndices = getSubcellIndices(itr->first,itr->second);

      // concatenate the two vectors
      indices.insert(indices.end(),subcellIndices.begin(),subcellIndices.end());
   }
}

int IntrepidFieldPattern::getDimension() const
{
   const shards::CellTopology ct = intrepidBasis_->getBaseCellTopology();
   return ct.getDimension();
}

shards::CellTopology IntrepidFieldPattern::getCellTopology() const
{
   return intrepidBasis_->getBaseCellTopology();
}

void IntrepidFieldPattern::getSubcellNodes(const shards::CellTopology & cellTopo,unsigned dim,unsigned subCell,
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

void IntrepidFieldPattern::findContainedSubcells(const shards::CellTopology & cellTopo,unsigned dim,
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
   if(dim==0) return;

   // find subcells in next sub dimension
   findContainedSubcells(cellTopo,dim-1,nodes,subCells);
}

void IntrepidFieldPattern::buildSubcellClosure(const shards::CellTopology & cellTopo,unsigned dim,unsigned subCell,
                                               std::set<std::pair<unsigned,unsigned> > & closure)
{
#if 0
   // get all nodes on requested sub cell
   std::vector<unsigned> nodes;
   getSubcellNodes(cellTopo,dim,subCell,nodes);

   closure.insert(std::make_pair(dim,subCell));
   findContainedSubcells(cellTopo,dim-1,nodes,closure);
#else
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
#endif
}

/** Get the local coordinates for this field. This is independent of element
  * locations.
  *
  * \param[in,out] coords   Coordinates associated with this field type.
  */
void IntrepidFieldPattern::getInterpolatoryCoordinates(Intrepid::FieldContainer<double> & coords) const
{
   typedef Intrepid::DofCoordsInterface<Intrepid::FieldContainer<double> > CoordsInterface;

   using Teuchos::RCP;
   using Teuchos::rcp_dynamic_cast;

   bool throwOnFail = true;

   // cast basis object to DofCoordsInterface: throw on failure
   RCP<CoordsInterface> coordsInterface
         = rcp_dynamic_cast<CoordsInterface>(intrepidBasis_,throwOnFail);

   // resize coordinates
   coords.resize(intrepidBasis_->getCardinality(),getDimension());
   coordsInterface->getDofCoords(coords);
}

/** Get the local coordinates for this field. This is independent of element
  * locations.
  *
  * \param[in,out] coords   Coordinates associated with this field type.
  */
void IntrepidFieldPattern::getInterpolatoryCoordinates(const Intrepid::FieldContainer<double> & cellVertices,
                                                       Intrepid::FieldContainer<double> & coords) const
{
   TEUCHOS_ASSERT(cellVertices.rank()==3);

   int numCells = cellVertices.dimension(0);

   // grab the local coordinates
   Intrepid::FieldContainer<double> localCoords;
   getInterpolatoryCoordinates(localCoords);

   // resize the coordinates field container
   coords.resize(numCells,localCoords.dimension(0),getDimension());

   if(numCells>0) {
      // map to phsyical coordinates
      Intrepid::CellTools<double> cellTools;
      cellTools.mapToPhysicalFrame(coords,localCoords,cellVertices,intrepidBasis_->getBaseCellTopology());
   }
}

}
