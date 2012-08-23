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

#ifndef __Panzer_IntrepidFieldPattern_hpp__
#define __Panzer_IntrepidFieldPattern_hpp__

#include "Panzer_FieldPattern.hpp"

// Trilinos includes
#include "Intrepid_Basis.hpp"
#include "Intrepid_FieldContainer.hpp"

#include "Teuchos_RCP.hpp"

namespace panzer {

/** This is a derived class that specializes based
  * on a single intrepid basis function.
  */
class IntrepidFieldPattern : public FieldPattern {
public:
   IntrepidFieldPattern(const Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > > & intrepidBasis)
      : intrepidBasis_(intrepidBasis) {}

   virtual int getSubcellCount(int dim) const;
   virtual const std::vector<int> & getSubcellIndices(int dim,int cellIndex) const;
   virtual int getDimension() const;
   virtual shards::CellTopology getCellTopology() const;

   virtual void getSubcellClosureIndices(int dim,int cellIndex,std::vector<int> & indices) const;

   // static functions for examining shards objects

   /** For a given sub cell find the set of sub cells at all dimensions contained
     * internally.  This is inclusive, so that (dim,subCell) will be in the set.
     * 
     * \param[in] cellTopo Parent cell topology being used.
     * \param[in] dim      Dimension of sub cell
     * \param[in] subCell  Ordinal of sub cell at specified dimension
     * \param[in,out] closure Set of sub cells associated with specified sub cell.
     *
     * \note Sub cell dimension and ordinals are inserted into <code>closure</code>.
     *       Previous information will not be removed.
     */
   static void buildSubcellClosure(const shards::CellTopology & cellTopo,unsigned dim,unsigned subCell,
                                   std::set<std::pair<unsigned,unsigned> > & closure);

   /** Search a cell topology for sub cells containing a specfic set of nodes.
     * This is a downward search (inclusive) from a user specified dimension.
     *
     * \param[in] cellTopo Parent cell topology being used.
     * \param[in] dim      Dimension of sub cell
     * \param[in] nodes    Nodes forming the super set
     * \param[in,out] subCells Specific sub cells containing the nodes.
     *
     * \note Sub cell dimension and ordinals are inserted into <code>subCells</code>.
     *       Previous information will not be removed.
     */
   static void findContainedSubcells(const shards::CellTopology & cellTopo,unsigned dim,
                                     const std::vector<unsigned> & nodes,
                                     std::set<std::pair<unsigned,unsigned> > & subCells);

   /** Get the set of nodes making up the user specified sub cells.
     *
     * \param[in] cellTopo Parent cell topology being used.
     * \param[in] dim      Dimension of sub cell
     * \param[in] subCell  Ordinal of sub cell at specified dimension
     * \param[in,out] nodes Nodes associated with sub cell.
     */
   static void getSubcellNodes(const shards::CellTopology & cellTopo,unsigned dim,unsigned subCell,
                               std::vector<unsigned> & nodes);

   /** Get the local coordinates for this field. This is independent of element
     * locations.
     *
     * \param[in,out] coords   Coordinates associated with this field type.
     */
   void getInterpolatoryCoordinates(Intrepid::FieldContainer<double> & coords) const;

   /** Get the local coordinates for this field.
     *
     * \param[in] cellVertices   Coordinates associated with this field type.
     * \param[in,out] coords   Coordinates associated with this field type.
     */
   void getInterpolatoryCoordinates(const Intrepid::FieldContainer<double> & cellVertices,
                                    Intrepid::FieldContainer<double> & coords) const;

protected:
   Teuchos::RCP< Intrepid::Basis<double,Intrepid::FieldContainer<double> > >
      intrepidBasis_;
   std::vector<int> empty_;
};

}

#endif 
