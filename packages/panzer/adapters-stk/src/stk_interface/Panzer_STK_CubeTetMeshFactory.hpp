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

#ifndef __Panzer_STK_CubeTetMeshFactory_hpp__
#define __Panzer_STK_CubeTetMeshFactory_hpp__

#include <Panzer_STK_MeshFactory.hpp>
#include <Panzer_STK_Interface.hpp>

namespace panzer_stk {

class STK_Interface;

/** This builds a parallel mesh object. Note that the
  * local IDs are ordered by going up the z, then y axis and
  * across the X-axis (in that order). See the SquareQuad mesh
  * factory for more information.
  */
class CubeTetMeshFactory : public STK_MeshFactory {
public:
   //! Constructor
   CubeTetMeshFactory();

   //! Destructor
   virtual ~CubeTetMeshFactory();

   //! Build the mesh object
   Teuchos::RCP<STK_Interface> buildMesh(stk::ParallelMachine parallelMach) const;

   virtual Teuchos::RCP<STK_Interface> buildUncommitedMesh(stk::ParallelMachine parallelMach) const;
   virtual void completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const;

   //! From ParameterListAcceptor
   void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList);

   //! From ParameterListAcceptor
   Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

   //! what is the 3D tuple describe this processor distribution
   Teuchos::Tuple<std::size_t,3> procRankToProcTuple(std::size_t procRank) const;

protected: 
   void initializeWithDefaults();

   void buildMetaData(stk::ParallelMachine parallelMach,STK_Interface & mesh) const;
   void buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const;
   void buildBlock(stk::ParallelMachine machRank,int xBlock,int yBlock,int zBlock,STK_Interface & mesh) const;

   std::pair<int,int> determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int rank) const;
   std::pair<int,int> determineYElemSizeAndStart(int yBlock,unsigned int size,unsigned int rank) const;
   std::pair<int,int> determineZElemSizeAndStart(int zBlock,unsigned int size,unsigned int rank) const;

   void addSideSets(STK_Interface & mesh) const;
   void addNodeSets(STK_Interface & mesh) const;

   void buildTetsOnHex(const Teuchos::Tuple<int,3> & meshDesc,
                       const Teuchos::Tuple<int,3> & element,
                       stk::mesh::Part * block,
                       const std::vector<stk::mesh::EntityId> & h_nodes, 
                       STK_Interface & mesh) const;

   double x0_, y0_, z0_;
   double xf_, yf_, zf_;

   int xBlocks_, yBlocks_, zBlocks_;

   int nXElems_, nYElems_, nZElems_;

   mutable int xProcs_, yProcs_, zProcs_;

   mutable unsigned int machRank_, machSize_;

   mutable Teuchos::Tuple<std::size_t,3> procTuple_;
};

}

#endif
