// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_CubeHexMeshFactory_hpp__
#define __Panzer_STK_CubeHexMeshFactory_hpp__

#include <Panzer_Traits.hpp> // for panzer::GlobalOrdinal
#include <Panzer_STK_MeshFactory.hpp>
#include <Panzer_STK_Interface.hpp>

namespace panzer_stk {

class STK_Interface;

/** This builds a parallel mesh object. Note that the
  * local IDs are ordered by going up the z, then y axis and
  * across the X-axis (in that order). See the SquareQuad mesh
  * factory for more information.
  */
class CubeHexMeshFactory : public STK_MeshFactory {
public:
   //! Constructor
   CubeHexMeshFactory();

   //! Destructor
   virtual ~CubeHexMeshFactory();

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

   std::pair<panzer::GlobalOrdinal,panzer::GlobalOrdinal> determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int rank) const;
   std::pair<panzer::GlobalOrdinal,panzer::GlobalOrdinal> determineYElemSizeAndStart(int yBlock,unsigned int size,unsigned int rank) const;
   std::pair<panzer::GlobalOrdinal,panzer::GlobalOrdinal> determineZElemSizeAndStart(int zBlock,unsigned int size,unsigned int rank) const;

   void addSides(STK_Interface & mesh) const; // this adds side entities only (does not inject them into side sets)
   void addSideSets(STK_Interface & mesh) const;
   void addNodeSets(STK_Interface & mesh) const;
   void addEdgeBlocks(STK_Interface & mesh) const;
   void addFaceBlocks(STK_Interface & mesh) const;

   double x0_, y0_, z0_;
   double xf_, yf_, zf_;

   int xBlocks_, yBlocks_, zBlocks_;

   panzer::GlobalOrdinal nXElems_, nYElems_, nZElems_;

   mutable int xProcs_, yProcs_, zProcs_;

   mutable unsigned int machRank_, machSize_;

   bool buildInterfaceSidesets_;
   bool buildSubcells_;
   bool createEdgeBlocks_;
   bool createFaceBlocks_;

   mutable Teuchos::Tuple<std::size_t,3> procTuple_;

   std::string edgeBlockName_;
   std::string faceBlockName_;
};

}

#endif
