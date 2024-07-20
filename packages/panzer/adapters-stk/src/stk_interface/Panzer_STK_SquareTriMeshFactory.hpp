// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_SquareTriMeshFactory_hpp__
#define __Panzer_STK_SquareTriMeshFactory_hpp__

#include <Panzer_STK_MeshFactory.hpp>
#include <Panzer_STK_Interface.hpp>

namespace panzer_stk {

class STK_Interface;

/** This builds a parallel mesh object. Note that the
  * local IDs are ordered by going up the y axis and
  * across the X-axis (in that order). For a mesh with
  * two X blocks and one Y-block, with each block composed
  * of 3x2 (x2) elements the numbering looks like:
  \verbatim
   8  9 10 11
   4  5  6  7
   0  1  2  3
  \endverbatim
  */
class SquareTriMeshFactory : public STK_MeshFactory {
public:
   //! Constructor
   SquareTriMeshFactory();

   //! Destructor
   ~SquareTriMeshFactory();

   //! Build the mesh object
   Teuchos::RCP<STK_Interface> buildMesh(stk::ParallelMachine parallelMach) const;

   virtual Teuchos::RCP<STK_Interface> buildUncommitedMesh(stk::ParallelMachine parallelMach) const;
   virtual void completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const;

   //! From ParameterListAcceptor
   void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList);

   //! From ParameterListAcceptor
   Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

   //! what is the 2D tuple describe this processor distribution
   Teuchos::Tuple<std::size_t,2> procRankToProcTuple(std::size_t procRank) const;

protected:
   void initializeWithDefaults();

   void buildMetaData(stk::ParallelMachine parallelMach,STK_Interface & mesh) const;
   void buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const;
   void buildBlock(stk::ParallelMachine machRank,int xBlock,int yBlock,STK_Interface & mesh) const;

   std::pair<int,int> determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int rank) const;
   std::pair<int,int> determineYElemSizeAndStart(int yBlock,unsigned int size,unsigned int rank) const;

   void addSideSets(STK_Interface & mesh) const;
   void addNodeSets(STK_Interface & mesh) const;
   void addEdgeBlocks(STK_Interface & mesh) const;

   double x0_, y0_;
   double xf_, yf_;

   int xBlocks_, yBlocks_;

   int nXElems_, nYElems_;
   mutable int xProcs_, yProcs_;

   mutable unsigned int machRank_, machSize_;
   mutable Teuchos::Tuple<std::size_t,2> procTuple_;

   bool createEdgeBlocks_;

   std::string edgeBlockName_;
};

}

#endif
