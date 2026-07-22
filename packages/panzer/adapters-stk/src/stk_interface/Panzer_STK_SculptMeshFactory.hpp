// @HEADER
// *****************************************************************************
//           Panzer: A partial differential equation assembly
//       engine for strongly coupled complex multiphysics systems
//
// Copyright 2011 NTESS and the Panzer contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef __Panzer_STK_SculptMeshFactory_hpp__
#define __Panzer_STK_SculptMeshFactory_hpp__

#include <Panzer_STK_MeshFactory.hpp>
#include <Panzer_STK_Interface.hpp>


namespace panzer_stk {

class STK_Interface;

/** This builds a parallel mesh object. Sculpt is an all-hex meshing algorithm to handle general solids.
  * Common geometry formats for general solids include STL, SAT, and Diatom.  Major steps in Sculpt include 
  * (1) create a overlay grid  on the general solid
  * (2) remove the grid cells outside the solid
  * (3) project the stair-step grid nodes on the boundary surface, pillow, and smooth.
  */
class SculptMeshFactory : public STK_MeshFactory {
public:
   //! Constructor
   SculptMeshFactory();

   //! Destructor
   ~SculptMeshFactory();

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

   int callSculptor(stk::ParallelMachine parallelMach, char *diatom_file  )const ;
   int writeDiatomFile( std::string stl_path, std::string stl_filename, char *diatom_file ) const;

   void buildNodes(stk::ParallelMachine parallelMach, STK_Interface &mesh ) const;
   void buildMetaData(stk::ParallelMachine parallelMach,STK_Interface & mesh) const;
   void buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const;
   void buildBlock(stk::ParallelMachine parallelMach,STK_Interface & mesh, int block_index, int *block_id, int elem_start, int *elements, int *nodes_per_elem, int *elem_attributes, int **elm_node_linkage ) const;

   std::pair<int,int> determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int rank) const;

   void addSideSets(STK_Interface & mesh) const;
   void addNodeSets(STK_Interface & mesh) const;
   void addEdgeBlocks(STK_Interface & mesh) const;
   void addFaceBlocks(STK_Interface & mesh) const;

   // search through relations for the one matching the ID
   const stk::mesh::Relation * getRelationByID(unsigned ID,stk::mesh::PairIterRelation edges) const;

   std::string stlFileName_, stlFileDir_;
   int xInterval_, yInterval_, zInterval_;
   double xMin_, yMin_, zMin_;
   double xMax_, yMax_, zMax_;

   mutable unsigned int machRank_, machSize_;
   mutable Teuchos::Tuple<std::size_t,2> procTuple_;

};

}

#endif
