#ifndef __Panzer_STK_CubeHexMeshFactory_hpp__
#define __Panzer_STK_CubeHexMeshFactory_hpp__

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

protected: 
   void initializeWithDefaults();

   void buildMetaData(stk::ParallelMachine parallelMach,STK_Interface & mesh) const;
   void buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const;
   void buildBlock(stk::ParallelMachine machRank,int xBlock,int yBlock,int zBlock,STK_Interface & mesh) const;

   std::pair<int,int> determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int rank) const;
   std::pair<int,int> determineYElemSizeAndStart(int yBlock,unsigned int size,unsigned int rank) const;
   std::pair<int,int> determineZElemSizeAndStart(int zBlock,unsigned int size,unsigned int rank) const;

   void addSideSets(STK_Interface & mesh) const;

   // search through relations for the one matching the ID: for use with addSideSets
   const stk::mesh::Relation * getRelationByID(unsigned ID,stk::mesh::PairIterRelation edges) const;

   double x0_, y0_, z0_;
   double xf_, yf_, zf_;

   int xBlocks_, yBlocks_, zBlocks_;

   int nXElems_, nYElems_, nZElems_;

   mutable unsigned int machRank_, machSize_;
};

}

#endif
