#ifndef __Panzer_STK_SquareQuadMeshFactory_hpp__
#define __Panzer_STK_SquareQuadMeshFactory_hpp__

#include <Panzer_STK_MeshFactory.hpp>
#include <Panzer_STK_Interface.hpp>

namespace panzer_stk {

class STK_Interface;

/** This builds a parallel mesh object. Note that the
  * local IDs are ordered by going left to right
  * across the X-axis. 
  */
class LineMeshFactory : public STK_MeshFactory {
public:
   //! Constructor
   LineMeshFactory();

   //! Destructor
   ~LineMeshFactory();

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
   void buildBlock(stk::ParallelMachine machRank,int xBlock,STK_Interface & mesh) const;

   std::pair<int,int> determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int rank) const;

   void addSideSets(STK_Interface & mesh) const;

   // search through relations for the one matching the ID
   const stk::mesh::Relation * getRelationByID(unsigned ID,stk::mesh::PairIterRelation edges) const;

   double x0_;
   double xf_;

   int xBlocks_;

   int nXElems_;

   mutable unsigned int machRank_, machSize_;
   mutable Teuchos::Tuple<std::size_t,2> procTuple_;
};

}

#endif
