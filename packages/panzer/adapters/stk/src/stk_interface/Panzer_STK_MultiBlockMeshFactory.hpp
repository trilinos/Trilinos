#ifndef __Panzer_STK_MultiBlockMeshFactory_hpp__
#define __Panzer_STK_MultiBlockMeshFactory_hpp__

#include <Panzer_STK_MeshFactory.hpp>
#include <Panzer_STK_Interface.hpp>

namespace panzer_stk {

class STK_Interface;

class MultiBlockMeshFactory : public STK_MeshFactory {
public:
   //! Constructor
   MultiBlockMeshFactory();

   //! Destructor
   ~MultiBlockMeshFactory();

   //! Build the mesh object
   Teuchos::RCP<STK_Interface> buildMesh(stk::ParallelMachine parallelMach) const;

   virtual Teuchos::RCP<STK_Interface> buildUncommitedMesh(stk::ParallelMachine parallelMach) const;
   virtual void completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const;

   //! From ParameterListAcceptor
   void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList);

   //! From ParameterListAcceptor
   Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

   void initializeWithDefaults();

protected: 
   void buildMetaData(stk::ParallelMachine parallelMach, STK_Interface & mesh) const;
   void buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const;
   void buildBlock(stk::ParallelMachine parallelMach,int xBlock,int yBlock,STK_Interface & mesh) const;
   std::pair<int,int> determineXElemSizeAndStart(int xBlock,unsigned int size,unsigned int rank) const;
   std::pair<int,int> determineYElemSizeAndStart(int yBlock,unsigned int size,unsigned int rank) const;
   const stk::mesh::Relation * getRelationByID(unsigned ID,stk::mesh::PairIterRelation relations) const;
   void addSideSets(STK_Interface & mesh) const;

   int nXElems_;
   int nYElems_;
   mutable unsigned int machRank_, machSize_;
};

}

#endif
