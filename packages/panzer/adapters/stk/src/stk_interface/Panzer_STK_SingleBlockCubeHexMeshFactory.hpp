#ifndef __Panzer_STK_SingleBlockCubeHexMeshFactory_hpp__
#define __Panzer_STK_SingleBlockCubeHexMeshFactory_hpp__

#include <Panzer_STK_MeshFactory.hpp>
#include <Panzer_STK_Interface.hpp>

#include <Teuchos_Tuple.hpp>

namespace panzer_stk {

class STK_Interface;

/** This builds a parallel mesh object. Note that the
  * local IDs are ordered by going up the z, then y axis and
  * across the X-axis (in that order). See the SquareQuad mesh
  * factory for more information.
  */
class SingleBlockCubeHexMeshFactory : public STK_MeshFactory {
public:
   //! Constructor
   SingleBlockCubeHexMeshFactory();

   //! Destructor
   virtual ~SingleBlockCubeHexMeshFactory();

   //! Build the mesh object
   Teuchos::RCP<STK_Interface> buildMesh(stk::ParallelMachine parallelMach) const;

   virtual Teuchos::RCP<STK_Interface> buildUncommitedMesh(stk::ParallelMachine parallelMach) const;
   virtual void completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const;

   //! From ParameterListAcceptor
   void setParameterList(const Teuchos::RCP<Teuchos::ParameterList> & paramList);

   //! From ParameterListAcceptor
   Teuchos::RCP<const Teuchos::ParameterList> getValidParameters() const;

   //! Convert processor rank to a tuple
   Teuchos::Tuple<std::size_t,3> procRankToProcTuple(std::size_t procRank) const; 

protected: 
   void initializeWithDefaults();

   void buildMetaData(stk::ParallelMachine parallelMach,STK_Interface & mesh) const;
   void buildElements(stk::ParallelMachine parallelMach,STK_Interface & mesh) const;
   void buildBlock(stk::ParallelMachine machRank,STK_Interface & mesh) const;

   std::pair<int,int> determineXElemSizeAndStart() const;
   std::pair<int,int> determineYElemSizeAndStart() const;
   std::pair<int,int> determineZElemSizeAndStart() const;

   void addSideSets(STK_Interface & mesh) const;

   // search through relations for the one matching the ID: for use with addSideSets
   const stk::mesh::Relation * getRelationByID(unsigned ID,stk::mesh::PairIterRelation edges) const;

   double x0_, y0_, z0_;
   double xf_, yf_, zf_;

   std::size_t xProcs_, yProcs_, zProcs_;

   std::size_t nXElems_, nYElems_, nZElems_;

   mutable unsigned int machRank_, machSize_;

   mutable Teuchos::Tuple<std::size_t,3> procTuple_;
};

}

#endif
