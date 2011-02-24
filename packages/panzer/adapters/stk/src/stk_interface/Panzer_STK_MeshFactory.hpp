#ifndef Panzer_STK_MeshFactory_hpp__
#define Panzer_STK_MeshFactory_hpp__

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>

#include <stk_util/parallel/Parallel.hpp>

#include "Panzer_STK_PeriodicBC_Parser.hpp"

namespace panzer_stk {

class STK_Interface;

/** Pure virtual interface that constructs a 
  * STK_Mesh interface object.
  */
class STK_MeshFactory : public Teuchos::ParameterListAcceptorDefaultBase {
public:
   /** Construct a STK_Inteface object described
     * by this factory.
     *
     * \param[in] parallelMach Descriptor for machine to build this mesh on.
     *
     * \returns Pointer to <code>STK_Interface</code> object with 
     *          <code>isModifiable()==false</code>.
     */ 
   virtual Teuchos::RCP<STK_Interface> buildMesh(stk::ParallelMachine parallelMach) const = 0;

   /** This builds all the meta data of the mesh. Does not call metaData->commit.
     * Allows user to add solution fields and other pieces. The mesh can be "completed"
     * by calling <code>completeMeshConstruction</code>.
     */
   virtual Teuchos::RCP<STK_Interface> buildUncommitedMesh(stk::ParallelMachine parallelMach) const = 0;

   /** Finishes building a mesh object started by <code>buildUncommitedMesh</code>.
     */
   virtual void completeMeshConstruction(STK_Interface & mesh,stk::ParallelMachine parallelMach) const = 0;

   /** Parse the periodic boundary condition parameter list and build a vector of periodic boundary
     * conditions (a convenience function)
     */
   static void parsePeriodicBCList(const Teuchos::RCP<Teuchos::ParameterList> & pl,
                                   std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > & periodicBC)
   {
      panzer_stk::PeriodicBC_Parser parser;
      parser.setParameterList(pl);
      periodicBC = parser.getMatchers();
   }

protected:
   // vector of periodic boundary condition objects
   std::vector<Teuchos::RCP<const PeriodicBC_MatcherBase> > periodicBCVec_; 
};

}

#endif
