#ifndef Panzer_STK_MeshFactory_hpp__
#define Panzer_STK_MeshFactory_hpp__

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterListAcceptorDefaultBase.hpp>

#include <stk_util/parallel/Parallel.hpp>

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
};

}

#endif
