#ifndef __Panzer_ResponseEvaluatorFactoryBase_hpp__
#define __Panzer_ResponseEvaluatorFactoryBase_hpp__

#include <string>

#include "Teuchos_RCP.hpp"

#include "Panzer_config.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ResponseBase.hpp"
#include "Panzer_WorksetDescriptor.hpp"

#include "Phalanx_FieldTag.hpp"

namespace panzer {

/** This pure virtual class defines the methods to be used
  * to construct a response using field managers. It is assumed
  * that multiple field managers can be populated with this
  * response, thus the object itself should remain stateless.
  * 
  * \note Users should not derive directly off of this object,
  *       but should instead derive from the ResponseEvaluatorFactory
  *       which is the templated version with a specific evaluation
  *       type.
  */
class ResponseEvaluatorFactoryBase {
public:

   ResponseEvaluatorFactoryBase() {}

   virtual ~ResponseEvaluatorFactoryBase() {}

   /** Build the response object used by this factory. This object
     * assumes the role of the scatter target and will be accessible
     * by all the evaluators in the field managers. This is the sideset
     * version of the buildResponseObject function.
     *
     * \param[in] responseName Name of response to be built. This
     *                         name will be used for looking up
     *                         the response in the <code>GlobalEvaluationDataContainer</code>
     *                         object.
     * \param[in] wkstdescs A vector of descriptors for the elements this response is over.
     */
   virtual Teuchos::RCP<ResponseBase> buildResponseObject(const std::string & responseName,
                                                          const std::vector<WorksetDescriptor> & wkstdescs) const = 0; 
   
   /** Build and register evaluators for a response on a particular physics
     * block.
     *
     * \param[in] responseName The name of the response to be constructed
     *                         by these evaluators.
     * \param[in,out] fm Field manager to be fuild with the evaluators.
     * \param[in] physicsBlock What physics block is being used for constructing
     *                         the evaluators
     * \param[in] user_data The user data parameter list, this stores things
     *                      that the user may find useful.
     */
   virtual void buildAndRegisterEvaluators(const std::string & responseName,
                                           PHX::FieldManager<panzer::Traits> & fm,
                                           const panzer::PhysicsBlock & physicsBlock,
                                           const Teuchos::ParameterList & user_data) const = 0;

   /** Is this evaluation type supported by the factory. This is used to determine cases
     * where a response may support a particular evaluation type, however at runtime the user
     * decides not to enable the (say) Jacobian evaluation of this response.
     *
     * Note that use of this mechanism is complementary to having the builder return 
     * <code>Teuchos::null</code> for a particular evaluation type.
     *
     * \note In this contect the "type" does not make a lot of sense. But in the dervied
     *       interface <code>ResponseEvaluatorFactory<EvalT></code> the type is the <code>EvalT</code>.
     *       Inclusion of this method here simply makes dynamic access to this method
     *       possible with out a cast. In the end it cleans up the code.
     */
   virtual bool typeSupported() const = 0;
};

}

#endif
