#ifndef __Panzer_ResponseEvaluatorFactory_hpp__
#define __Panzer_ResponseEvaluatorFactory_hpp__

#include <string>

#include "Panzer_config.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ResponseEvaluatorFactoryBase.hpp"

namespace panzer {

/** This pure virtual class defines the methods to be used
  * to construct a response using field managers. It is assumed
  * that multiple field managers can be populated with this
  * response, thus the object itself should remain stateless.
  *
  * \note Users should derive off of this class directly. The
  *       <code>ReponseEvaluatorFactory_TemplateManager</code> 
  *       should be used to access these objects.
  */
template <typename EvalT> 
class ResponseEvaluatorFactory : public ResponseEvaluatorFactoryBase {
public:

   ResponseEvaluatorFactory() {}

   virtual ~ResponseEvaluatorFactory() {}
   
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
     * block. Note that it is assumed that a field has been marked required
     * during this method call.
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
     */
   virtual bool typeSupported() const = 0;
};

}

#endif
