#ifndef __Panzer_ResponseEvaluatorFactoryBase_hpp__
#define __Panzer_ResponseEvaluatorFactoryBase_hpp__

#include <string>

#include "Teuchos_RCP.hpp"

#include "Panzer_config.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ResponseBase.hpp"

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
     * by all the evaluators in the field managers. 
     *
     * \param[in] responseName Name of response to be built. This
     *                         name will be used for looking up
     *                         the response in the <code>GlobalEvaluationDataContainer</code>
     *                         object.
     */
   virtual Teuchos::RCP<ResponseBase> buildResponseObject(const std::string & responseName) const = 0; 
   
   /** Build and register evaluators for a response on a particular physics
     * block. Note the required returned field tag is needed for specifying
     * a scatter evaluator. This does not prevent a developer from marking
     * something as required within the buildAndRegisterEvaluators method.
     *
     * \param[in] responseName The name of the response to be constructed
     *                         by these evaluators.
     * \param[in,out] fm Field manager to be fuild with the evaluators.
     * \param[in] physicsBlock What physics block is being used for constructing
     *                         the evaluators
     * \param[in] user_data The user data parameter list, this stores things
     *                      that the user may find useful.
     *
     * \returns Field tag that corresponds to the required scatter field for
     *          this response. 
     */
   virtual Teuchos::RCP<const PHX::FieldTag> buildAndRegisterEvaluators(const std::string & responseName,
                                            PHX::FieldManager<panzer::Traits> & fm,
                                            const panzer::PhysicsBlock & physicsBlock,
                                            const Teuchos::ParameterList & user_data) const = 0;

   /** Build and register evaluators for a response on a particular side set.
     * Note the required returned field tag is needed for specifying
     * a scatter evaluator. This does not prevent a developer from marking
     * something as required within the buildAndRegisterEvaluators method.
     *
     * \param[in] responseName The name of the response to be constructed
     *                         by these evaluators.
     * \param[in,out] fm Field manager to be fuild with the evaluators.
     * \param[in] bc Boundary condition object for the sideset being used.
     * \param[in] physicsBlock What physics block is being used for constructing
     *                         the evaluators
     * \param[in] user_data The user data parameter list, this stores things
     *                      that the user may find useful.
     *
     * \returns Field tag that corresponds to the required scatter field for
     *          this response. 
     */
   virtual Teuchos::RCP<const PHX::FieldTag> buildAndRegisterEvaluators(const std::string & responseName,
                                            PHX::FieldManager<panzer::Traits> & fm,
                                            const panzer::BC & bc,
                                            const panzer::PhysicsBlock & physicsBlock,
                                            const Teuchos::ParameterList & user_data) const = 0;
};

}

#endif
