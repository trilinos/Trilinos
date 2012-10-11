#ifndef __Panzer_ResponseEvaluatorFactory_Functional_hpp__
#define __Panzer_ResponseEvaluatorFactory_Functional_hpp__

#include <string>

#include "Panzer_config.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ResponseEvaluatorFactory.hpp"

namespace panzer {

/** This class defines a response based on a functional.
  */
template <typename EvalT> 
class ResponseEvaluatorFactory_Functional : public ResponseEvaluatorFactory<EvalT> {
public:

   ResponseEvaluatorFactory_Functional() {}

   virtual ~ResponseEvaluatorFactory_Functional() {}
 
   /** Build the response object used by this factory. This object
     * assumes the role of the scatter target and will be accessible
     * by all the evaluators in the field managers. 
     *
     * \param[in] responseName Name of response to be built. This
     *                         name will be used for looking up
     *                         the response in the <code>GlobalEvaluationDataContainer</code>
     *                         object.
     */
   virtual Teuchos::RCP<ResponseBase> buildResponseObject(const std::string & responseName) const;
   
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
                                            const Teuchos::ParameterList & user_data) const;

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
                                            const Teuchos::ParameterList & user_data) const 
   { return Teuchos::null; }
};

}

#include "Panzer_ResponseEvaluatorFactory_Functional_impl.hpp"

#endif
