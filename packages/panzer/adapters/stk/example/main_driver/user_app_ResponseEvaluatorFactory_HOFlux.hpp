#ifndef __user_app_ResponseEvaluatorFactory_HOFlux_hpp__
#define __user_app_ResponseEvaluatorFactory_HOFlux_hpp__

#include <string>

#include "Panzer_config.hpp"
#include "Panzer_ResponseEvaluatorFactory_Functional.hpp"

namespace user_app {

/** This class defines a response based on a functional. */
  
template <typename EvalT> 
class ResponseEvaluatorFactory_HOFlux : public panzer::ResponseEvaluatorFactory_Functional<EvalT> 
{
public:

   ResponseEvaluatorFactory_HOFlux(MPI_Comm comm, int cubatureDegree)
     : panzer::ResponseEvaluatorFactory_Functional<EvalT>(comm,cubatureDegree,false)
   {}

   virtual ~ResponseEvaluatorFactory_HOFlux() {}
   
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
                                           const Teuchos::ParameterList & user_data) const;

};

struct HOFluxResponse_Builder {
  MPI_Comm comm;
  int cubatureDegree;

  template <typename T>
  Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> build() const
  { return Teuchos::rcp(new ResponseEvaluatorFactory_HOFlux<T>(comm,cubatureDegree)); }
};

}

#include "user_app_ResponseEvaluatorFactory_HOFlux_impl.hpp"

#endif
