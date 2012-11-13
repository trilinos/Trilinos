#ifndef __Panzer_STK_ResponseEvaluatorFactory_SolutionWriter_hpp__
#define __Panzer_STK_ResponseEvaluatorFactory_SolutionWriter_hpp__

#include <string>

#include "Panzer_config.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ResponseEvaluatorFactory.hpp"

#include "Panzer_STK_Interface.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

#include <mpi.h>

namespace panzer_stk {

/** This class defines a response based solution writer.
  */
template <typename EvalT> 
class ResponseEvaluatorFactory_SolutionWriter : public panzer::ResponseEvaluatorFactory<EvalT> {
public:

   ResponseEvaluatorFactory_SolutionWriter(const Teuchos::RCP<STK_Interface> & mesh)
     : mesh_(mesh) {}

   virtual ~ResponseEvaluatorFactory_SolutionWriter() {}
 
   /** Build the response object used by this factory. This object
     * assumes the role of the scatter target and will be accessible
     * by all the evaluators in the field managers. 
     *
     * \param[in] responseName Name of response to be built. This
     *                         name will be used for looking up
     *                         the response in the <code>GlobalEvaluationDataContainer</code>
     *                         object.
     */
   virtual Teuchos::RCP<panzer::ResponseBase> buildResponseObject(const std::string & responseName) const;

   virtual Teuchos::RCP<panzer::ResponseBase> buildResponseObject(const std::string & responseName,
                                                          const std::vector<std::string> & eBlocks) const 
   { return buildResponseObject(responseName); }

   virtual Teuchos::RCP<panzer::ResponseBase> buildResponseObject(const std::string & responseName,
                                                          const std::vector<std::pair<std::string,std::string> > & sideset_blocks) const
   { return buildResponseObject(responseName); }
   
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

   /** Take a vector of (std::string (field name), RCP<PureBasis>) pairs and bucket them
     * by basis name. What is returned is a map pairing the basis to a vector of field names.
     */
   static void bucketByBasisType(const std::vector<panzer::StrPureBasisPair> & providedDofs,
                                 std::map<std::string,std::vector<std::string> > & basisBucket);


private:
   void computeReferenceCentroid(const std::map<std::string,Teuchos::RCP<panzer::PureBasis> > & bases,
                                 int baseDimension,
                                 Intrepid::FieldContainer<double> & centroid) const;

   Teuchos::RCP<STK_Interface> mesh_;
};

/** A simple builder for this the SolutionWriter response factory, simply set the mesh and MPI_Comm
  * and this will build the response factories for you. (Pass into ResponseLibrary::addResponse)
  */
struct RespFactorySolnWriter_Builder {
  Teuchos::RCP<panzer_stk::STK_Interface> mesh;

  template <typename T>
  Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> build() const
  { return Teuchos::rcp(new panzer_stk::ResponseEvaluatorFactory_SolutionWriter<T>(mesh)); }
};

}

#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter_impl.hpp"

#endif
