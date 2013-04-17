#ifndef __Panzer_ResponseEvaluatorFactory_Functional_hpp__
#define __Panzer_ResponseEvaluatorFactory_Functional_hpp__

#include <string>

#include "Panzer_config.hpp"
#include "Panzer_PhysicsBlock.hpp"
#include "Panzer_BC.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ResponseEvaluatorFactory.hpp"
#include "Panzer_LinearObjFactory.hpp"

#include <mpi.h>

namespace panzer {

/** This class defines a response based on a functional.
  */
template <typename EvalT,typename LO,typename GO> 
class ResponseEvaluatorFactory_Functional : public ResponseEvaluatorFactory<EvalT> {
public:

   ResponseEvaluatorFactory_Functional(MPI_Comm comm, int cubatureDegree=1,bool requiresCellIntegral=true,const std::string & quadPointField="",
                                       const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & linearObjFactory=Teuchos::null,
                                       const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & globalIndexer=Teuchos::null)
     : comm_(comm), cubatureDegree_(cubatureDegree), requiresCellIntegral_(requiresCellIntegral)
     , quadPointField_(quadPointField), linearObjFactory_(linearObjFactory), globalIndexer_(globalIndexer)
   {
     TEUCHOS_ASSERT((linearObjFactory==Teuchos::null && globalIndexer==Teuchos::null) ||
                    (linearObjFactory!=Teuchos::null && globalIndexer!=Teuchos::null));
   }

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

   virtual Teuchos::RCP<ResponseBase> buildResponseObject(const std::string & responseName,
                                                          const std::vector<WorksetDescriptor> & wkstDesc) const 
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

   /** Is this evaluation type supported by the factory. This is used to determine cases
     * where a response may support a particular evaluation type, however at runtime the user
     * decides not to enable the (say) Jacobian evaluation of this response.
     *
     * Note that use of this mechanism is complementary to having the builder return 
     * <code>Teuchos::null</code> for a particular evaluation type.
     */
   virtual bool typeSupported() const;

protected:
   //! Accessor method for Cubature degree (can be used by sub classes)
   int getCubatureDegree() const { return cubatureDegree_; }

private:
   MPI_Comm comm_;
   int cubatureDegree_;
   bool requiresCellIntegral_;
   std::string quadPointField_;
   Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linearObjFactory_;
   Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > globalIndexer_;
};

template <typename LO,typename GO> 
struct FunctionalResponse_Builder {
  MPI_Comm comm;
  int cubatureDegree;
  bool requiresCellIntegral;
  std::string quadPointField;
  void setDerivativeInformation(const Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > & in_linearObjFactory,
                                const Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > & in_globalIndexer)
  {
    linearObjFactory = in_linearObjFactory;
    globalIndexer = in_globalIndexer;

    TEUCHOS_ASSERT((linearObjFactory==Teuchos::null && globalIndexer==Teuchos::null) ||
                   (linearObjFactory!=Teuchos::null && globalIndexer!=Teuchos::null));
  }

  template <typename T>
  Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> build() const
  { return Teuchos::rcp(new ResponseEvaluatorFactory_Functional<T,LO,GO>(comm,cubatureDegree,requiresCellIntegral,quadPointField,linearObjFactory,globalIndexer)); }

private:
  Teuchos::RCP<panzer::LinearObjFactory<panzer::Traits> > linearObjFactory;
  Teuchos::RCP<const panzer::UniqueGlobalIndexer<LO,GO> > globalIndexer;
};


}

#include "Panzer_ResponseEvaluatorFactory_Functional_impl.hpp"

#endif
