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

#include <boost/unordered_map.hpp>
#include <boost/unordered_set.hpp>

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
                                                          const std::vector<panzer::WorksetDescriptor> & wkstDesc) const 
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

   /** Scale a field before writing out to STK by a prescribed value (the implicit default is 1.0). 
     * A warning will be printed if no field is found of that name when <code>buildAndRegisterEvaluators</code>
     * is called.
     *
     * \note Currently only works for HGRAD fields.
     *
     * \param[in] fieldName HGRAD field to scale.
     * \param[in] fieldScalar Value to scale the field by.
     */
   void scaleField(const std::string & fieldName,double fieldScalar);


private:
   void computeReferenceCentroid(const std::map<std::string,Teuchos::RCP<panzer::PureBasis> > & bases,
                                 int baseDimension,
                                 Intrepid::FieldContainer<double> & centroid) const;

   Teuchos::RCP<STK_Interface> mesh_;

   boost::unordered_map<std::string,double> fieldToScalar_;
   boost::unordered_set<std::string> scaledFieldsHash_; // used to print the warning about unused scaling
};

/** A simple builder for this the SolutionWriter response factory, simply set the mesh 
  * and this will build the response factories for you. (Pass into ResponseLibrary::addResponse)
  */
struct RespFactorySolnWriter_Builder {
  Teuchos::RCP<panzer_stk::STK_Interface> mesh;

  void scaleField(const std::string & fieldName,double fieldScalar)
  { fieldToScalar_[fieldName] = fieldScalar; }

  template <typename T>
  Teuchos::RCP<panzer::ResponseEvaluatorFactoryBase> build() const
  { 
    Teuchos::RCP<ResponseEvaluatorFactory_SolutionWriter<T> > ref = 
        Teuchos::rcp(new panzer_stk::ResponseEvaluatorFactory_SolutionWriter<T>(mesh)); 

    // set all scaled field values
    for(boost::unordered_map<std::string,double>::const_iterator itr=fieldToScalar_.begin();
        itr!=fieldToScalar_.end();++itr) 
      ref->scaleField(itr->first,itr->second);

    return ref;
  }

private:
  boost::unordered_map<std::string,double> fieldToScalar_;
};

}

#include "Panzer_STK_ResponseEvaluatorFactory_SolutionWriter_impl.hpp"

#endif
