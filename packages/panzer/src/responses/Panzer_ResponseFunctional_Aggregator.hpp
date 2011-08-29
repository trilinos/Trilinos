#ifndef __Panzer_ResponseFunctional_Aggregator_hpp__
#define __Panzer_ResponseFunctional_Aggregator_hpp__

#include "Panzer_config.hpp"

#include <string>
#include <vector>

#include "Panzer_ResponseAggregatorBase.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace panzer {

/** Abstract data object for response storage. Also manages
  * setting up and transfering to the responses. This works
  * in concert with the response aggregator.
  */
template <typename EvalT,typename TraitsT> class ResponseFunctional_Data; 
template <typename EvalT,typename TraitsT> class ResponseFunctional_Aggregator; 

template <typename TraitsT> 
class ResponseFunctional_Data<panzer::Traits::Residual,TraitsT> 
   : public ResponseData<TraitsT> {
public:
   ResponseFunctional_Data() {}

   /** \defgroup Required from ResponseData
     * @{
     */

   //! Allocate and initialize required storage based on the fields passed in
   virtual void allocateAndInitializeData(const std::vector<std::string> & fields)
   { 
      fields_ = fields;
      data_.resize(fields.size());
      reinitializeData();
   }

   /** Reinitialize the data based on the fields original requested by
     * <code>allocateAndInitializeData</code>.
     */
   virtual void reinitializeData()
   { for(std::size_t i=0;i<data_.size();i++) data_[i] = 0.0; }
   

   /** Fill a response object with data from a particular field.
     */
   virtual void fillResponse(const std::string & field,Response<TraitsT> & response) const
   {
      std::size_t index = fieldIndex(field);
      TEST_FOR_EXCEPTION(index>=fields_.size(),std::logic_error,
                         "Cannnot find field \""+field+"\" in ReponseFunctional_Data object.");
      response.setValue(data_[index]);
   }

   //! Is this field contain in the response data
   virtual bool contains(const std::string & field) const
   { return this->fieldIndex(field)<fields_.size(); }

   //! @}
 
   //! Get field index
   std::size_t fieldIndex(const std::string & field) const
   {
      std::vector<std::string>::const_iterator itr 
            = std::find(fields_.begin(),fields_.end(),field);
      return itr-fields_.begin();
   }

   //! Get vector of strings describing owned fields
   const std::vector<std::string> & getFields() const
   { return fields_; }

   //! get the vector with data (writable)
   std::vector<typename TraitsT::RealType> & getData() 
   { return data_; }

   //! get the vector with data (non-writable)
   const std::vector<typename TraitsT::RealType> & getData() const
   { return data_; }

private:
   //! what fields this stores
   std::vector<std::string> fields_;

   //! values stored
   std::vector<typename TraitsT::RealType> data_;
};

template <typename TraitsT>
class ResponseFunctional_Aggregator<panzer::Traits::Residual,TraitsT>
   : public ResponseAggregatorBase<TraitsT> {
public:
   // useful for cloning and the factory mechanism
   ResponseFunctional_Aggregator();

   ResponseFunctional_Aggregator(const Teuchos::ParameterList & p);

   /** \defgroup Methods required by ResponseAggregatorBase
     * @{
     */

   //! Clone this aggregator
   Teuchos::RCP<ResponseAggregatorBase<TraitsT> > clone(const Teuchos::ParameterList & p) const;

   //! Build response data for a specified set of fields (ResponseAggregator decides data layout)
   Teuchos::RCP<ResponseData<TraitsT> > buildResponseData(const std::vector<std::string> & fields) const;

   //! Register and build evaluator required by this aggregator and this set of data.
   virtual void registerAndRequireEvaluators(PHX::FieldManager<TraitsT> & fm,Teuchos::RCP<ResponseData<TraitsT> > & data, 
                                             const Teuchos::ParameterList & p) const;

   //! perform global reduction on this set of response data
   void globalReduction(const Teuchos::Comm<int> & comm,ResponseData<TraitsT>  & rd) const;

   //! @}

   //! Build an evaluator of the set of fields to be aggregated (calculated) together
   Teuchos::RCP<PHX::Evaluator<TraitsT> > buildEvaluator(Teuchos::RCP<ResponseData<TraitsT> > & data,
                                                         const Teuchos::ParameterList & p) const;

   //! Aggregate fields into a specific data object
   void evaluateFields(panzer::Workset & wkst,ResponseFunctional_Data<panzer::Traits::Residual,TraitsT> & data,
                                              const std::vector<PHX::MDField<panzer::Traits::Residual::ScalarT,Cell> > & fields) const;

};

}

#include "Panzer_ResponseFunctional_AggregatorT.hpp"

#endif
