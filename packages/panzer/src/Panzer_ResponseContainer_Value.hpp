#ifndef __Panzer_ResponseContainer_Value_hpp__
#define __Panzer_ResponseContainer_Value_hpp__

#include <vector>
#include <string>
#include <list>

#include "Panzer_Traits.hpp"

#include "Panzer_ResponseContainer.hpp"

#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"

namespace panzer {

//! Value type response
template < >
class Response <panzer::Traits::Value> {
public:
   typedef panzer::Traits::Value RespType;

   Response(Sacado::ScalarType<RespType::ScalarT>::type value)
      : value_(value) {}

   //! get value of this response
   Sacado::ScalarType<RespType::ScalarT>::type getValue() const
   { return value_; }

private:
   // hide these constructors
   Response();
   Response(const Response<panzer::Traits::Value> &);
    
   Sacado::ScalarType<RespType::ScalarT>::type value_;
};

// Value type responses
//////////////////////////////////////////

/** A response container that handles summed values.
  */
template < >
class ResponseContainer<panzer::Traits::Value> 
   : public ResponseContainerBase {
public:
   typedef panzer::Traits::Value RespType;

   ResponseContainer();

   virtual const std::string & getResponseType() const
   { return PHX::TypeString<RespType>::value; }

   void registerResponses(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,int worksetSize,
                          PHX::FieldManager<panzer::Traits> & fm);

   void clear() 
   { for(std::size_t i=0;i<responseVector_.size();i++) responseVector_[i] = 0.0; }

   void addDependentFields(PHX::EvaluatorWithBaseImpl<panzer::Traits> & eval) const;
   void setFieldData(PHX::FieldManager<panzer::Traits> & eval);
   void aggregateFieldsLocally(const panzer::Workset & wkst);
   void reduceFieldsGlobally();

   //! This returns the response for a given field
   Teuchos::RCP<const Response<RespType> > getResponse(const std::string & name) const;

   static Teuchos::RCP<const Response<RespType> > 
   aggregateResponses(std::list<Teuchos::RCP<const Response<RespType> > > & responses);

private:
   //! build response fields, called by registerResponses
   void buildResponseFields(int worksetSize);
   
   bool responseFieldsBuilt_;
   std::vector<PHX::MDField<RespType::ScalarT,Cell> > responseFields_;
   std::vector<Sacado::ScalarType<RespType::ScalarT>::type> responseVector_; 
   std::vector<Sacado::ScalarType<RespType::ScalarT>::type> localResponseVector_; 

   Teuchos::RCP<const Teuchos::Comm<int> > comm_;
};

}

#endif
