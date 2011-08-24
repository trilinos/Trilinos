#ifndef __Panzer_ResponseContainer_hpp__
#define __Panzer_ResponseContainer_hpp__

#include <vector>
#include <string>

#include "Panzer_Traits.hpp"

#include "Panzer_Response.hpp"

#include "Teuchos_Comm.hpp"
#include "Teuchos_RCP.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_Evaluator_WithBaseImpl.hpp"

namespace panzer {

/** Base class for use with the Phalanx TemplateManager.
  * This object will serve as an accessor to the various
  * response objects.
  */
class ResponseContainerBase {
public:
   virtual ~ResponseContainerBase() {}

   //! What is the evaluation type for this response
   virtual const std::string & getResponseType() const = 0;

   //! Reserve a particular response
   virtual void reserve(const std::string & name)
   { names_.push_back(name); }

   //! get a list of all reserved responses
   virtual void getReserved(std::vector<std::string> & names) const
   { names = names_; }

   std::size_t getFieldIndex(const std::string & name) const
   {
      std::vector<std::string>::const_iterator itr = std::find(names_.begin(),names_.end(),name);
      TEST_FOR_EXCEPTION(itr==names_.end(),std::logic_error,
                         "Could not find response field \""+name+"\" in response container");
      return itr-names_.begin();
   }

   //! Register responses with a particular field manager
   virtual void registerResponses(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,int worksetSize,
                          PHX::FieldManager<panzer::Traits> & fm) = 0;

   //! Clear response containers data
   virtual void clear() = 0;

   virtual void addDependentFields(PHX::EvaluatorWithBaseImpl<panzer::Traits> & eval) const = 0;
   virtual void setFieldData(PHX::FieldManager<panzer::Traits> & eval) = 0;
   virtual void aggregateFieldsLocally(const panzer::Workset & wkst) = 0;
   virtual void reduceFieldsGlobally() = 0;
  
private:

   std::vector<std::string> names_;
};

// forward declaration
template <typename RespT> class ResponseContainer;

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

private:
   //! build response fields, called by registerResponses
   void buildResponseFields(int worksetSize);
   
   bool responseFieldsBuilt_;
   std::vector<PHX::MDField<RespType::ScalarT,Cell> > responseFields_;
   std::vector<Sacado::ScalarType<RespType::ScalarT>::type> responseVector_; 
   std::vector<Sacado::ScalarType<RespType::ScalarT>::type> localResponseVector_; 

   Teuchos::RCP<const Teuchos::Comm<int> > comm_;
};

// Derivative type responses
//////////////////////////////////////////

template < >
class ResponseContainer<panzer::Traits::Derivative> 
   : public ResponseContainerBase {
public:
   typedef panzer::Traits::Derivative RespType;

   ResponseContainer();

   virtual const std::string & getResponseType() const
   { return PHX::TypeString<RespType>::value; }

   //! Register responses with a particular field manager
   void registerResponses(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,int worksetSize,
                          PHX::FieldManager<panzer::Traits> & fm) {}

   void clear() { }

   void addDependentFields(PHX::EvaluatorWithBaseImpl<panzer::Traits> & eval) const {}
   void setFieldData(PHX::FieldManager<panzer::Traits> & eval) {}
   void aggregateFieldsLocally(const panzer::Workset & wkst) {}
   void reduceFieldsGlobally() {}

   //! This returns the response for a given field
   Teuchos::RCP<const Response<RespType> > getResponse(const std::string & name) const
   { return Teuchos::null; }
};

}

#endif
