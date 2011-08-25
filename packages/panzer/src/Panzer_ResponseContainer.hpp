#ifndef __Panzer_ResponseContainer_hpp__
#define __Panzer_ResponseContainer_hpp__

#include <vector>
#include <string>
#include <list>

#include "Panzer_Traits.hpp"

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

   //! Does this container have a response of a specified name
   virtual bool contains(const std::string & name) const
   {
      std::vector<std::string>::const_iterator itr = std::find(names_.begin(),names_.end(),name);
      return itr!=names_.end();
   }

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

// Forward declaration
template <typename RespType> class Response;

}

#endif
