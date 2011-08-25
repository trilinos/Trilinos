#ifndef __Panzer_ResponseContainer_Derivative_hpp__
#define __Panzer_ResponseContainer_Derivative_hpp__

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

//! Derivative response type
template < >
class Response <panzer::Traits::Derivative> { };

// Derivative type responses
//////////////////////////////////////////

template < >
class ResponseContainer<panzer::Traits::Derivative> 
   : public ResponseContainerBase {
public:
   typedef panzer::Traits::Derivative RespType;

   ResponseContainer() {}

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

   static Teuchos::RCP<const Response<RespType> > 
   aggregateResponses(std::list<Teuchos::RCP<const Response<RespType> > > & responses)
   { return Teuchos::null; }
 
};

}

#endif
