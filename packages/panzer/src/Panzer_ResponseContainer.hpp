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

   /** \defgroup Methods common to all ResponseContainers
     * These methods are implementations of a basic set of functionality
     * common to all ResponseContainer classes.
     * @{
     */

   //! Reserve a particular response
   void reserve(const std::string & name)
   { names_.push_back(name); }

   //! Does this container have a response of a specified name
   bool contains(const std::string & name) const
   {
      std::vector<std::string>::const_iterator itr = std::find(names_.begin(),names_.end(),name);
      return itr!=names_.end();
   }

   //! get a list of all reserved responses
   void getReserved(std::vector<std::string> & names) const
   { names = names_; }

   std::size_t getFieldIndex(const std::string & name) const
   {
      std::vector<std::string>::const_iterator itr = std::find(names_.begin(),names_.end(),name);
      TEST_FOR_EXCEPTION(itr==names_.end(),std::logic_error,
                         "Could not find response field \""+name+"\" in response container");
      return itr-names_.begin();
   }

   /** @} */

   /** \defgroup Pure virtual funcions for setting up field managers
     * @{
     */

   //! What is the evaluation type for this response
   virtual const std::string & getResponseType() const = 0;

   /** Register responses with a particular field manager.
     * This essentially builds a ScatterResponse object 
     * and registers it as a required field in the field manager.
     */
   virtual void registerResponses(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,int worksetSize,
                          PHX::FieldManager<panzer::Traits> & fm) = 0;

   /** Register dependent fields (responses) defined by ResponseContainer with an evaluator.
     * This should register the reserved fields.
     */
   virtual void addDependentFields(PHX::EvaluatorWithBaseImpl<panzer::Traits> & eval) const = 0;

   /** @} */

   /** \defgroup Pure viritual for modification and aggregation of data.
     * @{
     */

   //! Clear response containers data
   virtual void clear() = 0;

   //! Load up dependent field data into internally stored fields
   virtual void setFieldData(PHX::FieldManager<panzer::Traits> & eval) = 0;

   //! Aggregate locally stored fields, for this workset.
   virtual void aggregateFieldsLocally(const panzer::Workset & wkst) = 0;

   //! Aggregate over other (similar) responses on other processors.
   virtual void reduceFieldsGlobally() = 0;

   /** @} */

   //
   // This is a static member sub classes must implement. This takes a list
   // of responses and does an aggregation step (max,min,sum, etc...) over
   // many Response objects
   //
   // static Teuchos::RCP<const Response<RespType> > 
   // aggregateResponses(std::list<Teuchos::RCP<const Response<RespType> > > & responses);
   //
  
private:

   std::vector<std::string> names_;
};

// forward declaration
template <typename RespT> class ResponseContainer;

// Forward declaration
template <typename RespType> class Response;

}

#endif
