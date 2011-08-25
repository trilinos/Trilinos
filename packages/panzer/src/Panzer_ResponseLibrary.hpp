#ifndef __Panzer_ResponseLibrary_hpp__
#define __Panzer_ResponseLibrary_hpp__

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>

#include "Teuchos_ParameterList.hpp"

#include "Phalanx_FieldTag.hpp"
#include "Phalanx_FieldManager.hpp"
#include "Phalanx_TemplateManager.hpp"

#include "Panzer_BC.hpp"
#include "Panzer_Traits.hpp"
#include "Panzer_ResponseContainer.hpp"
#include "Panzer_ResponseContainer_Value.hpp"
#include "Panzer_ResponseContainer_Derivative.hpp"

namespace panzer {

/** This contains, collects and serves as a resource for
  * responses computed by panzer. This functions as a library
  * where there are many "responses" maintained (as many as
  * a user adds).  When a response is maintained that simply
  * means there is a mechansim to "check it out". A response is
  * not required by any field manager until a user "checks it
  * out". The checkout process is done by response name and 
  * can be specified by a parameter list. It is assume that
  * each response parameter is defined on a 
  * <code>MDALayout<ScalarT,Cell></code> data type. The field
  * tag specified in the <code>addResponse</code> is used only
  * for the "name" field, however that use of <code>PHX::FieldTag</code>
  * reminds the user that something better be in the evaluation tree.
  *
  * The last aspect of the response library is that when a user
  * wants to access the response (typically a model evaluator) they
  * must specify the mechanism to use. TBD
  */
class ResponseLibrary {
public:
   ResponseLibrary();

   //////////////////////////////////////////////////////////////////////
   /** \defgroup BC and Volume methods
     * Methods that act on both boundary condition and volume responses.
     * @{
     */

   //! Print available reponses
   void printAvailableResponses(std::ostream & os) const;

   //! Print available reponses
   void printReservedResponses(std::ostream & os) const;

   //! Print available and checked out responses
   void print(std::ostream & os) const;

   /** @} */

   //////////////////////////////////////////////////////////////////////
   /** \defgroup Volume methods
     * Methods that act on volume responses.
     * @{
     */

   //! Include a volume type response in the library
   void addResponse(const PHX::FieldTag & ft, const std::string & blockId);

   //! Get the names of the available volume responses 
   void getAvailableVolumeResponses(std::vector<std::pair<std::string,std::set<std::string> > > & responsePairs) const;

    //! Reserve a response for actual calculation by name and element block.
   template <typename RespT>
   void reserveVolumeResponse(const std::string & name,const std::string & eBlock);

   template <typename RespT>
   void registerReservedResponses(const std::string & eBlock,
                                  const Teuchos::RCP<const Teuchos::Comm<int> > & comm,int worksetSize,
                                  PHX::FieldManager<panzer::Traits> & fm);

   template <typename RespT>
   Teuchos::RCP<const Response<RespT> > getVolumeResponse(const std::string & name,
                                                          const std::string & eBlock) const;

   template <typename RespT>
   Teuchos::RCP<const Response<RespT> > getVolumeResponse(const std::string & name) const;

   /** @} */

private:
   typedef std::map<std::string,std::set<std::string> > VolumeMap;
   typedef PHX::TemplateManager<panzer::Traits::RespTypes,ResponseContainerBase,ResponseContainer<_> > RespContManager;

   void initializeReservedVolumes();

   //! Contains a library of responses and the element blocks they correspond to: field name -> block ids
   VolumeMap volumeResponses_;

   std::map<std::string,Teuchos::RCP<RespContManager> > rsvdVolResp_;
   bool isInitialized_;
};

inline std::ostream & operator<<(std::ostream & os,const ResponseLibrary & rl)
{ rl.print(os); return os; }

template <typename RespT>
void ResponseLibrary::
reserveVolumeResponse(const std::string & name,const std::string & eBlock)
{
   using Teuchos::RCP;

   // look up in available responses, and build reserved volume reponses map
   initializeReservedVolumes();

   std::map<std::string,Teuchos::RCP<RespContManager> >::iterator itr = rsvdVolResp_.find(eBlock);
   TEST_FOR_EXCEPTION(itr==rsvdVolResp_.end(),std::logic_error,"Could not find element block \""+eBlock+"\"");

   Teuchos::RCP<ResponseContainer<RespT> > respContainer = itr->second->getAsObject<RespT>();
   TEUCHOS_ASSERT(respContainer!=Teuchos::null);
   respContainer->reserve(name);
}

template <typename RespT>
void ResponseLibrary::
registerReservedResponses(const std::string & eBlock,
                          const Teuchos::RCP<const Teuchos::Comm<int> > & comm,int worksetSize,
                          PHX::FieldManager<panzer::Traits> & fm)
{
   std::map<std::string,Teuchos::RCP<RespContManager> >::iterator itr = rsvdVolResp_.find(eBlock);
   TEST_FOR_EXCEPTION(itr==rsvdVolResp_.end(),std::logic_error,"Could not find element block \""+eBlock+"\"");

   Teuchos::RCP<ResponseContainer<RespT> > respContainer = itr->second->getAsObject<RespT>();
   TEUCHOS_ASSERT(respContainer!=Teuchos::null);
   
   respContainer->registerResponses(comm,worksetSize,fm);
}

template <typename RespT>
Teuchos::RCP<const Response<RespT> > ResponseLibrary::
getVolumeResponse(const std::string & name, const std::string & eBlock) const
{
   std::map<std::string,Teuchos::RCP<RespContManager> >::const_iterator itr = rsvdVolResp_.find(eBlock);
   TEST_FOR_EXCEPTION(itr==rsvdVolResp_.end(),std::logic_error,
                      "panzer::ResponseLibrary::getVolumeResponse could not find element block \""+eBlock+"\"");

   Teuchos::RCP<const ResponseContainer<RespT> > respContainer = itr->second->getAsObject<RespT>();
   TEUCHOS_ASSERT(respContainer!=Teuchos::null);

   // if response container contains desired response, return
   if(respContainer->contains(name))
      return respContainer->getResponse(name);

   TEST_FOR_EXCEPTION(true,std::logic_error,
                      "panzer::ResponseLibrary::getVolumeResponse could not find reserved field \""+name+
                      "\" in element block \""+eBlock+"\"");

   return Teuchos::null; // never execueted!
}

template <typename RespT>
Teuchos::RCP<const Response<RespT> > ResponseLibrary::
getVolumeResponse(const std::string & name) const
{
   using Teuchos::RCP;
 
   std::list<RCP<const Response<RespT> > > responses;

   // loop over all element blocks
   for(std::map<std::string,Teuchos::RCP<RespContManager> >::const_iterator itr=rsvdVolResp_.begin();
       itr!=rsvdVolResp_.end();++itr) {
      // extract response container of specified type
      Teuchos::RCP<const ResponseContainer<RespT> > respContainer = itr->second->getAsObject<RespT>();
      TEUCHOS_ASSERT(respContainer!=Teuchos::null); // sanity check

      // if response container contains desired response, add to responses list
      if(respContainer->contains(name))
         responses.push_back(respContainer->getResponse(name));
   }

   if(responses.size()>0)
      return ResponseContainer<RespT>::aggregateResponses(responses);

   TEST_FOR_EXCEPTION(true,std::logic_error,
                      "panzer::ResponseLibrary::getVolumeResponse could not find reserved field \""+name+"\"");

   return Teuchos::null;
}

}

#endif
