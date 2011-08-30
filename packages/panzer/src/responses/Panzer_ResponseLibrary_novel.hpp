#ifndef __Panzer_ResponseLibrary_hpp__
#define __Panzer_ResponseLibrary_hpp__

#include <iostream>
#include <string>
#include <vector>
#include <map>

#include "Teuchos_ParameterList.hpp"

#include "Phalanx_FieldManager.hpp"

#include "Panzer_Response.hpp"
#include "Panzer_ResponseAggregatorBase.hpp"
#include "Panzer_ResponseFunctional_Aggregator.hpp"

namespace panzer {
namespace novel {

template <typename TraitsT> class ResponseContainerBase;
template <typename EvalT,typename TraitsT> class ResponseContainer;

/** This contains, collects and serves as a resource for
  * responses computed by panzer. This functions as a library
  * where there are many "responses" maintained (as many as
  * a user adds).  When a response is maintained that simply
  * means there is a mechansim to "reserve" it. A response is
  * not required by any field manager until a user "reserves"
  * it. The reservation process is done by response name and 
  * the element block or BC it is associated with. The field
  * tag specified in the <code>addResponse</code> is used only
  * for the "name" field, however that use of <code>PHX::FieldTag</code>
  * reminds the user that something better be in the evaluation tree.
  */
template <typename TraitsT>
class ResponseLibrary {
public:
   typedef typename TraitsT::EvalTypes TypeSeq;

   ResponseLibrary() {}

   /** Asks, does this string correspond to a response type
     * in this library?
     */
   template <typename EvalT>
   bool isResponseType(const std::string & type) const
   { return true; }

   //! Get a particular aggregator 
   template <typename EvalT>
   const ResponseAggregatorBase<TraitsT> & getAggregator(const std::string & type) const
   {
      static ResponseFunctional_Aggregator<EvalT,TraitsT> rfa;
      return rfa;
   }

   /** \defgroup volume
     * Volume methods for volumetric reponses
     * @{
     */

   /** User access to a volume response object. This will throw if called before
     * the volume field managers (all blocks) are evaluated.
     */
   Teuchos::RCP<const Response<TraitsT> > getVolumeResponse(const ResponseId & rid,
                                                            const std::string & eBlock) const;

   //! Reserve a response for actual calculation (by response id and element block).
   template <typename EvalT>
   void reserveVolumeResponse(const ResponseId & rid,const std::string & eBlock);

   /** Veryify that this response and element block are actual valid choices
     * for the evaluation type. This is optional error checking but makes debugging
     * simplier.
     */
   template <typename EvalT>
   bool validateResponseIdInElementBlock(const ResponseId & rid,const std::string & eBlock) const 
   { return true; }

   /** This method builds the volume field managers from the reserved
     * responses. It also registers a number of evaluators, using the closure
     * model and equation set factories. Unlike in the assembly engine only gather
     * evaluators, DOF, and Gradient evaluators are automaticlly included. This
     * method also stores the worksets passed in to be used when the evaluate method
     * is called.
     */
   void buildVolumeFieldManagersFromResponses(const Teuchos::ParameterList & pl);

   /** Evaluate all the volume field managers of a particular evaluator type.
     */
   template <typename EvalT>
   void evaluateVolumeFieldManagers();

   /** @} */

   /** Returns the set of element blocks required by the ResponseLibrary to 
     * build any responses. This only details those element blocks the library is
     * currently aware of, if new responses are reserved in new element blocks
     * those will not be included. Additionally, this defines the set of worksets
     * required for the evaluating the responses.
     */
   void getRequiredElementBlocks(std::vector<std::string> & eBlocks) const;

protected:
   //! Access a container field for a specified element block
   template <typename EvalT>
   Teuchos::RCP<ResponseContainerBase<TraitsT> > getVolumeContainer(const std::string & eBlock);

private:
   // This could be a template manager, but this turns out to be more in line with
   // what is desired here. Direct access to the vector.
   typedef std::vector<Teuchos::RCP<ResponseContainerBase<TraitsT> > > RespContVector;
   typedef std::vector<Teuchos::RCP<PHX::FieldManager<TraitsT> > > FMVector;

   std::map<std::string,Teuchos::RCP<RespContVector> > rsvdVolResp_;
   std::map<std::string,Teuchos::RCP<FMVector> > volFieldManagers_;
};

}
}

#include "Panzer_ResponseLibrary_novelT.hpp"

#endif
