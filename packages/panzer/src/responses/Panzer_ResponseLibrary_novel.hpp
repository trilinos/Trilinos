#ifndef __Panzer_ResponseLibrary_hpp__
#define __Panzer_ResponseLibrary_hpp__

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>

#include "Teuchos_ParameterList.hpp"

#include "Phalanx_FieldManager.hpp"
#include "Phalanx_TemplateManager.hpp"

#include "Panzer_ResponseAggregatorBase.hpp"
#include "Panzer_ResponseFunctional_Aggregator.hpp"

namespace panzer {
namespace novel {

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

   ResponseLibrary() {}

   /** Asks, does this string correspond to a response type
     * in this library?
     */
   template <typename EvalT>
   bool isResponseType(const std::string & type) const
   { return true; }

   template <typename EvalT>
   const ResponseAggregatorBase<TraitsT> & getAggregator(const std::string & type) const
   {
      static ResponseFunctional_Aggregator<EvalT,TraitsT> rfa;
      return rfa;
   }
#if 0
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

   /** Register reserved responses with a field manager. This only registers those
     * responses in a user specified element block. The workset size is the maximum
     * number of cells in any workset.
     *
     * \param[in] eBlock Registered fields will belong to the element block with this ID.
     * \param[in] comm Parallel communicator (needed for parallel reduce type operations)
     * \param[in] worksetSize Maximum number of cells in any given workset to be evaluated
     *                        by the field manager
     * \param[in] fm Field manager to register responses with.
     */
   template <typename RespT>
   void registerReservedResponses(const std::string & eBlock,
                                  const Teuchos::RCP<const Teuchos::Comm<int> > & comm,
                                  int worksetSize,
                                  PHX::FieldManager<panzer::Traits> & fm);


   //! Register volume responses with the specified field manager. Template specifies response type.
   /** Get volume response associated with the response type (template parameter) and
     * the field name and element block.
     *
     * \param[in] name Field name
     * \param[in] eBlock Element block id
     *
     * \return Response object associated with the field on the specified element block.
     *
     * \note If no element block and field pair is found then a std::logic_error exception is thrown.
     */
   template <typename RespT>
   Teuchos::RCP<const Response<RespT> > getVolumeResponse(const std::string & name,
                                                          const std::string & eBlock) const;

   /** Get volume response associated with the response type (template parameter) and
     * the field name. The response is "aggregated" (dependeing on field type) over all
     * element blocks.
     *
     * \param[in] name Field name
     *
     * \return Response object associated with the field on the specified element block.
     *
     * \note If no element block and field pair is found then a std::logic_error exception is thrown.
     */
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
#endif
};

}
}

#endif
