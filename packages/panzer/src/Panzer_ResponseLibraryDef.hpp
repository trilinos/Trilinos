#ifndef __Panzer_ResponseLibraryDef_hpp__
#define __Panzer_ResponseLibraryDef_hpp__

#include <iostream>
#include <string>
#include <vector>
#include <set>
#include <map>

#include "Teuchos_ParameterList.hpp"

#include "Phalanx_FieldTag.hpp"
#include "Phalanx_FieldManager.hpp"

#include "Panzer_BC.hpp"

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
   
   /** Checkout several responses. This method can only
     * be called once. This sets up the response library
     * to register responses. Basically this says I care
     * about these responses.
     */
   void checkoutResponses(const Teuchos::ParameterList & p);

   //! Print available reponses
   void printAvailableResponses(std::ostream & os) const;

   //! Print checked out reponses
   void printCheckedOutResponses(std::ostream & os) const;

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

   /** Register repsonses with a particular volume field manager. 
     */
   template <typename Traits,typename ScalarT>
   void registerResponses(const Teuchos::RCP<const Teuchos::Comm<int> > & comm,int worksetSize,
                          const std::string & blockId,PHX::FieldManager<Traits> & fm) const;

   //! Get the names of the available volume responses
   void getAvailableVolumeResponses(std::vector<std::string> & name) const;

   //! Get the names of the available volume responses (primarily for debugging)
   void getAvailableVolumeResponses(std::vector<std::pair<std::string,std::set<std::string> > > & responsePairs) const;

   //! Get the names of the checked out volume responses
   void getCheckedOutVolumeResponses(std::vector<std::string> & name) const;

   //! Get the names of the checked out volume responses, using the element block
   void getCheckedOutVolumeResponses(const std::string & eBlock,std::vector<std::string> & name) const;

   /** @} */

   //////////////////////////////////////////////////////////////////////
   /** \defgroup BC methods
     * Methods that act on boundary responses.
     * @{
     */

   //! Include an edge type response in the library
   void addResponse(const PHX::FieldTag & ft, const panzer::BC & bc);

   //! Get the names of the available BC responses
   void getAvailableBCResponses(std::vector<std::string> & name) const;

   //! Get the names of the checked out BC responses
   void getCheckedOutBCResponses(std::vector<std::string> & name) const;

   /** @} */

private:
   typedef std::map<std::string,std::set<std::string> > VolumeMap;

   //! Contains a library of responses and the element blocks they correspond to: field name -> block ids
   VolumeMap volumeResponses_;

   //! Contains a library of checkout responses ordered by element block: block id -> field names
   VolumeMap coVolumeResponses_;

   //! Has checkoutResponses been previously called?
   bool previousCheckout_;
};

inline std::ostream & operator<<(std::ostream & os,const ResponseLibrary & rl)
{ rl.print(os); return os; }

}

#endif
