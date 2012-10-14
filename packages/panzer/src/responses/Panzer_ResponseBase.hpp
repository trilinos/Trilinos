#ifndef __Panzer_ResponseBase_hpp__
#define __Panzer_ResponseBase_hpp__

#include <string>

namespace panzer {

// *** Derive, do not modify this interface ***

/** This object serves as the base class which
  * the user derives from to build there response.
  * Note that by itself the interface for this is 
  * not useful for storing data.
  */
class ResponseBase : public GlobalEvaluationData_Default {
public:

   /** Only available constructor for this object.
     * This gurantees that the reponse object is
     * instantiatied with a name, and gives us basic
     * capability for name lookups.
     */
   ResponseBase(const std::string & responseName) :
    responseName_(responseName) {}

   virtual ~ResponseBase() {}

   /** Get the unmodified name for this response.
     */
   std::string getName() const { return responseName_; }

   /** Get the name of this response useful
     * for looking up data containers.
     */
   std::string getLookupName() const { return buildLookupName(responseName_); }

   /** Static member to build consisitent look up names
     * based on a response name.
     */
   static std::string buildLookupName(const std::string & responseName)
   { return "RESPONSE_"+responseName; }

   //! Inherited from GlobalEvaluationData, 
   virtual void ghostToGlobal(int) 
   { scatterResponse(); }

   virtual void intializeData()
   { initializeResponse(); }

   //! Prepare the response for access by the user (do global communication)
   virtual void scatterResponse() = 0;

   virtual void initializeResponse() = 0;

private:

   std::string responseName_; 
 
   // hide these methods
   ResponseBase();
   ResponseBase(const ResponseBase &);
};

}

#endif
