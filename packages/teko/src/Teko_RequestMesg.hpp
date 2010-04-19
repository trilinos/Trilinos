#ifndef __Teko_RequestMesg_hpp__
#define __Teko_RequestMesg_hpp__

#include <string>

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"

namespace Teko {

class RequestMesg {
public:
   /** Construct a request messages specifing the
     * details of the request.
     *
     * \param[in] name Name of request to be satisfied.
     * \param[in] tag Optional tag describing other information
     *                about this tag.
     */
   explicit RequestMesg(const std::string & name,unsigned int tag=0)
      : name_(name), tag_(tag) {}

   /** Construct a parameter list message. This sets the 
     * name to "Parameter List" and the tag to 0. But now the
     * parameter list can completely describe the request.
     *
     * \param[in] pl Parameter list describing the request
     */
   explicit RequestMesg(const Teuchos::RCP<Teuchos::ParameterList> & pl)
      : name_("Parameter List"), tag_(0), paramList_(pl) {}

   //! Simple base class destructor
   virtual ~RequestMesg() {}

   //! Get the name for this request
   std::string getName() const
   { return name_; }

   //! Get the tag for this request
   unsigned int getTag() const
   { return tag_; }

   //! Get parameter list for this request
   const Teuchos::ParameterList & getParameterList() const
   { return *paramList_; }
   
protected:
   std::string name_;
   unsigned int tag_;
   Teuchos::RCP<Teuchos::ParameterList> paramList_;
};

} // end namespace Teko

#endif
