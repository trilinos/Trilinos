#ifndef __Teko_RequestHandler_hpp__
#define __Teko_RequestHandler_hpp__

#include <vector>

#include "Teko_RequestCallback.hpp"

#include "Teuchos_RCP.hpp"

namespace Teko {

/** Classes that handles and distrubutes requests.
  * This has two types of users. Those that register
  * callbacks to handle the requests, and those that
  * make requests.
  *
  * This class is passive-aggressive.  It calls required
  * data a "request", however it actually requires 
  * a handler to satisfy the request.
  *
  * \note The implemntation is based on the listener pattern.
  */
class RequestHandler {
public:
   explicit RequestHandler();


   /** Add a call back object to handle requests
     *
     * \param[in] callback Ref-count-pointer to a call back object.
     */
   void addRequestCallback(const Teuchos::RCP<RequestCallbackBase> & callback);

   /** Get the data for a particular request.
     *
     * \param[in] rm The message describing the request.
     */
   template <typename DataT>
   DataT request(const RequestMesg & rm) const;

private:
   // stores the callbacks to be used by this handler.
   mutable std::vector<Teuchos::RCP<RequestCallbackBase> > callbacks_;
};

#include "Teko_RequestHandler_impl.hpp"

}

#endif
