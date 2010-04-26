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

   /** Send a pre-request message to the callback
     * allowing them to do some work up front and
     * ahead of time. This is meant to be called
     * at construction/initialization.
     *
     * \param[in] rm The message describing the request.
     */
   template <typename DataT>
   void preRequest(const RequestMesg & rm) const;

private:
   // stores the callbacks to be used by this handler.
   mutable std::vector<Teuchos::RCP<RequestCallbackBase> > callbacks_;

   // hidden from the user
   RequestHandler(const RequestHandler & rh);
   RequestHandler & operator=(const RequestHandler &);
   const RequestHandler & operator=(const RequestHandler &) const;
};

#include "Teko_RequestHandler_impl.hpp"

}

#endif
