#ifndef __Teko_RequestCallback_hpp__
#define __Teko_RequestCallback_hpp__

#include "Teko_RequestMesg.hpp"

namespace Teko {

class RequestCallbackBase {
public:
   virtual ~RequestCallbackBase() {}

   virtual bool handlesRequest(const RequestMesg &) = 0;
};

template <typename DataT>
class RequestCallback : public RequestCallbackBase {
public:
   virtual ~RequestCallback() {}

   virtual DataT request(const RequestMesg &) = 0;
   virtual bool handlesRequest(const RequestMesg &) = 0;
};

}

#endif
