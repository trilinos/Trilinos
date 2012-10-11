#ifndef __Panzer_Response_Functional_hpp__
#define __Panzer_Response_Functional_hpp__

#include <string>

#include "Panzer_ResponseBase.hpp"

namespace panzer {

/** This class provides a response storage for
  * simple functionals of the solution (i.e. scalar
  * values).
  */
template <typename ScalarT>
class Response_Functional : public ResponseBase {
public:

   Response_Functional(const std::string & responseName)
     : ResponseBase(responseName) {}

   //! provide direct access, this thing is pretty simple
   ScalarT value;

private:
   // hide these methods
   Response_Functional();
   Response_Functional(const Response_Functional &);
};

}

#endif
