#ifndef __Panzer_Response_hpp__
#define __Panzer_Response_hpp__

#include "Panzer_config.hpp"
#include "Panzer_Traits.hpp"

#include "Sacado_Traits.hpp"

namespace panzer {

// Forward declaration
template <typename RespType> class Response;

template < >
class Response <panzer::Traits::Value> {
public:
   typedef panzer::Traits::Value RespType;

   Response(Sacado::ScalarType<RespType::ScalarT>::type value)
      : value_(value) {}

   //! get value of this response
   Sacado::ScalarType<RespType::ScalarT>::type getValue() const
   { return value_; }

private:
   // hide these constructors
   Response();
   Response(const Response<panzer::Traits::Value> &);
    
   Sacado::ScalarType<RespType::ScalarT>::type value_;
};

template < >
class Response <panzer::Traits::Derivative> { };

}

#endif
