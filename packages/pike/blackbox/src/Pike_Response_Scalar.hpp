#ifndef PIKE_RESPONSE_SCALAR_HPP
#define PIKE_RESPONSE_SCALAR_HPP

#include "Teuchos_RCP.hpp"
#include "Pike_Response.hpp"
#include <string>
#include <iostream>

namespace pike {

  template<typename Scalar>
  class ScalarResponse : public pike::Response {

  public:

    ScalarResponse(const std::string name) : name_(name) { }

    std::string name() const
    { return name_; } 

    Scalar get() const
    { return value_; }
    
    void set(const Scalar& value)
    { value_ = value; }

  private:

    std::string name_;
    Scalar value_;

  };

  /** \brief non-member ctor
      \relates ScalarResponse
  */
  template<typename Scalar>
  Teuchos::RCP<ScalarResponse<Scalar> > scalarResponse(const std::string name = "")
  {
    return Teuchos::rcp(new pike::ScalarResponse<Scalar>(name));
  }

  /** \brief non-member accessor
      \relates ScalarResponse
  */
  template<typename Scalar>
  void getScalarResponse(Scalar& value,
			 const pike::Response& response) {
    value = (dynamic_cast<const pike::ScalarResponse<Scalar>& >(response)).get();
  }

  /** \brief non-member accessor
      \relates ScalarResponse
  */
  template<typename Scalar>
  Scalar getScalarResponse(const pike::Response& response) {
    return (dynamic_cast<const pike::ScalarResponse<Scalar>& >(response)).get();
  }

}

#endif
