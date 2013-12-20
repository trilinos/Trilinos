#ifndef PIKE_RESPONSE_HPP
#define PIKE_RESPONSE_HPP

#include "Pike_BlackBox_config.hpp"
#include "Teuchos_Describable.hpp"
#include "Teuchos_VerboseObject.hpp"
#include <string>
#include <iostream>

namespace pike {

  class Response : public Teuchos::Describable, public Teuchos::VerboseObject<pike::Response> {

  public:

    virtual ~Response() {}

    virtual std::string name() const = 0;

  };

  std::ostream& operator<<(std::ostream& os, const pike::Response& r);

}

inline 
std::ostream& pike::operator<<(std::ostream& os, const pike::Response& r)
{
  r.describe(os); 
  return os;
}

#endif
