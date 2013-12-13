#ifndef PIKE_RESPONSE_HPP
#define PIKE_RESPONSE_HPP

#include <string>
#include <iostream>

namespace pike {

  class Response {

    virtual ~Response();

    virtual std::string name() const = 0;

    virtual std::ostream& print(std::ostream& os) const = 0;

  };

  std::ostream& operator<<(std::ostream& os);

}

#endif
