#ifndef SIMPLECXX_HELLO_WOLRD_HPP
#define SIMPLECXX_HELLO_WOLRD_HPP


#include <ostream>
#include "SimpleCxx_config.h"


namespace SimpleCxx {


/* \brief Simple hello world class.
 */
class HelloWorld {
public:
  /** \brief . */
  void printHelloWorld(std::ostream &out) const;
};


} // namespace SimpleCxx


#endif // SIMPLECXX_HELLO_WOLRD_HPP
