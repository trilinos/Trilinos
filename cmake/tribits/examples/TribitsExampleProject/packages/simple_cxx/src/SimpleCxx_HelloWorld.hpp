#ifndef SIMPLECXX_HELLO_WOLRD_HPP
#define SIMPLECXX_HELLO_WOLRD_HPP


#include <ostream>
#include "SimpleCxx_config.h"


namespace SimpleCxx {

/** \brief . */
std::string deps();

/* \brief Simple hello world class.
 */
class HelloWorld {
public:
  /** \brief. */
  HelloWorld();
  /** \brief . */
  void printHelloWorld(std::ostream &out) const;
  /** \brief Deprecated. */
  SIMPLECXX_DEPRECATED int someOldFunc() const; 
  /** \brief Deprecated. */
  SIMPLECXX_DEPRECATED_MSG("Just don't call this function at all please!")
  int someOldFunc2() const; 
};


} // namespace SimpleCxx


#endif // SIMPLECXX_HELLO_WOLRD_HPP
