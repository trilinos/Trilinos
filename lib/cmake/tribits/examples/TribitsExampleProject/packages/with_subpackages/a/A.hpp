#ifndef PACKAGEWITHSUBPACKAGES_A_HPP_
#define PACKAGEWITHSUBPACKAGES_A_HPP_

#include <string>

namespace WithSubpackages {

  // return a string containing "A"
  std::string getA();

  // return a string describing the dependencies of "A", recursively
  std::string depsA();

  // return special value
  int specialValue();

}


#endif /* PACKAGEWITHSUBPACKAGES_A_HPP_ */
