#ifndef PACKAGEWITHSUBPACKAGES_B_HPP_
#define PACKAGEWITHSUBPACKAGES_B_HPP_

#include "WithSubpackagesB_config.h"

#include <string>

namespace WithSubpackages {

  // return a string containing "B"
  std::string getB();

  // return a string describing the dependencies of "B", recursively
  std::string depsB();

}


#endif /* PACKAGEWITHSUBPACKAGES_B_HPP_ */
