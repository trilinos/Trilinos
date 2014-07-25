#ifndef PACKAGEWITHSUBPACKAGES_C_HPP_
#define PACKAGEWITHSUBPACKAGES_C_HPP_

#include <string>

namespace PackageWithSubpackages {

  // return a string containing "C"
  std::string getC();

  // return a string describing the dependencies of "C", recursively
  std::string depsC();

}


#endif /* PACKAGEWITHSUBPACKAGES_C_HPP_ */
