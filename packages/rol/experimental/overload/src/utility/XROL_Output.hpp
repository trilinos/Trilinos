#pragma once 

#include "cxxstd.hpp"

namespace XROL {

template<class Arg>
void to_ostream( std::ostream &os, const Arg& arg ) {
  os << arg;
}

template<class First, class ...Args>
void to_ostream( std::ostream &os, const First& first, Args&&... args) {
  os << first;
  to_ostream(os,args...);
}









} // namespace XROL
