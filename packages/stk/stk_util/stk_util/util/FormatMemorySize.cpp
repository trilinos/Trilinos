/*------------------------------------------------------------------------*/
/*                 Copyright 2010, 2011 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#include <sstream>
#include <iomanip>
#include <cmath>

#include <boost/lexical_cast.hpp>

#include <stk_util/util/FormatMemorySize.hpp>

namespace stk {

std::string
formatMemorySize(
  double                size)
{
  std::string           result;
  
  static const double kb = 1024.0;
  // static const double mb = kb * kb;
  // static const double gb = kb * kb * kb;

  if (size < 0.0) {
    result = "-";
    size = -size;
  }

  // output size in kilo bytes
  result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size / kb));
  result += " KB";
  // if (size < kb) {
  //   // output size in bytes
  //   result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size));
  //   result += " B";
  // }
  // else if (size < mb) {
  //   // output size in kilo bytes
  //   result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size / kb));
  //   result += " KB";
  // }
  // else if (size < gb) {
  //   // output size in mega bytes
  //   result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size / mb));
  //   result += " MB";
  // }
  // else {
  //   // everything else output in giga bytes
  //   result += boost::lexical_cast<std::string>(static_cast<unsigned long>(size / gb));
  //   result += " GB";
  // }
  
  return result;
}


std::string
formatMemorySize(
  MemorySize            size)
{
  std::string           result;
  
  static const MemorySize kb = 1024;
  // static const MemorySize mb = kb * kb;
  // static const MemorySize gb = kb * kb * kb;

  // output size in kilo bytes
  result = boost::lexical_cast<std::string>(size / kb);
  result += " KB";
  
  // if (size < kb) {
  //   // output size in bytes
  //   result = boost::lexical_cast<std::string>(size);
  //   result += " B";
  // }
  // else if (size < mb) {
  //   // output size in kilo bytes
  //   result = boost::lexical_cast<std::string>(size / kb);
  //   result += " KB";
  // }
  // else if (size < gb) {
  //   // output size in mega bytes
  //   result = boost::lexical_cast<std::string>(size / mb);
  //   result += " MB";
  // }
  // else {
  //   // everything else output in giga bytes
  //   result = boost::lexical_cast<std::string>(size / gb);
  //   result += " GB";
  // }
  
  return result;
}

} // namespace stk
