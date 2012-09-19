/**   ------------------------------------------------------------
 *    Copyright 2001 - 2009 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 *    ------------------------------------------------------------
 */

#include <iostream>
#include <stk_util/util/IndentStreambuf.hpp>

namespace sierra {

std::ostream &
out() {
  static std::ostream s_out(std::cout.rdbuf());

  return s_out;
}


std::ostream &
pout() {
  static std::ostream s_pout(std::cout.rdbuf());

  return s_pout;
}


std::ostream &
dout() {
  static std::ostream s_dout(std::cout.rdbuf());

  return s_dout;
}


std::ostream &
tout() {
  static std::ostream s_tout(std::cout.rdbuf());

  return s_tout;
}


std::ostream &
dwout() {
  static stk::indent_streambuf s_dwoutStreambuf(std::cout.rdbuf());
  static std::ostream s_dwout(&s_dwoutStreambuf);
  
  return s_dwout;
}

} // namespace sierra
