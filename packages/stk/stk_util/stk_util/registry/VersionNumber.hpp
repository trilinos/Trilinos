/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef STK_UTIL_INCLUDE_VERSION_H_
#define STK_UTIL_INCLUDE_VERSION_H_

#include <ostream>  // for ostream
#include <string>   // for string

#undef major
#undef minor

namespace stk 
{
namespace util
{

class VersionNumber
{
public:
  VersionNumber(int maj, int min) : major_val(maj), minor_val(min) {}

  int major() const { return major_val; }
  int minor() const { return minor_val; }
  bool is_release() const { return minor_val % 2 == 0; }

  static VersionNumber current_version();
  static void set_current_version(const VersionNumber & v);
  static void set_current_version(const std::string & version_string);

private:
  int major_val;
  int minor_val;
};

bool operator<(const VersionNumber & lhs, const VersionNumber & rhs);
bool operator==(const VersionNumber & lhs, const VersionNumber & rhs);
std::ostream & operator<<(std::ostream & os, const VersionNumber & v);
}
}

#endif /* STK_UTIL_INCLUDE_VERSION_H_ */
