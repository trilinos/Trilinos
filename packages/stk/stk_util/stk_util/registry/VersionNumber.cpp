/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_util/registry/VersionNumber.hpp>
#include <stk_util/registry/ProductRegistry.hpp>
#include <stdexcept>
#include <string>

namespace stk
{
namespace util
{

namespace
{
static VersionNumber current_version_number(0, 0);
}

VersionNumber VersionNumber::current_version()
{
  static bool versionIsSet = false;
  if(!versionIsSet) {
    VersionNumber::set_current_version(stk::ProductRegistry::version());
    versionIsSet = true;
  }
  return current_version_number;
}

void VersionNumber::set_current_version(const VersionNumber & v) { current_version_number = v; }

void VersionNumber::set_current_version(const std::string & version_string)
{
  std::runtime_error err("Failed to parse version string " + version_string);

  std::string delimiter = ".";
  unsigned indexOfFirst = version_string.find(delimiter);
  unsigned indexOfSecond = version_string.find(delimiter, indexOfFirst);

  const int major = std::stoi(version_string.substr(0, indexOfFirst));
  const int minor = std::stoi(version_string.substr(indexOfFirst+1, indexOfSecond));

  set_current_version(VersionNumber{major, minor});
}

bool operator<(const VersionNumber & lhs, const VersionNumber & rhs)
{
  if (lhs.major() == rhs.major())
  {
    return lhs.minor() < rhs.minor();
  }
  return lhs.major() < rhs.major();
}

bool operator==(const VersionNumber & lhs, const VersionNumber & rhs)
{
  return lhs.major() == rhs.major() && lhs.minor() == rhs.minor();
}

std::ostream & operator<<(std::ostream & os, const VersionNumber & v)
{
  os << "Version " << v.major() << "." << v.minor();
  return os;
}
}
}
