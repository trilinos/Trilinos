/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/
#ifndef STK_UTIL_INCLUDE_DEPRECATIONWARNING_H_
#define STK_UTIL_INCLUDE_DEPRECATIONWARNING_H_

#include "stk_util/registry/VersionNumber.hpp"  // for VersionNumber
#include <sstream>                              // for ostringstream

namespace stk
{
namespace util
{

/* This class is used to add a deprecation warning in the code. Because Sierrra
 * uses even release number (4.48, 4.50, etc.) and odd development numbers
 * (4.49, 4.51, etc.), it is suggested to add deprecation warnings for the
 * odd development number preceeding the release in which it will be deprecated,
 * e.g., add 4.49 for a feature that will be deprecated (throw) in 4.50.
 *
 * Deletion of the warning and code is suggested one release ater, so
 * a 4.49 deprecation can be deleted in 4.51 versions.
 */
class DeprecationWarning
{
public:
  DeprecationWarning(const VersionNumber & removal_version);
  ~DeprecationWarning();

  template <class T> DeprecationWarning & operator<<(const T & t)
  {
    message << t;
    return *this;
  }

private:
  VersionNumber my_removal_version;
  std::ostringstream message;
};
}
}

#endif /* STK_UTIL_INCLUDE_DEPRECATIONWARNING_H_ */
