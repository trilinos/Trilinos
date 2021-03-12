/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#ifndef STK_COUPLING_impl_VERSION_UTILS_HPP
#define STK_COUPLING_impl_VERSION_UTILS_HPP

#include <stk_util/parallel/Parallel.hpp>
#include <string>
#include <cstdint>

namespace stk
{
namespace coupling
{
namespace impl
{

class OfficialCouplingVersion
{
public:
  OfficialCouplingVersion();
  void set_version();
protected:
  int m_couplingVersion;
};

enum CouplingCompatibilityMode : uint8_t {
  NotYetQueried = 0,
  Current,
  BackwardsCompatible,
  Deprecated,
  Incompatible
};

CouplingCompatibilityMode get_coupling_compatibility_mode(MPI_Comm comm);

class OfficialCompatibilityMode
{
public:
  OfficialCompatibilityMode();
  void exchange_app_versions(MPI_Comm comm);
  void set_compatibility_mode();
protected:
  OfficialCompatibilityMode(CouplingCompatibilityMode mode);
  CouplingCompatibilityMode m_couplingCompatibilityMode;
};

}
}
}

#endif /* STK_COUPLING_impl_VERSION_UTILS_HPP */
