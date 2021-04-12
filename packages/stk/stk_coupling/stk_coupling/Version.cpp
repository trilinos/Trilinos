/*--------------------------------------------------------------------*/
/*    Copyright 2002 - 2008, 2010, 2011 National Technology &         */
/*    Engineering Solutions of Sandia, LLC (NTESS). Under the terms   */
/*    of Contract DE-NA0003525 with NTESS, there is a                 */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include <stk_coupling/Version.hpp>
#include <stk_coupling/SyncInfo.hpp>
#include <stk_coupling/Constants.hpp>
#include <stk_coupling/impl_VersionUtils.hpp>
#include <stk_util/parallel/ParallelReduce.hpp>

namespace stk
{
namespace coupling
{

namespace
{
  static constexpr int OFFICIAL_VERSION = 1;
  static constexpr int BAD_VERSION = -99;
  static int staticCouplingVersion = BAD_VERSION;

  static impl::CouplingCompatibilityMode staticCouplingCompatibilityMode = impl::NotYetQueried;
}

int version()
{
  bool versionIsSet = !(staticCouplingVersion == BAD_VERSION);
  if (!versionIsSet) {
    impl::OfficialCouplingVersion().set_version();
  }

  return staticCouplingVersion;
}

namespace impl {

CouplingCompatibilityMode get_coupling_compatibility_mode(MPI_Comm comm)
{
  bool modeHasBeenSet = !(staticCouplingCompatibilityMode == NotYetQueried);
  if (!modeHasBeenSet) {
    OfficialCompatibilityMode officialCompatibilityMode;
    officialCompatibilityMode.exchange_app_versions(comm);
    officialCompatibilityMode.set_compatibility_mode();
  }
  return staticCouplingCompatibilityMode;
}

OfficialCouplingVersion::OfficialCouplingVersion()
 : m_couplingVersion(OFFICIAL_VERSION)
{}

void OfficialCouplingVersion::set_version()
{
  staticCouplingVersion = m_couplingVersion;
}

OfficialCompatibilityMode::OfficialCompatibilityMode()
: m_couplingCompatibilityMode(NotYetQueried)
{}

OfficialCompatibilityMode::OfficialCompatibilityMode(CouplingCompatibilityMode mode)
: m_couplingCompatibilityMode(mode)
{}

void OfficialCompatibilityMode::exchange_app_versions(MPI_Comm comm)
{
    int minVersion = version();
    int maxVersion = version();
    stk::all_reduce(comm, stk::ReduceMin<1>(&minVersion) & stk::ReduceMax<1>(&maxVersion));
    int difference = maxVersion - minVersion;
    switch(difference) {
      case 0:
        m_couplingCompatibilityMode = Current;
        break;
      case 1:
        m_couplingCompatibilityMode = ((minVersion==version()) ? Deprecated : BackwardsCompatible);
        break;
      default:
        m_couplingCompatibilityMode = Incompatible;
        break;
    }
}

void OfficialCompatibilityMode::set_compatibility_mode()
{
  staticCouplingCompatibilityMode = m_couplingCompatibilityMode;
}

} // namespace impl
} // namespace coupling
} // namespace stk
