// Copyright(C) 1999-2021 National Technology & Engineering Solutions
// of Sandia, LLC (NTESS).  Under the terms of Contract DE-NA0003525 with
// NTESS, the U.S. Government retains certain rights in this software.
//
// See packages/seacas/LICENSE for details

#ifndef IOSS_IOVS_CATALYST_VERSION_H
#define IOSS_IOVS_CATALYST_VERSION_H

#ifndef __CATALYST_PLUGIN_BUILD
#include "iovs_export.h"
#else
#define IOVS_EXPORT
#endif

#include <string>

namespace Iovs {

  class IOVS_EXPORT CatalystVersion
  {
  public:
    const std::string iossCatalystInterfaceVersion = "2.0.0";
    std::string       getIOSSCatalystInterfaceVersion();

    // The IOSS interface and IOSS Catalyst plugin are versioned
    // using semantic versioning: MAJOR.MINOR.PATCH
    //
    // The interface and plugin versions are incompatible if:
    //
    //    The MAJOR versions are not equal.
    //                     (or)
    //    The MAJOR versions are equal, and the interface MINOR version
    //    is greater than the plugin MINOR version.
    //                     (or)
    //    The MAJOR versions are equal, the MINOR versions are equal,
    //    and the interface PATCH version is greater than the plugin PATCH
    //    version.

    bool isInterfaceCompatibleWithPlugin(const std::string &interface_version,
                                         const std::string &plugin_version);

  private:
    const unsigned int SEMANTIC_VERSION_LENGTH = 3;
    const unsigned int MAJOR_INDEX             = 0;
    const unsigned int MINOR_INDEX             = 1;
    const unsigned int PATCH_INDEX             = 2;
  };

} // namespace Iovs

#endif
