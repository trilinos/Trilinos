#include "visualization/utils/Iovs_CatalystVersion.h"
#include <iostream>
#include <sstream>
#include <vector>

namespace Iovs {

  std::string CatalystVersion::getIOSSCatalystInterfaceVersion()
  {
    return iossCatalystInterfaceVersion;
  }

  bool CatalystVersion::isInterfaceCompatibleWithPlugin(const std::string &interface_version,
                                                        const std::string &plugin_version)
  {
    bool               retVal = false;
    std::istringstream interface_version_parser(interface_version);
    std::istringstream plugin_version_parser(plugin_version);
    std::vector<int>   iver(SEMANTIC_VERSION_LENGTH, 0);
    std::vector<int>   pver(SEMANTIC_VERSION_LENGTH, 0);

    for (unsigned int i = 0; i < SEMANTIC_VERSION_LENGTH; i++) {
      interface_version_parser >> iver[i];
      plugin_version_parser >> pver[i];
      interface_version_parser.get();
      plugin_version_parser.get();
    }

    if (iver[MAJOR_INDEX] == pver[MAJOR_INDEX]) {
      if (iver[MINOR_INDEX] < pver[MINOR_INDEX]) {
        retVal = true;
      }
      else if (iver[MINOR_INDEX] == pver[MINOR_INDEX]) {
        if (iver[PATCH_INDEX] <= pver[PATCH_INDEX]) {
          retVal = true;
        }
      }
    }

    return retVal;
  }

} // namespace Iovs
