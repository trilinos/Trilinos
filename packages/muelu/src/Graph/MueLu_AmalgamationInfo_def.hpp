/*
 * MueLu_AmalgamationInfo_def.hpp
 *
 *  Created on: Mar 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AMALGAMATIONINFO_DEF_HPP_
#define MUELU_AMALGAMATIONINFO_DEF_HPP_

#include <Xpetra_MapFactory.hpp>

#include "MueLu_AmalgamationInfo_decl.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetAmalgamationParams(RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > nodegid2dofgids) const {
    nodegid2dofgids_ = nodegid2dofgids;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetGlobalAmalgamationParams() const {
    return nodegid2dofgids_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  void AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetNodeGIDVector(RCP<std::vector<GlobalOrdinal> > nodegids) const {
    gNodeIds_ = nodegids;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<std::vector<GlobalOrdinal> > AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetNodeGIDVector() const {
    return gNodeIds_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    return "AmalgamationInfo";
  }
}


#endif /* MUELU_AMALGAMATIONINFO_DEF_HPP_ */
