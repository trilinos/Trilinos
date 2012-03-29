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
  void AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::SetAmalgamationParams(RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > globalamalblockid2myrowid,RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > globalamalblockid2globalrowid) const {
    globalamalblockid2myrowid_ = globalamalblockid2myrowid;
    globalamalblockid2globalrowid_ = globalamalblockid2globalrowid;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<std::map<GlobalOrdinal,std::vector<LocalOrdinal> > > AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetMyAmalgamationParams() const {
    return globalamalblockid2myrowid_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::GetGlobalAmalgamationParams() const {
    return globalamalblockid2globalrowid_;
  }

  template <class LocalOrdinal, class GlobalOrdinal, class Node, class LocalMatOps>
  std::string AmalgamationInfo<LocalOrdinal, GlobalOrdinal, Node, LocalMatOps>::description() const {
    return "AmalgamationInfo";
  }
}


#endif /* MUELU_AMALGAMATIONINFO_DEF_HPP_ */
