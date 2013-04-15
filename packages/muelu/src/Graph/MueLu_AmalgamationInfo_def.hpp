// @HEADER
//
// ***********************************************************************
//
//        MueLu: A package for multigrid based preconditioning
//                  Copyright 2012 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
/*
 * MueLu_AmalgamationInfo_def.hpp
 *
 *  Created on: Mar 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AMALGAMATIONINFO_DEF_HPP_
#define MUELU_AMALGAMATIONINFO_DEF_HPP_

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
