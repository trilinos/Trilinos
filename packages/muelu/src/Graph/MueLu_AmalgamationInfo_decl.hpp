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
 * MueLu_AmalgamationInfo_decl.hpp
 *
 *  Created on: Mar 28, 2012
 *      Author: wiesner
 */

#ifndef MUELU_AMALGAMATIONINFO_DECL_HPP_
#define MUELU_AMALGAMATIONINFO_DECL_HPP_

#include <Xpetra_ConfigDefs.hpp>   // global_size_t
#include <Xpetra_Map_fwd.hpp>
#include <Xpetra_MapFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"

#include "MueLu_BaseClass.hpp"
#include "MueLu_AmalgamationInfo_fwd.hpp"

namespace MueLu {

/*!
  @class AmalgamationInfo
  @brief minimal container class for storing amalgamation information


  stores map of global node id on current processor to global DOFs ids on current processor
  nodegid2dofgids_

  that is used for unamalgamation
*/

  template <class LocalOrdinal  = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class AmalgamationInfo
    : public BaseClass {
#undef MUELU_AMALGAMATIONINFO_SHORT
#include "MueLu_UseShortNamesOrdinal.hpp"

  public:

    AmalgamationInfo() {
      nodegid2dofgids_ = Teuchos::null;
    }

    virtual ~AmalgamationInfo() {}

    /// Return a simple one-line description of this object.
    std::string description() const;

    //! Print the object with some verbosity level to an FancyOStream object.
    //using MueLu::Describable::describe; // overloading, not hiding
    //void describe(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;;
    void print(Teuchos::FancyOStream &out, const VerbLevel verbLevel = Default) const;

    void SetAmalgamationParams(RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > nodegid2dofgids) const;

    RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > GetGlobalAmalgamationParams() const;

    void SetNodeGIDVector(RCP<std::vector<GlobalOrdinal> > nodegids) const;

    RCP<std::vector<GlobalOrdinal> > GetNodeGIDVector() const;

    void SetNumberOfNodes(GlobalOrdinal cntNodes) const { cntNodes_ = cntNodes; }
    GlobalOrdinal GetNumberOfNodes() const { return cntNodes_; }

  private:

    //! @name amalgamation information variables
    //@{

    // map of global node id on current processor to global DOFs ids on current processor
    mutable RCP<std::map<GlobalOrdinal,std::vector<GlobalOrdinal> > > nodegid2dofgids_; //< used for building overlapping ImportDofMap

    // contains global node ids on current proc (used by CoalesceDropFactory to build nodeMap)
    mutable RCP<std::vector<GlobalOrdinal> > gNodeIds_;

    /// number of nodes on current proc
    mutable GlobalOrdinal cntNodes_;

    //@}

  };

} // namespace MueLu

#define MUELU_AMALGAMATIONINFO_SHORT
#endif /* MUELU_AMALGAMATIONINFO_DECL_HPP_ */
