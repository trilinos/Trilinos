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
#ifndef MUELU_PARAMETERLISTINTERPRETER_DECL_HPP
#define MUELU_PARAMETERLISTINTERPRETER_DECL_HPP

#include <string>

#include <Teuchos_Array.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Xpetra_MultiVectorFactory_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_HierarchyFactory_fwd.hpp"
#include "MueLu_HierarchyManager.hpp"

namespace MueLu {

  template <class Scalar = double, class LocalOrdinal = int, class GlobalOrdinal = LocalOrdinal, class Node = Kokkos::DefaultNode::DefaultNodeType, class LocalMatOps = typename Kokkos::DefaultKernels<void,LocalOrdinal,Node>::SparseOps>
  class ParameterListInterpreter : public HierarchyManager<Scalar, LocalOrdinal, GlobalOrdinal, Node, LocalMatOps> {
#undef MUELU_PARAMETERLISTINTERPRETER_SHORT
#include "MueLu_UseShortNames.hpp"

  public:

    //@{

    //!
    ParameterListInterpreter() { }

    //!
    ParameterListInterpreter(Teuchos::ParameterList & paramList);

    //!
    ParameterListInterpreter(const std::string & xmlFileName);

    //! Destructor.
    virtual ~ParameterListInterpreter() { }

    //@}

    //@{

    void SetParameterList(const Teuchos::ParameterList & paramList);

    //@}

  private:
    typedef std::map<std::string, RCP<const FactoryBase> > FactoryMap; //TODO: remove this line

    void BuildFactoryMap(const Teuchos::ParameterList & paramList, const FactoryMap & factoryMapIn, FactoryMap & factoryMapOut) const;

    //@{ Matrix configuration

    //! Setup Matrix object
    virtual void SetupMatrix(Matrix & Op) const;

    //! Setup extra data
    virtual void SetupExtra(Hierarchy & H) const;

    // Matrix configuration storage
    Teuchos::ParameterList operatorList_; //TODO: should it be stored in another format to avoid xml parsing in SetupMatrix()?

    //@}

  }; // class

} // namespace MueLu

#define MUELU_PARAMETERLISTINTERPRETER_SHORT
#endif // MUELU_PARAMETERLISTINTERPRETER_DECL_HPP

// TODO:
// - parameter list validator
// - SetParameterList
// - Set/Get directly Level manager
// - build per level
// - comments/docs
// - use FactoryManager instead of FactoryMap
