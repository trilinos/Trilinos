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
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Andrey Prokopenko (aprokop@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_BLACKBOXPFACTORY_DEF_HPP
#define MUELU_BLACKBOXPFACTORY_DEF_HPP

#include <stdlib.h>
#include <iomanip>


// #include <Teuchos_LAPACK.hpp>
#include <Teuchos_SerialDenseMatrix.hpp>
#include <Teuchos_SerialDenseVector.hpp>
#include <Teuchos_SerialDenseSolver.hpp>

#include <Xpetra_CrsMatrixWrap.hpp>
#include <Xpetra_ImportFactory.hpp>
#include <Xpetra_Matrix.hpp>
#include <Xpetra_MapFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_VectorFactory.hpp>

#include <Xpetra_IO.hpp>

#include "MueLu_FactoryManagerBase.hpp"
#include "MueLu_BlackBoxPFactory_decl.hpp"

#include "MueLu_MasterList.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    // Coarsen can come in two forms, either a single char that will be interpreted as an integer which is used as the coarsening
    // rate in every spatial dimentions, or it can be a longer string that will then be interpreted as an array of integers.
    // By default coarsen is set as "{2}", hence a coarsening rate of 2 in every spatial dimension is the default setting!
    validParamList->set<std::string >           ("Coarsen",                 "{2}", "Coarsening rate in all spatial dimensions");
    validParamList->set<RCP<const FactoryBase> >("A",                         Teuchos::null, "Generating factory of the matrix A");
    validParamList->set<RCP<const FactoryBase> >("Nullspace",                 Teuchos::null, "Generating factory of the nullspace");
    validParamList->set<RCP<const FactoryBase> >("Coordinates",               Teuchos::null, "Generating factory for coorindates");
    validParamList->set<RCP<const FactoryBase> >("gNodesPerDim",              Teuchos::null, "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");
    validParamList->set<RCP<const FactoryBase> >("lNodesPerDim",              Teuchos::null, "Number of nodes per spatial dimmension provided by CoordinatesTransferFactory.");
    validParamList->set<std::string >           ("axisPermutation",         "", "Assuming a global (x,y,z) orientation, local might be (z,y,x). This vector gives a permutation from global to local orientation.");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& fineLevel, Level& coarseLevel) const {
    Input(fineLevel, "A");
    Input(fineLevel, "Nullspace");
    Input(fineLevel, "Coordinates");
    // Request the global number of nodes per dimensions
    if(fineLevel.GetLevelID() == 0) {
      if(fineLevel.IsAvailable("gNodesPerDim", NoFactory::get())) {
        fineLevel.DeclareInput("gNodesPerDim", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(fineLevel.IsAvailable("gNodesPerDim", NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "gNodesPerDim was not provided by the user on level0!");
      }
    } else {
      Input(fineLevel, "gNodesPerDim");
    }

    // Request the local number of nodes per dimensions
    if(fineLevel.GetLevelID() == 0) {
      if(fineLevel.IsAvailable("lNodesPerDim", NoFactory::get())) {
        fineLevel.DeclareInput("lNodesPerDim", NoFactory::get(), this);
      } else {
        TEUCHOS_TEST_FOR_EXCEPTION(fineLevel.IsAvailable("lNodesPerDim", NoFactory::get()),
                                   Exceptions::RuntimeError,
                                   "lNodesPerDim was not provided by the user on level0!");
      }
    } else {
      Input(fineLevel, "lNodesPerDim");
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& fineLevel, Level& coarseLevel) const {
    return BuildP(fineLevel, coarseLevel);

  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlackBoxPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BuildP(Level& fineLevel, Level& coarseLevel) const {
    FactoryMonitor m(*this, "Build", coarseLevel);

  }

} //namespace MueLu

#define MUELU_BLACKBOXPFACTORY_SHORT
#endif // MUELU_BLACKBOXPFACTORY_DEF_HPP
