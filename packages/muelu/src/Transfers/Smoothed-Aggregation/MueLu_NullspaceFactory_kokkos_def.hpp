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
#ifndef MUELU_NULLSPACEFACTORY_KOKKOS_DEF_HPP
#define MUELU_NULLSPACEFACTORY_KOKKOS_DEF_HPP

#include "MueLu_NullspaceFactory_kokkos_decl.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_Level.hpp"
#include "MueLu_Monitor.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  RCP<const ParameterList> NullspaceFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A",          Teuchos::null, "Generating factory of the fine level matrix (only needed if default null space is generated)");
    validParamList->set< RCP<const FactoryBase> >("Nullspace",  Teuchos::null, "Generating factory of the fine level null space");

    // TODO not very elegant.
    validParamList->set< std::string >("Fine level nullspace", "Nullspace", "Variable name which is used to store null space multi vector on the finest level (default=\"Nullspace\"). For block matrices also \"Nullspace1\" to \"Nullspace9\" are accepted to describe the null space vectors for the (i,i) block (i=1..9).");

    // 1/20/2016: we could add a sublist (e.g. "Nullspaces" which is excluded from parameter validation)
    validParamList->set< RCP<const FactoryBase> >("Nullspace1", Teuchos::null, "Generating factory of the fine level null space associated with the (1,1) block in your n x n block matrix.");
    validParamList->set< RCP<const FactoryBase> >("Nullspace2", Teuchos::null, "Generating factory of the fine level null space associated with the (2,2) block in your n x n block matrix.");
    validParamList->set< RCP<const FactoryBase> >("Nullspace3", Teuchos::null, "Generating factory of the fine level null space associated with the (3,3) block in your n x n block matrix.");
    validParamList->set< RCP<const FactoryBase> >("Nullspace4", Teuchos::null, "Generating factory of the fine level null space associated with the (4,4) block in your n x n block matrix.");
    validParamList->set< RCP<const FactoryBase> >("Nullspace5", Teuchos::null, "Generating factory of the fine level null space associated with the (5,5) block in your n x n block matrix.");
    validParamList->set< RCP<const FactoryBase> >("Nullspace6", Teuchos::null, "Generating factory of the fine level null space associated with the (6,6) block in your n x n block matrix.");
    validParamList->set< RCP<const FactoryBase> >("Nullspace7", Teuchos::null, "Generating factory of the fine level null space associated with the (7,7) block in your n x n block matrix.");
    validParamList->set< RCP<const FactoryBase> >("Nullspace8", Teuchos::null, "Generating factory of the fine level null space associated with the (8,8) block in your n x n block matrix.");
    validParamList->set< RCP<const FactoryBase> >("Nullspace9", Teuchos::null, "Generating factory of the fine level null space associated with the (9,9) block in your n x n block matrix.");

    return validParamList;
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void NullspaceFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::DeclareInput(Level& currentLevel) const {
    const ParameterList& pL = GetParameterList();
    std::string nspName = pL.get<std::string>("Fine level nullspace");

    // only request "A" in DeclareInput if
    // 1) there is not nspName (e.g. "Nullspace") is available in Level, AND
    // 2) it is the finest level (i.e. LevelID == 0)
    if (currentLevel.IsAvailable(nspName, NoFactory::get()) == false && currentLevel.GetLevelID() == 0)
      Input(currentLevel, "A");

    if (currentLevel.GetLevelID() != 0) {
      // validate nullspaceFact_
      // 1) The factory for "Nullspace" (or nspName) must not be Teuchos::null, since the default factory
      //    for "Nullspace" is a NullspaceFactory
      // 2) The factory for "Nullspace" (or nspName) must be a TentativePFactory or any other TwoLevelFactoryBase derived object
      //    which generates the variable "Nullspace" as output
      TEUCHOS_TEST_FOR_EXCEPTION(GetFactory(nspName).is_null(), Exceptions::RuntimeError,
        "MueLu::NullspaceFactory::DeclareInput(): You must declare an existing factory which "
        "produces the variable \"Nullspace\" in the NullspaceFactory (e.g. a TentativePFactory).");
      currentLevel.DeclareInput("Nullspace", GetFactory(nspName).get(), this); /* ! "Nullspace" and nspName mismatch possible here */
    }
  }

  template<class NullspaceType, class LO>
  class NullspaceFunctor {
  private:
    NullspaceType                               nullspace;
    LO                                          numPDEs;
    typedef typename NullspaceType::value_type  SC;
    typedef Kokkos::ArithTraits<SC>             ATS;

  public:
    NullspaceFunctor(NullspaceType nullspace_, LO numPDEs_) :
      nullspace(nullspace_),
      numPDEs(numPDEs_)
    { }

    KOKKOS_INLINE_FUNCTION
    void operator()(const LO j) const {
      SC one = ATS::one();
      for  (LO i = 0; i < numPDEs; i++)
        nullspace(j*numPDEs + i, i) = one;
    }
  };
  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class DeviceType>
  void NullspaceFactory_kokkos<Scalar,LocalOrdinal,GlobalOrdinal,Kokkos::Compat::KokkosDeviceWrapperNode<DeviceType>>::Build(Level& currentLevel) const {
    FactoryMonitor m(*this, "Nullspace factory", currentLevel);

    RCP<MultiVector> nullspace;

    //TEUCHOS_TEST_FOR_EXCEPTION(currentLevel.GetLevelID() != 0, Exceptions::RuntimeError, "MueLu::NullspaceFactory::Build(): NullspaceFactory can be used for finest level (LevelID == 0) only.");
    const ParameterList& pL = GetParameterList();
    std::string nspName = pL.get<std::string>("Fine level nullspace");

    if (currentLevel.GetLevelID() == 0) {
      if (currentLevel.IsAvailable(nspName, NoFactory::get())) {
        // When a fine nullspace have already been defined by user using Set("Nullspace", ...) or
        // Set("Nullspace1", ...), we use it.
        nullspace = currentLevel.Get< RCP<MultiVector> >(nspName, NoFactory::get());
        GetOStream(Runtime1) << "Use user-given nullspace " << nspName << ":"
            << " nullspace dimension=" << nullspace->getNumVectors()
            << " nullspace length="    << nullspace->getGlobalLength() << std::endl;

      } else {
        // "Nullspace" (nspName) is not available
        auto A = Get< RCP<Matrix> >(currentLevel, "A");

        // Determine numPDEs
        LO numPDEs = 1;
        if (A->IsView("stridedMaps") == true) {
          Xpetra::viewLabel_t oldView = A->SwitchToView("stridedMaps"); // note: "stridedMaps are always non-overlapping (correspond to range and domain maps!)
          TEUCHOS_TEST_FOR_EXCEPTION(rcp_dynamic_cast<const StridedMap>(A->getRowMap()).is_null(), Exceptions::BadCast,
                                     "MueLu::CoalesceFactory::Build: cast to strided row map failed.");
          numPDEs = rcp_dynamic_cast<const StridedMap>(A->getRowMap())->getFixedBlockSize();
          oldView = A->SwitchToView(oldView);
        }

        GetOStream(Runtime1) << "Generating canonical nullspace: dimension = " << numPDEs << std::endl;

        nullspace = MultiVectorFactory::Build(A->getDomainMap(), numPDEs);
        auto nullspaceView = nullspace->getDeviceLocalView(Xpetra::Access::OverwriteAll);

        int numBlocks = nullspace->getLocalLength() / numPDEs;

        NullspaceFunctor<decltype(nullspaceView), LO> nullspaceFunctor(nullspaceView, numPDEs);
        Kokkos::parallel_for("MueLu:NullspaceF:Build:for", range_type(0,numBlocks), nullspaceFunctor);
        // TODO extend null space factory for blocked multi vectors
      }

    } else {
      // On coarser levels always use "Nullspace" as variable name, since it is expected by
      // tentative P factory to be "Nullspace"
      nullspace = currentLevel.Get< RCP<MultiVector> >("Nullspace", GetFactory(nspName).get()); // NOTE: "Nullspace" and nspName mismatch possible here
    }

    // provide "Nullspace" variable on current level (used by TentativePFactory)
    Set(currentLevel, "Nullspace", nullspace);

  } // Build

} //namespace MueLu

#endif // MUELU_NULLSPACEFACTORY_KOKKOS_DEF_HPP
