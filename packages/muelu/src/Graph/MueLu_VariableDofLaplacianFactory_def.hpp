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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef PACKAGES_MUELU_SRC_GRAPH_MUELU_VARIABLEDOFLAPLACIANFACTORY_DEF_HPP_
#define PACKAGES_MUELU_SRC_GRAPH_MUELU_VARIABLEDOFLAPLACIANFACTORY_DEF_HPP_


#include "MueLu_Monitor.hpp"

#include "MueLu_VariableDofLaplacianFactory_decl.hpp"

namespace MueLu {

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());
/*
#define SET_VALID_ENTRY(name) validParamList->setEntry(name, MasterList::getEntry(name))
    SET_VALID_ENTRY("aggregation: drop tol");
    SET_VALID_ENTRY("aggregation: Dirichlet threshold");
    SET_VALID_ENTRY("aggregation: drop scheme");
    {
      typedef Teuchos::StringToIntegralParameterEntryValidator<int> validatorType;
      validParamList->getEntry("aggregation: drop scheme").setValidator(
        rcp(new validatorType(Teuchos::tuple<std::string>("classical", "distance laplacian"), "aggregation: drop scheme")));
    }
#undef  SET_VALID_ENTRY
    validParamList->set< bool >                  ("lightweight wrap",           true, "Experimental option for lightweight graph access");
*/

    validParamList->set< RCP<const FactoryBase> >("A",                  Teuchos::null, "Generating factory of the matrix A");
    validParamList->set< RCP<const FactoryBase> >("Coordinates",        Teuchos::null, "Generating factory for Coordinates");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::VariableDofLaplacianFactory() { }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    //const ParameterList& pL = GetParameterList();
    Input(currentLevel, "A");
    Input(currentLevel, "Coordinates");

    //if (currentLevel.GetLevelID() == 0) // TODO check for finest level (special treatment)
    if (currentLevel.IsAvailable("DofPresent", NoFactory::get())) {
      currentLevel.DeclareInput("DofPresent", NoFactory::get(), this);
    }
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level &currentLevel) const {
    FactoryMonitor m(*this, "Build", currentLevel);
    typedef Teuchos::ScalarTraits<SC> STS;

    const ParameterList  & pL = GetParameterList();

    RCP<Matrix> A = Get< RCP<Matrix> >(currentLevel, "A");

    typedef Xpetra::MultiVector<double,LO,GO,NO> dxMV;
    RCP<dxMV> Coords = Get< RCP<Xpetra::MultiVector<double,LO,GO,NO> > >(currentLevel, "Coordinates");

    int maxDofPerNode = 3;    // TODO add parameter from parameter list for this...
    Scalar dirDropTol = 1e-5; // TODO parameter from parameter list ("ML advnaced Dirichlet: threshold")

    bool bHasZeroDiagonal = false;
    Teuchos::ArrayRCP<const bool> dirOrNot = MueLu::Utilities<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DetectDirichletRowsExt(*A,bHasZeroDiagonal,dirDropTol);

    // TODO check availability + content (length)
    Teuchos::ArrayRCP<const bool> dofPresent = currentLevel.Get< Teuchos::ArrayRCP<const bool> >("DofPresent", NoFactory::get());

    // map[k] indicates that the kth dof in the variable dof matrix A would
    // correspond to the map[k]th dof in the padded system. If, i.e., it is
    // map[35] = 39 then dof no 35 in the variable dof matrix A corresponds to
    // row map id 39 in an imaginary padded matrix Apadded.
    // The padded system is never built but would be the associated matrix if
    // every node had maxDofPerNode dofs.
    std::vector<LocalOrdinal> map(A->getNodeNumRows());
    this->buildPaddedMap(dofPresent, map, A->getNodeNumRows());

    //GetOStream(Parameters0) << "lightweight wrap = " << doExperimentalWrap << std::endl;


    std::vector<LocalOrdinal> myLocalNodeIds(A->getRowMap()->getNodeNumElements()); // we reserve the absolute possible maximum

    // calculate local node ids for local dof ids.
    for(size_t i = 0; i < A->getNodeNumRows(); i++)
      myLocalNodeIds[i] = std::floor<LocalOrdinal>( map[i] / maxDofPerNode );

    // assign the local node ids for the ghosted nodes

  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::buildPaddedMap(const Teuchos::ArrayRCP<const bool> & dofPresent, std::vector<LocalOrdinal> & map, size_t nDofs) const {
    size_t count = 0;
    for (size_t i = 0; i < dofPresent.size(); i++)
      if(dofPresent[i] == true) map[count++] = Teuchos::as<LocalOrdinal>(i);
    TEUCHOS_TEST_FOR_EXCEPTION(nDofs != count, MueLu::Exceptions::RuntimeError, "VariableDofLaplacianFactory::buildPaddedMap: #dofs in dofPresent does not match the expected value (number of rows of A): " << nDofs << " vs. " << count);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void VariableDofLaplacianFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::assignGhostLocalNodeIds(const std::vector<LocalOrdinal> & myLocalNodeIds, size_t nLocalDofs, size_t nLocalPlusGhostDofs, size_t& nLocalNodes, size_t& nLocalPlusGhostNodes) const {

  }

} /* MueLu */


#endif /* PACKAGES_MUELU_SRC_GRAPH_MUELU_VARIABLEDOFLAPLACIANFACTORY_DEF_HPP_ */
