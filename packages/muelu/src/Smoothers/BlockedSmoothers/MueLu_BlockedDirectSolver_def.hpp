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
/*
 * MueLu_BlockedDirectSolver_def.hpp
 *
 *  Created on: 09.02.2014
 *      Author: tobias
 */

#ifndef MUELU_BLOCKEDDIRECTSOLVER_DEF_HPP_
#define MUELU_BLOCKEDDIRECTSOLVER_DEF_HPP_

#include "Teuchos_ArrayViewDecl.hpp"
#include "Teuchos_ScalarTraits.hpp"

#include "MueLu_ConfigDefs.hpp"

#include <Xpetra_Matrix.hpp>
#include <Xpetra_BlockedCrsMatrix.hpp>
#include <Xpetra_MultiVectorFactory.hpp>

#include "MueLu_BlockedDirectSolver.hpp"
#include "MueLu_MergedBlockedMatrixFactory.hpp"
#include "MueLu_Level.hpp"
#include "MueLu_Utilities.hpp"
#include "MueLu_Monitor.hpp"
#include "MueLu_HierarchyHelpers.hpp"
#include "MueLu_SmootherBase.hpp"
#include "MueLu_DirectSolver.hpp"

namespace MueLu {

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::BlockedDirectSolver()
    : type_("blocked direct solver")
  {
    MergedAFact_ = Teuchos::rcp(new MergedBlockedMatrixFactory());

    Teuchos::ParameterList params;
    s_ = Teuchos::rcp(new DirectSolver("Klu", params));
    //s_->SetFactory("A", this); // use this factory as generating factory of the merged matrix A
  }

  template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<const ParameterList> BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
    RCP<ParameterList> validParamList = rcp(new ParameterList());

    validParamList->set< RCP<const FactoryBase> >("A", null, "Generating factory of the matrix A");

    return validParamList;
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level &currentLevel) const {
    //this->Input(currentLevel, "A");
    // TODO: check me: why is this->Input not freeing properly A in release mode?
    currentLevel.DeclareInput("A",this->GetFactory("A").get());

    MergedAFact_->SetFactory("A", this->GetFactory("A"));
    //MergedAFact_->DeclareInput(currentLevel);

    s_->SetFactory("A",MergedAFact_);
    s_->DeclareInput(currentLevel);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Setup(Level &currentLevel) {
    RCP<Teuchos::FancyOStream> out = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));

    FactoryMonitor m(*this, "Setup BlockedDirectSolver", currentLevel);
    if (this->IsSetup() == true)
      this->GetOStream(Warnings0) << "MueLu::BlockedDirectSolver::Setup(): Setup() has already been called";

    // extract blocked operator A from current level
    A_ = Factory::Get< RCP<Matrix> >(currentLevel, "A"); // A needed for extracting map extractors
    RCP<BlockedCrsMatrix> bA = Teuchos::rcp_dynamic_cast<BlockedCrsMatrix>(A_);
    TEUCHOS_TEST_FOR_EXCEPTION(bA.is_null(), Exceptions::BadCast,
                               "MueLu::BlockedDirectSolver::Build: input matrix A is not of type BlockedCrsMatrix.");

    s_->Setup(currentLevel);

    this->IsSetup(true);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Apply(MultiVector &X, const MultiVector& B, bool InitialGuessIsZero) const {
    TEUCHOS_TEST_FOR_EXCEPTION(this->IsSetup() == false, Exceptions::RuntimeError,
                               "MueLu::BlockedDirectSolver::Apply(): Setup() has not been called");

    s_->Apply(X, B, InitialGuessIsZero);
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  RCP<MueLu::SmootherPrototype<Scalar, LocalOrdinal, GlobalOrdinal, Node> >
  BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Copy() const {
    return rcp (new BlockedDirectSolver (*this));
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  std::string BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::description() const {
    std::ostringstream out;
    out << SmootherPrototype::description();
    out << "{type = " << type_ << "}";
    return out.str();
  }

  template <class Scalar,class LocalOrdinal, class GlobalOrdinal, class Node>
  void BlockedDirectSolver<Scalar, LocalOrdinal, GlobalOrdinal, Node>::print(Teuchos::FancyOStream &out, const VerbLevel verbLevel) const {
    MUELU_DESCRIBE;

    if (verbLevel & Parameters0)
      out0 << "Prec. type: " << type_ << std::endl;

    if (verbLevel & Debug)
      out0 << "IsSetup: " << Teuchos::toString(SmootherPrototype::IsSetup()) << std::endl;
  }

} // namespace MueLu



#endif /* MUELU_BLOCKEDDIRECTSOLVER_DEF_HPP_ */
