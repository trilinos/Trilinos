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
#ifndef MUELU_TRANSPFACTORY_DEF_HPP
#define MUELU_TRANSPFACTORY_DEF_HPP

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_Time.hpp>

#include <Xpetra_Matrix.hpp>

#include "MueLu_TransPFactory_decl.hpp"

#include "MueLu_Monitor.hpp"
#include "MueLu_PerfUtils.hpp"
#include "MueLu_Utilities.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
RCP<const ParameterList> TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::GetValidParameterList() const {
  RCP<ParameterList> validParamList = rcp(new ParameterList());
  validParamList->set<RCP<const FactoryBase> >("P", Teuchos::null, "Generating factory of the matrix P");

  // Make sure we don't recursively validate options for the matrixmatrix kernels
  ParameterList norecurse;
  norecurse.disableRecursiveValidation();
  validParamList->set<ParameterList>("matrixmatrix: kernel params", norecurse, "MatrixMatrix kernel parameters");

  return validParamList;
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::DeclareInput(Level& /* fineLevel */, Level& coarseLevel) const {
  Input(coarseLevel, "P");
}

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void TransPFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node>::Build(Level& /* fineLevel */, Level& coarseLevel) const {
  FactoryMonitor m(*this, "Transpose P", coarseLevel);
  std::string label = "MueLu::TransP-" + Teuchos::toString(coarseLevel.GetLevelID());

  RCP<Matrix> P = Get<RCP<Matrix> >(coarseLevel, "P");
  // If we failed to create a valid P (e.g., # of global aggregates is zero), then we just bail here
  //  This level will ultimately be removed in MueLu_Hierarchy_defs.h via a resize()
  if (P == Teuchos::null) return;

  const Teuchos::ParameterList& pL = GetParameterList();

  // Reuse pattern if available (multiple solve)
  RCP<ParameterList> Tparams;
  if (pL.isSublist("matrixmatrix: kernel params"))
    Tparams = rcp(new ParameterList(pL.sublist("matrixmatrix: kernel params")));
  else
    Tparams = rcp(new ParameterList);

  // By default, we don't need global constants for transpose
  Tparams->set("compute global constants: temporaries", Tparams->get("compute global constants: temporaries", false));
  Tparams->set("compute global constants", Tparams->get("compute global constants", false));

  RCP<Matrix> R = Utilities::Transpose(*P, true, label, Tparams);

  if (IsPrint(Statistics2)) {
    RCP<ParameterList> params = rcp(new ParameterList());
    params->set("printLoadBalancingInfo", true);
    params->set("printCommInfo", true);
    GetOStream(Statistics2) << PerfUtils::PrintMatrixInfo(*R, "R", params);
  }

  Set(coarseLevel, "R", R);

  ///////////////////////// EXPERIMENTAL
  if (P->IsView("stridedMaps"))
    R->CreateView("stridedMaps", P, true);
  ///////////////////////// EXPERIMENTAL
}

}  // namespace MueLu

#endif  // MUELU_TRANSPFACTORY_DEF_HPP
