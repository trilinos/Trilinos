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
#ifndef MUELU_SMOOTHERBASECLONER_HPP
#define MUELU_SMOOTHERBASECLONER_HPP

#include <Xpetra_MultiVector_fwd.hpp>

#include "MueLu_ConfigDefs.hpp"
#include "MueLu_BaseClass.hpp"

#include "MueLu_SmootherBase.hpp"
#if defined(HAVE_MUELU_IFPACK2)
#include "MueLu_Ifpack2Smoother.hpp"
#endif
#include "MueLu_TrilinosSmoother.hpp"

namespace MueLu {

template <class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node1, class Node2>
Teuchos::RCP<SmootherBase<Scalar,LocalOrdinal,GlobalOrdinal,Node2> >
clone (const Teuchos::RCP<SmootherBase<Scalar,LocalOrdinal,GlobalOrdinal,Node1> >& SB,
       const Teuchos::RCP<const Xpetra::Matrix<Scalar,LocalOrdinal,GlobalOrdinal,Node2> >& cloneA,
       const RCP<Node2>& node2)
{
#if defined(HAVE_MUELU_IFPACK2)
  Teuchos::RCP<SmootherBase<Scalar, LocalOrdinal, GlobalOrdinal, Node2> >  cloneSB;
  Teuchos::RCP<TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node1> > trilSmoother =
    Teuchos::rcp_dynamic_cast<TrilinosSmoother<Scalar, LocalOrdinal, GlobalOrdinal, Node1> >(SB);
  if (trilSmoother != Teuchos::null){
    cloneSB = trilSmoother->clone (node2, cloneA);
    return cloneSB;
  }

  Teuchos::RCP<Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node1> > ifSmoother =
    Teuchos::rcp_dynamic_cast<Ifpack2Smoother<Scalar, LocalOrdinal, GlobalOrdinal, Node1> >(SB);
  if (ifSmoother != Teuchos::null){
    cloneSB = ifSmoother->clone(node2, cloneA);
    return cloneSB;
  }

  else {
    TEUCHOS_TEST_FOR_EXCEPTION(
        true, std::invalid_argument, "MueLu::SmootherClone: "
        "Invalid smoother type to clone (not TrilinosSmoother or Ifpack2 ) \"");
  }
#else
  TEUCHOS_TEST_FOR_EXCEPTION(
      true, std::invalid_argument, "MueLu::SmootherClone: "
      "clone() only available with IFPACK2 enabled.");
#endif

}
} //namespace MueLu

#define MUELU_SMOOTHERBASECLONER_SHORT
#endif //ifndef MUELU_SMOOTHERBASECLONER_HPP

