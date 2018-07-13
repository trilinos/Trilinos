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
//                    Tobias Wiesner    (tawiesn@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef STRATIMIKOS_MUELU_TPETRA_HELPERS_HPP
#define STRATIMIKOS_MUELU_TPETRA_HELPERS_HPP

#include "Stratimikos_DefaultLinearSolverBuilder.hpp"

#include "Thyra_MueLuPreconditionerFactory.hpp"
#include "Thyra_MueLuRefMaxwellPreconditionerFactory.hpp"

#if defined(HAVE_MUELU_EXPERIMENTAL) && defined(HAVE_MUELU_TEKO)
#include "Thyra_MueLuTpetraQ2Q1PreconditionerFactory.hpp"
#endif

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include <string>

namespace Stratimikos {

  // Dynamically register MueLu Tpetra adapters in Stratimikos
  // Note: No Scalar template argument is available because Stratimikos
  // does not support types beyond double
  template <typename LocalOrdinal = int, typename GlobalOrdinal = int, typename Node = KokkosClassic::DefaultNode::DefaultNodeType>
  void enableMueLu(DefaultLinearSolverBuilder& builder, const std::string& stratName = "MueLu")
  {
    const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

    TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                               "Stratimikos::enableMueLu cannot add \"" + stratName +"\" because it is already included in builder!");

    typedef Thyra::PreconditionerFactoryBase<double>                                     Base;
    typedef Thyra::MueLuPreconditionerFactory<double, LocalOrdinal, GlobalOrdinal, Node> Impl;

    builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
  }

  // Dynamically register MueLu RefMaxwell adapters in Stratimikos
  // Note: No Scalar template argument is available because Stratimikos
  // does not support types beyond double
  template <typename LocalOrdinal = int, typename GlobalOrdinal = int, typename Node = KokkosClassic::DefaultNode::DefaultNodeType>
  void enableMueLuRefMaxwell(DefaultLinearSolverBuilder& builder, const std::string& stratName = "MueLuRefMaxwell")
  {
    const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

    TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                               "Stratimikos::enableMueLuRefMaxwell cannot add \"" + stratName +"\" because it is already included in builder!");

    typedef Thyra::PreconditionerFactoryBase<double>                                     Base;
    typedef Thyra::MueLuRefMaxwellPreconditionerFactory<double, LocalOrdinal, GlobalOrdinal, Node> Impl;

    builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
  }

#if defined(HAVE_MUELU_EXPERIMENTAL) && defined(HAVE_MUELU_TEKO)
#if 0
  // Dynamically register MueLu Tpetra adapters in Stratimikos
  void enableMueLuTpetraQ2Q1(DefaultLinearSolverBuilder &builder, const std::string &stratName = "MueLu");
#endif

  // Dynamically register MueLu Tpetra adapters in Stratimikos
  // Note: No Scalar template argument is available because Stratimikos
  // does not support types beyond double
  template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
  void enableMueLuTpetraQ2Q1(DefaultLinearSolverBuilder& builder, const std::string &stratName = "MueLu") {
    const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

    TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                               "Stratimikos::enableMueLuTpetraQ2Q1 cannot add \"" + stratName +"\" because it is already included in builder!");

    typedef Thyra::PreconditionerFactoryBase<double>                                              Base;
    typedef Thyra::MueLuTpetraQ2Q1PreconditionerFactory<double, LocalOrdinal, GlobalOrdinal,Node> Impl;

    builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
  }
#endif

} // namespace Stratimikos

#endif
