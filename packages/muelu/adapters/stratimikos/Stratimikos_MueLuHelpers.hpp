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

#include "MueLu_Details_DefaultTypes.hpp"
#include "Stratimikos_LinearSolverBuilder.hpp"

#include "Thyra_MueLuPreconditionerFactory.hpp"
#include "Thyra_MueLuRefMaxwellPreconditionerFactory.hpp"
#include "Thyra_MueLuMaxwell1PreconditionerFactory.hpp"

#if defined(HAVE_MUELU_EXPERIMENTAL) && defined(HAVE_MUELU_TEKO)
#include "Thyra_MueLuTpetraQ2Q1PreconditionerFactory.hpp"
#endif

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include <string>

namespace Stratimikos {

template <typename Scalar = MueLu::DefaultScalar, typename LocalOrdinal = MueLu::DefaultLocalOrdinal, typename GlobalOrdinal = MueLu::DefaultGlobalOrdinal, typename Node = MueLu::DefaultNode>
void enableMueLu(LinearSolverBuilder<Scalar>& builder, const std::string& stratName = "MueLu") {
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
  const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

  TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                             "Stratimikos::enableMueLu cannot add \"" + stratName + "\" because it is already included in builder!");

  typedef Thyra::PreconditionerFactoryBase<Scalar> Base;
  typedef Thyra::MueLuPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> Impl;

  builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
#endif
}

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
MUELU_DEPRECATED void enableMueLu(LinearSolverBuilder<double>& builder, const std::string& stratName = "MueLu") {
  enableMueLu<double, LocalOrdinal, GlobalOrdinal, Node>(builder, stratName);
}

template <typename Scalar = MueLu::DefaultScalar, typename LocalOrdinal = MueLu::DefaultLocalOrdinal, typename GlobalOrdinal = MueLu::DefaultGlobalOrdinal, typename Node = MueLu::DefaultNode>
void enableMueLuRefMaxwell(LinearSolverBuilder<Scalar>& builder, const std::string& stratName = "MueLuRefMaxwell") {
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
  const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

  TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                             "Stratimikos::enableMueLuRefMaxwell cannot add \"" + stratName + "\" because it is already included in builder!");

  typedef Thyra::PreconditionerFactoryBase<Scalar> Base;
  typedef Thyra::MueLuRefMaxwellPreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> Impl;

  builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
#endif
}

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
MUELU_DEPRECATED void enableMueLuRefMaxwell(LinearSolverBuilder<double>& builder, const std::string& stratName = "MueLuRefMaxwell") {
  enableMueLuRefMaxwell<double, LocalOrdinal, GlobalOrdinal, Node>(builder, stratName);
}

template <typename Scalar = MueLu::DefaultScalar, typename LocalOrdinal = MueLu::DefaultLocalOrdinal, typename GlobalOrdinal = MueLu::DefaultGlobalOrdinal, typename Node = MueLu::DefaultNode>
void enableMueLuMaxwell1(LinearSolverBuilder<Scalar>& builder, const std::string& stratName = "MueLuMaxwell1") {
#if defined(HAVE_MUELU_STRATIMIKOS) && defined(HAVE_MUELU_THYRA)
  const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

  TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                             "Stratimikos::enableMueLuRefMaxwell cannot add \"" + stratName + "\" because it is already included in builder!");

  typedef Thyra::PreconditionerFactoryBase<Scalar> Base;
  typedef Thyra::MueLuMaxwell1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> Impl;

  builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
#endif
}

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
MUELU_DEPRECATED void enableMueLuMaxwell1(LinearSolverBuilder<double>& builder, const std::string& stratName = "MueLuMaxwell1") {
  enableMueLuMaxwell1<double, LocalOrdinal, GlobalOrdinal, Node>(builder, stratName);
}

#if defined(HAVE_MUELU_EXPERIMENTAL) && defined(HAVE_MUELU_TEKO)
#if 0
  // Dynamically register MueLu Tpetra adapters in Stratimikos
  void enableMueLuTpetraQ2Q1(DefaultLinearSolverBuilder &builder, const std::string &stratName = "MueLu");
#endif

template <typename Scalar = MueLu::DefaultScalar, typename LocalOrdinal = MueLu::DefaultLocalOrdinal, typename GlobalOrdinal = MueLu::DefaultGlobalOrdinal, typename Node = MueLu::DefaultNode>
void enableMueLuTpetraQ2Q1(LinearSolverBuilder<Scalar>& builder, const std::string& stratName = "MueLu") {
  const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

  TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                             "Stratimikos::enableMueLuTpetraQ2Q1 cannot add \"" + stratName + "\" because it is already included in builder!");

  typedef Thyra::PreconditionerFactoryBase<Scalar> Base;
  typedef Thyra::MueLuTpetraQ2Q1PreconditionerFactory<Scalar, LocalOrdinal, GlobalOrdinal, Node> Impl;

  builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
}

template <typename LocalOrdinal, typename GlobalOrdinal, typename Node>
MUELU_DEPRECATED void enableMueLuTpetraQ2Q1(LinearSolverBuilder<double>& builder, const std::string& stratName = "MueLu") {
  enableMueLuTpetraQ2Q1<double, LocalOrdinal, GlobalOrdinal, Node>(builder, stratName);
}
#endif

}  // namespace Stratimikos

#endif
