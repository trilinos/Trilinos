// @HEADER
// ***********************************************************************
//
//         Stratimikos: Thyra-based strategies for linear solvers
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
// Questions? Contact Jennifer A. Loe (jloe@sandia.gov)
//
// ***********************************************************************
// @HEADER
#ifndef STRATIMIKOS_BELOS_PREC_TPETRA_HELPERS_HPP
#define STRATIMIKOS_BELOS_PREC_TPETRA_HELPERS_HPP

#include "Stratimikos_LinearSolverBuilder.hpp"
#include "Thyra_BelosTpetraPreconditionerFactory_decl.hpp"
#include "Thyra_BelosTpetraPreconditionerFactory_def.hpp"

#include "Teuchos_RCP.hpp"
#include "Teuchos_ParameterList.hpp"
#include "Teuchos_TestForException.hpp"
#include "Teuchos_AbstractFactoryStd.hpp"

#include <string>

namespace Stratimikos {

  template <typename MatrixType>
  void enableBelosPrecTpetra(LinearSolverBuilder<typename MatrixType::scalar_type>& builder, const std::string& stratName = "BelosPrecTpetra")
  {
    const Teuchos::RCP<const Teuchos::ParameterList> precValidParams = Teuchos::sublist(builder.getValidParameters(), "Preconditioner Types");

    TEUCHOS_TEST_FOR_EXCEPTION(precValidParams->isParameter(stratName), std::logic_error,
                               "Stratimikos::enableBelosPrecTpetra cannot add \"" + stratName +"\" because it is already included in builder!");

    typedef typename MatrixType::scalar_type scalar_type;
    typedef Thyra::PreconditionerFactoryBase<scalar_type>       Base;
    typedef Thyra::BelosTpetraPreconditionerFactory<MatrixType> Impl;

    builder.setPreconditioningStrategyFactory(Teuchos::abstractFactoryStd<Base, Impl>(), stratName);
  }

} // namespace Stratimikos

#endif
