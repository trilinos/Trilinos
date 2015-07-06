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

#ifndef MUELU_AMGXOPERATOR_DEF_HPP
#define MUELU_AMGXOPERATOR_DEF_HPP


#if defined (HAVE_MUELU_EXPERIMENTAL) and defined (HAVE_MUELU_AMGX)
#include "MueLu_AMGXOperator_decl.hpp"
#include <amgx_c.h>
#include <Teuchos_ArrayRCP.hpp>
#include "cuda_runtime.h"

namespace MueLu {

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> >
AMGXOperator::getDomainMap() const {
   return domainMap_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
Teuchos::RCP<const Tpetra::Map<LocalOrdinal,GlobalOrdinal,Node> > AMGXOperator::getRangeMap() const {
   return rangeMap_;
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
void AMGXOperator::apply(const Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& X,
                                                                               Tpetra::MultiVector<Scalar,LocalOrdinal,GlobalOrdinal,Node>& Y,
                                                                               Teuchos::ETransp mode, Scalar alpha, Scalar beta) const {
  try {

    Teuchos::ArrayRCP<double> xdata, ydata;
    for(int i = 0; i < X.getNumVecs(); i++){
	xdata = X.getDataNonConst(i);
	ydata = Y.getData(i); 
        AMGX_vector_upload(X_, N_, 1, &xdata[0]);
   	AMGX_vector_upload(Y_, N_, 1, &ydata[0]);
   	AMGX_solver_solve(Solver_, X_, Y_);
        AMGX_vector_download(yVec, &ydata[0]);
    }
	
    
  } catch (std::exception& e) {
      std::string errMsg = "Caught an exception in MueLu::AMGXOperator::Apply():\n" + e.what() + "\n";
      throw Exceptions::RuntimeError(errMsg);
  }
}

template<class Scalar, class LocalOrdinal, class GlobalOrdinal, class Node>
bool AMGXOperator::hasTransposeApply() const {
  return false;
}

} // namespace
#endif //ifdef HAVE_MUELU_EXPERIMENTAL and defined(HAVE_MUELU_AMGX)

#endif //ifdef MUELU_AMGXOPERATOR_DEF_HPP
