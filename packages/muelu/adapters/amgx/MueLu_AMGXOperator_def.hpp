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

#if defined(HAVE_MUELU_AMGX)
#include "MueLu_AMGXOperator_decl.hpp"

namespace MueLu {

template <class Node>
Teuchos::RCP<const Tpetra::Map<int, int, Node> >
AMGXOperator<double, int, int, Node>::getDomainMap() const {
  return domainMap_;
}

template <class Node>
Teuchos::RCP<const Tpetra::Map<int, int, Node> > AMGXOperator<double, int, int, Node>::getRangeMap() const {
  return rangeMap_;
}

template <class Node>
void AMGXOperator<double, int, int, Node>::apply(const Tpetra::MultiVector<double, int, int, Node>& X,
                                                 Tpetra::MultiVector<double, int, int, Node>& Y,
                                                 Teuchos::ETransp mode, double alpha, double beta) const {
  RCP<const Teuchos::Comm<int> > comm = Y.getMap()->getComm();

  ArrayRCP<const double> mueluXdata, amgxXdata;
  ArrayRCP<double> mueluYdata, amgxYdata;

  try {
    for (int i = 0; i < (int)Y.getNumVectors(); i++) {
      {
        vectorTimer1_->start();

        mueluXdata = X.getData(i);
        mueluYdata = Y.getDataNonConst(i);

        if (comm->getSize() == 1) {
          amgxXdata = mueluXdata;
          amgxYdata = mueluYdata;

        } else {
          int n = mueluXdata.size();

          amgxXdata.resize(n);
          amgxYdata.resize(n);

          ArrayRCP<double> amgxXdata_nonConst = Teuchos::arcp_const_cast<double>(amgxXdata);
          for (int j = 0; j < n; j++) {
            amgxXdata_nonConst[muelu2amgx_[j]] = mueluXdata[j];
            amgxYdata[muelu2amgx_[j]]          = mueluYdata[j];
          }
        }

        AMGX_vector_upload(X_, N_, 1, &amgxXdata[0]);
        AMGX_vector_upload(Y_, N_, 1, &amgxYdata[0]);

        vectorTimer1_->stop();
        vectorTimer1_->incrementNumCalls();
      }

      // Solve the system and time.
      solverTimer_->start();
      AMGX_solver_solve(Solver_, X_, Y_);
      solverTimer_->stop();
      solverTimer_->incrementNumCalls();

      {
        vectorTimer2_->start();

        AMGX_vector_download(Y_, &amgxYdata[0]);

        if (comm->getSize() > 1) {
          int n = mueluYdata.size();

          for (int j = 0; j < n; j++)
            mueluYdata[j] = amgxYdata[muelu2amgx_[j]];
        }

        vectorTimer2_->stop();
        vectorTimer2_->incrementNumCalls();
      }
    }

  } catch (std::exception& e) {
    std::string errMsg = std::string("Caught an exception in MueLu::AMGXOperator::Apply():\n") + e.what() + "\n";
    throw Exceptions::RuntimeError(errMsg);
  }
}

template <class Node>
bool AMGXOperator<double, int, int, Node>::hasTransposeApply() const {
  return false;
}

}  // namespace MueLu
#endif  // if defined(HAVE_MUELU_AMGX)

#endif  // ifdef MUELU_AMGXOPERATOR_DEF_HPP
