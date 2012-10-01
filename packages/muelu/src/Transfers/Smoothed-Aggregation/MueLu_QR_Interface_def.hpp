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
//                    Jeremie Gaidamour (jngaida@sandia.gov)
//                    Jonathan Hu       (jhu@sandia.gov)
//                    Ray Tuminaro      (rstumin@sandia.gov)
//
// ***********************************************************************
//
// @HEADER
#ifndef MUELU_QR_INTERFACE_DEF_HPP
#define MUELU_QR_INTERFACE_DEF_HPP

#include "MueLu_QR_Interface_decl.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

  //! Non-member templated function to handle extracting Q from QR factorization for different Scalar types.
  template <class Scalar, class LocalOrdinal>
  void LapackQR(Teuchos::LAPACK<LocalOrdinal,Scalar> &lapack, LocalOrdinal myAggSize,
                int intFineNSDim, ArrayRCP<Scalar> &localQR, ArrayRCP<Scalar> &tau,
                ArrayRCP<Scalar> &work, LocalOrdinal &workSize, LocalOrdinal &info)
  {
    lapack.ORGQR(myAggSize, intFineNSDim, intFineNSDim, localQR.getRawPtr(),
                 myAggSize, tau.getRawPtr(), work.getRawPtr(), workSize, &info );
  }

  //! Non-member specialized function to handle extracting Q from QR factorization for Scalar==complex
  template <class LocalOrdinal>
  void LapackQR(Teuchos::LAPACK<LocalOrdinal, std::complex<double> > &lapack,
                LocalOrdinal myAggSize, int intFineNSDim, ArrayRCP<std::complex<double> > &localQR,
                ArrayRCP<std::complex<double> > &tau, ArrayRCP<std::complex<double> > &work,
                LocalOrdinal &workSize, LocalOrdinal &info)
  {
    lapack.UNGQR(myAggSize, intFineNSDim, intFineNSDim, localQR.getRawPtr(),
                 myAggSize, tau.getRawPtr(), work.getRawPtr(), workSize, &info );
  }

  template <class Scalar, class LocalOrdinal>
  QR_Interface<Scalar,LocalOrdinal>::QR_Interface(const size_t NSDim) : workSize_(NSDim), NSDim_(NSDim), info_(0) {
    tau_ = ArrayRCP<Scalar>(NSDim);
    work_ = ArrayRCP<Scalar>(NSDim);
  }

  template <class Scalar, class LocalOrdinal>
  void QR_Interface<Scalar,LocalOrdinal>::Compute(LocalOrdinal const &myAggSize, ArrayRCP<Scalar> &localQR)
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
    if (NSDim_ == 1) {
      //only one nullspace vector, so normalize by hand
      Magnitude dtemp=0;
      //scalar type might be complex, so take absolute value.
      for (LocalOrdinal k=0; k<myAggSize; ++k) {dtemp += Teuchos::ScalarTraits<Scalar>::magnitude(localQR[k])*Teuchos::ScalarTraits<Scalar>::magnitude(localQR[k]);}
      dtemp = Teuchos::ScalarTraits<Magnitude>::squareroot(dtemp);
      tau_[0] = localQR[0];
      localQR[0] = dtemp;
    } else {
      lapack_.GEQRF( myAggSize, Teuchos::as<int>(NSDim_), localQR.getRawPtr(), myAggSize,
                    tau_.getRawPtr(), work_.getRawPtr(), workSize_, &info_ );
      if (info_ != 0) {
        std::string msg = "QR_Interface: dgeqrf (LAPACK QR routine) returned error code " + Teuchos::toString(info_);
        throw(Exceptions::RuntimeError(msg));
      }
      // LAPACK may have determined a better length for the work array.  Returns it in work[0],
      // so we cast to avoid compiler warnings.  Taking a look at the NETLIB reference implementation
      // CGEQRF (complex), work[0] is assigned an integer, so it's safe to take the magnitude.
      // Scalar type might be complex, so take absolute value.
      if ( Teuchos::ScalarTraits<Scalar>::magnitude(work_[0]) > workSize_) {
        workSize_ = Teuchos::as<int>(Teuchos::ScalarTraits<Scalar>::magnitude(work_[0]));
        work_ = ArrayRCP<Scalar>(workSize_);
      }
    }
  } //Compute()

  template <class Scalar, class LocalOrdinal>
  void QR_Interface<Scalar,LocalOrdinal>::ExtractQ(LocalOrdinal const &myAggSize, ArrayRCP<Scalar> &localQR)
  {
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType Magnitude;
    if (NSDim_ == 1) {
      //again, only one nullspace vector, so calculate Q by hand
      Magnitude dtemp = Teuchos::ScalarTraits<Scalar>::magnitude(localQR[0]);
      localQR[0] = tau_[0];
      dtemp = 1 / dtemp;
      for (LocalOrdinal i=0; i<myAggSize; ++i)
        localQR[i] *= dtemp;
    } else {
      //lapack_.ORGQR(myAggSize, Teuchos::as<int>(workSize_), Teuchos::as<int>(workSize_), localQR.getRawPtr(),
      //             myAggSize, tau_.getRawPtr(), work_.getRawPtr(), workSize_, &info_ );
      //call nonmember function (perhaps specialized)
      LapackQR( lapack_, myAggSize, Teuchos::as<int>(NSDim_), localQR, tau_, work_, workSize_, info_ );
      if (info_ != 0) {
        std::string msg = "QR_Interface: dorgqr (LAPACK auxiliary QR routine) returned error code " + Teuchos::toString(info_);
        throw(Exceptions::RuntimeError(msg));
      }

      // LAPACK may have determined a better length for the work array.  Returns it in work[0],
      // so we cast to avoid compiler warnings.
      // Scalar type might be complex, so take absolute value.
      if ( Teuchos::ScalarTraits<Scalar>::magnitude(work_[0]) > workSize_) {
        workSize_ = Teuchos::as<int>(Teuchos::ScalarTraits<Scalar>::magnitude(work_[0]));
        work_ = ArrayRCP<Scalar>(workSize_);
      }
    }
  } //ExtractQ()

} //namespace MueLu

#endif // MUELU_QR_INTERFACE_DEF_HPP
