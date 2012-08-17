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
#ifndef STOKHOS_MUELU_QR_INTERFACE_DEF_HPP
#define STOKHOS_MUELU_QR_INTERFACE_DEF_HPP

#include "MueLu_QR_Interface_decl.hpp"
#include "MueLu_Exceptions.hpp"

namespace MueLu {

#if defined(HAVE_STOKHOS_MUELU)
  //Specialization for polynomial chaos expansion (PCE) scalar types.
  template <class Scalar, class Storage, class LocalOrdinal>
  QR_Interface< Sacado::PCE::OrthogPoly<Scalar, Storage>, LocalOrdinal>::QR_Interface(const size_t NSDim) : workSize_(NSDim), NSDim_(NSDim), info_(0) {
        tau_ = ArrayRCP<Scalar>(NSDim);
        work_ = ArrayRCP<Scalar>(NSDim);
      }

  template <class Scalar, class Storage, class LocalOrdinal>
  void QR_Interface< Sacado::PCE::OrthogPoly<Scalar, Storage>, LocalOrdinal>::Compute(LocalOrdinal const &myAggSize, ArrayRCP<Sacado::PCE::OrthogPoly<Scalar, Storage> > &localQR) {
        if (localQR.size() > localQR_.size())
          localQR_.resize(localQR.size());
        //convert pce to pod scalar
        for (int i=0; i<localQR.size(); ++i) {
          localQR_[i] = (localQR[i]).coeff(0);
        }
        lapack_.GEQRF( myAggSize, Teuchos::as<int>(NSDim_), localQR_.getRawPtr(), myAggSize,
                      tau_.getRawPtr(), work_.getRawPtr(), workSize_, &info_ );
        if (info_ != 0) {
          std::string msg = "QR_Interface: dgeqrf (LAPACK QR routine) returned error code " + Teuchos::toString(info_);
          throw(Exceptions::RuntimeError(msg));
        }
        //promote POD scalar back to pce
        for (int i=0; i<localQR.size(); ++i) {
          localQR[i] = localQR_[i];
        }
        // LAPACK may have determined a better length for the work array.  Returns it in work[0],
        // so we cast to avoid compiler warnings.  Taking a look at the NETLIB reference implementation
        // CGEQRF (complex), work[0] is assigned an integer, so it's safe to take the magnitude.
        // Scalar type might be complex, so take absolute value.
        if ( std::abs(work_[0]) > workSize_) {
          workSize_ = (int) std::abs(work_[0]);
          work_ = ArrayRCP<Scalar>(workSize_);
        }
      } //Compute

  template <class Scalar, class Storage, class LocalOrdinal>
  void QR_Interface< Sacado::PCE::OrthogPoly<Scalar, Storage>, LocalOrdinal>::ExtractQ(LocalOrdinal const &myAggSize, ArrayRCP<Sacado::PCE::OrthogPoly<Scalar, Storage> > &localQR) {
        //call nonmember function (perhaps specialized)
        //Note: localQR_ already contains the proper data because of prior call to Compute, so there is no need to resize or copy.
        //      If Compute is called twice in a row, all bets are off.
        LapackQR( lapack_, myAggSize, Teuchos::as<int>(NSDim_), localQR_, tau_, work_, workSize_, info_ );
        if (info_ != 0) {
          std::string msg = "QR_Interface: dorgqr (LAPACK auxiliary QR routine) returned error code " + Teuchos::toString(info_);
          throw(Exceptions::RuntimeError(msg));
        }
        //promote POD scalar back to pce
        for (int i=0; i<localQR.size(); ++i) {
          localQR[i] = localQR_[i];
        }
  
        // LAPACK may have determined a better length for the work array.  Returns it in work[0],
        // so we cast to avoid compiler warnings.
        // Scalar type might be complex, so take absolute value.
        if ( std::abs(work_[0]) > workSize_) {
          workSize_ = (int) std::abs(work_[0]);
          work_ = ArrayRCP<Scalar>(workSize_);
        }
      } //ExtractQ
#endif //if defined(HAVE_STOKHOS_MUELU)


} //namespace MueLu

#endif // STOKHOS_MUELU_QR_INTERFACE_DEF_HPP
