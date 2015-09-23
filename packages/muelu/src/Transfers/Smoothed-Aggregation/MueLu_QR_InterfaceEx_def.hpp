/*
 * MueLu_QR_InterfaceEx_def.hpp
 *
 *  Created on: Jan 31, 2013
 *      Author: tobias
 */

#ifndef MUELU_QR_INTERFACEEX_DEF_HPP_
#define MUELU_QR_INTERFACEEX_DEF_HPP_

#include "MueLu_QR_InterfaceEx_decl.hpp"
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
  QR_InterfaceEx<Scalar,LocalOrdinal>::QR_InterfaceEx(const size_t NSDim) : workSize_(NSDim), NSDim_(NSDim), info_(0) {
    tau_ = ArrayRCP<Scalar>(NSDim);
    work_ = ArrayRCP<Scalar>(NSDim);
  }

  template <class Scalar, class LocalOrdinal>
  bool QR_InterfaceEx<Scalar,LocalOrdinal>::isZeroNspColumn(LocalOrdinal const &myAggSize, ArrayRCP<Scalar> &localQR, LocalOrdinal nspCol) {
    bool bIsZeroColumn = true;
    for(LocalOrdinal r = 0; r < myAggSize; ++r) {
      if(localQR[ myAggSize*nspCol + r ]!=0.0) {
        bIsZeroColumn = false;
        break;
      }
    }
    return bIsZeroColumn;
  }

  template <class Scalar, class LocalOrdinal>
  void QR_InterfaceEx<Scalar,LocalOrdinal>::Compute(LocalOrdinal const &myAggSize, ArrayRCP<Scalar> &localQR)
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

      // check whether there are zero columns in localQR
      nonZeroCols_.clear();
      for(LocalOrdinal j = 0; j<NSDim_; ++j) {
        if(!isZeroNspColumn(myAggSize, localQR, j))
            nonZeroCols_.push_back(j);
      }

      // set internal null space dimension used for QR decomposition
      internalNSDim_ = Teuchos::as<LocalOrdinal>(nonZeroCols_.size());

      // reset internal variables
      tau_ = ArrayRCP<Scalar>(internalNSDim_);
      work_ = ArrayRCP<Scalar>(internalNSDim_);
      workSize_ = internalNSDim_;
      info_ = 0;

      // build new input array
      internalLocalQR_ = Teuchos::ArrayRCP<Scalar>(myAggSize*internalNSDim_);
      for(size_t j = 0; j<internalNSDim_; ++j) {
        for(size_t r = 0; r<myAggSize; ++r) {
          internalLocalQR_[j*myAggSize+r] = localQR[nonZeroCols_[j]*myAggSize+r];
        }
      }

      // print new input matrix
      /*std::cout << " internal input matrix for GEQRF: " << std::endl;
      for (int k=0; k<myAggSize; ++k) {            // loop over rows
        for (size_t j=0; j<internalNSDim_; ++j) { // loop over columns
          std::cout << internalLocalQR_[j*myAggSize+k] << "\t";
        }
        std::cout << std::endl;
      }*/

      // perform lapack call
      lapack_.GEQRF( myAggSize, Teuchos::as<int>(internalNSDim_), internalLocalQR_.getRawPtr(), myAggSize,
                          tau_.getRawPtr(), work_.getRawPtr(), workSize_, &info_ );
      if (info_ != 0) {
        std::string msg = "QR_InterfaceEx: dgeqrf (LAPACK QR routine) returned error code " + Teuchos::toString(info_);
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

      // write back results to localQR
      // build new input array
      for(size_t j = 0; j<internalNSDim_; ++j) {
        for(size_t r = 0; r<internalNSDim_; ++r) {
          localQR[nonZeroCols_[j]*myAggSize+nonZeroCols_[r]] = internalLocalQR_[j*myAggSize+r];
        }
      }

      // print new input matrix
      /*std::cout << " localQR after internal GEQRF call: " << std::endl;
      for (int k=0; k<myAggSize; ++k) {            // loop over rows
        for (size_t j=0; j<NSDim_; ++j) { // loop over columns
          std::cout << localQR[j*myAggSize+k] << "\t";
        }
        std::cout << std::endl;
      }*/
    }
  } //Compute()

  template <class Scalar, class LocalOrdinal>
  void QR_InterfaceEx<Scalar,LocalOrdinal>::ExtractQ(LocalOrdinal const &myAggSize, ArrayRCP<Scalar> &localQR)
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
      // do LAPACK call on internalLocalQR_ array
      //lapack_.ORGQR(myAggSize, Teuchos::as<int>(workSize_), Teuchos::as<int>(workSize_), localQR.getRawPtr(),
      //             myAggSize, tau_.getRawPtr(), work_.getRawPtr(), workSize_, &info_ );
      //call nonmember function (perhaps specialized)
      LapackQR( lapack_, myAggSize, Teuchos::as<int>(internalNSDim_), internalLocalQR_, tau_, work_, workSize_, info_ );
      if (info_ != 0) {
        std::string msg = "QR_InterfaceEx: dorgqr (LAPACK auxiliary QR routine) returned error code " + Teuchos::toString(info_);
        throw(Exceptions::RuntimeError(msg));
      }
      // LAPACK may have determined a better length for the work array.  Returns it in work[0],
      // so we cast to avoid compiler warnings.
      // Scalar type might be complex, so take absolute value.
      if ( Teuchos::ScalarTraits<Scalar>::magnitude(work_[0]) > workSize_) {
        workSize_ = Teuchos::as<int>(Teuchos::ScalarTraits<Scalar>::magnitude(work_[0]));
        work_ = ArrayRCP<Scalar>(workSize_);
      }

      // write back results to localQR
      // build new input array
      for(size_t j = 0; j<internalNSDim_; ++j) {
        for(size_t r = 0; r<myAggSize; ++r) {
          localQR[nonZeroCols_[j]*myAggSize+r] = internalLocalQR_[j*myAggSize+r];
        }
      }

      // print new input matrix
      /*std::cout << " localQR after internal ORGQR call: " << std::endl;
      for (int k=0; k<myAggSize; ++k) {            // loop over rows
        for (size_t j=0; j<NSDim_; ++j) { // loop over columns
          std::cout << localQR[j*myAggSize+k] << "\t";
        }
        std::cout << std::endl;
      }*/
    }
  } //ExtractQ2()

} //namespace MueLu

#endif /* MUELU_QR_INTERFACEEX_DEF_HPP_ */
