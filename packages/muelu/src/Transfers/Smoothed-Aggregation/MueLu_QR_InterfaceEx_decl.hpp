/*
 * MueLu_QR_InterfaceEx_decl.hpp
 *
 *  Created on: Jan 31, 2013
 *      Author: tobias
 */

#ifndef MUELU_QR_INTERFACEEX_DECL_HPP_
#define MUELU_QR_INTERFACEEX_DECL_HPP_

#include <Teuchos_LAPACK.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include "MueLu_ConfigDefs.hpp"

namespace MueLu {

  /*!
    @class QR_InterfaceEx class.
    @brief Interface to the Teuchos wrappers for LAPACK's QR.

    Used in the computation of the tentative prolongator during smoothed aggregation.
    Allows for specializations for Scalar types such as std::complex<double> and PCE.
  */
  template <class Scalar, class LocalOrdinal>
  class QR_InterfaceEx {

    public:
      //! @name Constructors/Destructors.
      //@{
      /*! @brief Constructor
        @param nullSpaceDimension number of fine level nullspace vectors.
      */
      QR_InterfaceEx(const size_t nullSpaceDimension);

      //! Destructor.
      virtual ~QR_InterfaceEx() {};
      //@}

      //! @name Computational methods.
      //@{
      //! Compute the QR factorization.
      void Compute(LocalOrdinal const &myAggSize, ArrayRCP<Scalar> &localQR);

      /*! @brief Calculate the Q factor.

        @param[in] myAggSize   number of points in current aggregate
        @param[in,out] localQR array that on input implicitly contains the factorization calculated by Compute.  On
                               output, this array explicitly contains Q.
      */
      void ExtractQ(LocalOrdinal const &myAggSize, ArrayRCP<Scalar> &localQR);
      //@}

    private:
      bool isZeroNspColumn(LocalOrdinal const &myAggSize, ArrayRCP<Scalar> &localQR, LocalOrdinal nspCol);

      //! Teuchos LAPACK wrapper.
      Teuchos::LAPACK<LocalOrdinal,Scalar> lapack_;
      //! Length of work vectors. Must be at least dimension of nullspace.
      LocalOrdinal     workSize_;
      //! Dimension of nullspace
      LocalOrdinal     NSDim_; // null space dimension of original input

      LocalOrdinal     internalNSDim_; // internal null space dimension (after removing zero columns of input)
      //! (out) =0: success; =i, i<0: i-th argument has illegal value
      LocalOrdinal     info_;
      //! Internal LAPACK variables.
      ArrayRCP<Scalar> tau_;
      //! Temporary work space.
      ArrayRCP<Scalar> work_;

      //! internal local QR (without zero columns)
      ArrayRCP<Scalar> internalLocalQR_;

      std::vector<LocalOrdinal> nonZeroCols_;
  }; //class QR_Interface

  //! @brief Non-member templated function to handle extracting Q from QR factorization.
  template <class Scalar, class LocalOrdinal>
  void LapackQR(Teuchos::LAPACK<LocalOrdinal,Scalar> &lapack, LocalOrdinal myAggSize,
                int intFineNSDim, ArrayRCP<Scalar> &localQR, ArrayRCP<Scalar> &tau,
                ArrayRCP<Scalar> &work, LocalOrdinal &workSize, LocalOrdinal &info);

  //! @brief Non-member specialized function to handle extracting Q from QR factorization for Scalar==std::complex<double>
  template <class LocalOrdinal>
  void LapackQR(Teuchos::LAPACK<LocalOrdinal, std::complex<double> > &lapack,
                LocalOrdinal myAggSize, int intFineNSDim, ArrayRCP<std::complex<double> > &localQR,
                ArrayRCP<std::complex<double> > &tau, ArrayRCP<std::complex<double> > &work,
                LocalOrdinal &workSize, LocalOrdinal &info);

} //namespace MueLu

#define MUELU_QR_INTERFACEEX_SHORT


#endif /* MUELU_QR_INTERFACEEX_DECL_HPP_ */
