#ifndef MUELU_QR_INTERFACE_DECL_HPP
#define MUELU_QR_INTERFACE_DECL_HPP

#include <Teuchos_LAPACK.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include "MueLu_ConfigDefs.hpp"

#ifdef HAVE_MUELU_STOKHOS
#include "linear2d_diffusion_scalar_types.hpp"
#endif

namespace MueLu {

  /*!
    @class QR_Interface class.
    @brief Interface to the Teuchos wrappers for LAPACK's QR.

    Used in the computation of the tentative prolongator during smoothed aggregation.
    Allows for specializations for Scalar types such as std::complex<double> and PCE.
  */
  template <class Scalar, class Storage, class LocalOrdinal>
  class QR_Interface {

    public:
      //! @name Constructors/Destructors.
      //@{
      /*! @brief Constructor
        @param nullSpaceDimension number of fine level nullspace vectors.
      */
      QR_Interface(const size_t nullSpaceDimension);

      //! Destructor.
      virtual ~QR_Interface() {};
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
      //! Teuchos LAPACK wrapper.
      Teuchos::LAPACK<LocalOrdinal,Scalar> lapack_;
      //! Length of work vectors. Must be at least dimension of nullspace.
      LocalOrdinal     workSize_;
      //! (out) =0: success; =i, i<0: i-th argument has illegal value
      LocalOrdinal     info_;
      //! Internal LAPACK variables.
      ArrayRCP<Scalar> tau_;
      //! Temporary work space.
      ArrayRCP<Scalar> work_;
  }; //class QR_Interface

#if defined(HAVE_MUELU_STOKHOS) and defined(MUELU_SCALAR_IS_PCE_TYPE)
  /*!
    @brief Specialization for polynomial chaos expansion (PCE) scalar types.
  */
  template <class Scalar, class Storage, class LocalOrdinal, class GlobalOrdinal>
  class QR_Interface< Sacado::PCE::OrthogPoly<Scalar, Storage>, LocalOrdinal, GlobalOrdinal> {
    public:
      //! @name Constructors/Destructors.
      //@{
      /*! @brief Constructor
        @param nullSpaceDimension number of fine level nullspace vectors.
      */
      QR_Interface(const size_t nullSpaceDimension);

      //! Destructor.
      virtual ~QR_Interface() {};
      //@}

      //! Compute the QR factorization.
      void Compute(LocalOrdinal const &myAggSize, ArrayRCP<Sacado::PCE::OrthogPoly<Scalar, Storage> > &localQR);


      /*! @brief Calculate the Q factor.

        @param[in] myAggSize   number of points in current aggregate
        @param[in,out] localQR array that on input implicitly contains the factorization calculated by Compute.  On
                               output, this array explicitly contains Q.
      */
      void ExtractQ(LocalOrdinal const &myAggSize, ArrayRCP<Sacado::PCE::OrthogPoly<Scalar, Storage> > &localQR);

    private:
      //! Teuchos LAPACK wrapper.
      Teuchos::LAPACK<LocalOrdinal,Scalar> lapack_;
      //! Length of work vectors. Must be at least dimension of nullspace.
      LocalOrdinal     workSize_;
      //! (out) =0: success; =i, i<0: i-th argument has illegal value
      LocalOrdinal     info_;
      //! Internal LAPACK variables.
      ArrayRCP<Scalar> tau_;
      //! Temporary work space.
      ArrayRCP<Scalar> work_;
      //! POD (as opposed to PCE) local QR.
      ArrayRCP<Scalar> localQR_;
  };
#endif //ifdef HAVE_MUELU_STOKHOS

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

#define MUELU_QR_INTERFACE_SHORT

#endif //MUELU_QR_INTERFACE_DECL_HPP
