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
#ifndef MUELU_QR_INTERFACE_DECL_HPP
#define MUELU_QR_INTERFACE_DECL_HPP

#include <Teuchos_LAPACK.hpp>
#include <Teuchos_ArrayRCP.hpp>
#include "MueLu_ConfigDefs.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include "Teuchos_SerialQRDenseSolver.hpp"

namespace MueLu {

  /*!
    @class QR_Interface class.
    @brief Interface to the Teuchos wrappers for LAPACK's QR.

    Used in the computation of the tentative prolongator during smoothed aggregation.
    Allows for specializations for Scalar types such as std::complex<double> and PCE.
  */
  template <class Scalar, class LocalOrdinal>
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
      //! QR solver object.
      Teuchos::SerialQRDenseSolver<LocalOrdinal,Scalar> qrSolver_;
      //! Length of work vectors. Must be at least dimension of nullspace.
      LocalOrdinal     workSize_;
      //! Dimension of nullspace
      LocalOrdinal     NSDim_;
      //! (out) =0: success; =i, i<0: i-th argument has illegal value
      LocalOrdinal     info_;
      //! Internal LAPACK variables.
      ArrayRCP<Scalar> tau_;
      //! Temporary work space.
      ArrayRCP<Scalar> work_;
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

#define MUELU_QR_INTERFACE_SHORT

#endif //MUELU_QR_INTERFACE_DECL_HPP
