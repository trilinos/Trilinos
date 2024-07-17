// @HEADER
// *****************************************************************************
//                 Belos: Block Linear Solvers Package
//
// Copyright 2004-2016 NTESS and the Belos contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef BELOS_TPETRA_UPDATENEWTON_HPP
#define BELOS_TPETRA_UPDATENEWTON_HPP

#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_SerialDenseMatrix.hpp"
#include <complex>

namespace BelosTpetra {
namespace Impl {

template<class SC,
	 class MV,
	 const bool isComplex = Teuchos::ScalarTraits<SC>::isComplex>
struct UpdateNewton {};

template<class SC, class MV>
struct UpdateNewton<SC, MV, true> {
  using STS = Teuchos::ScalarTraits<SC>;
  using real_type = typename STS::magnitudeType;
  using complex_type = std::complex<real_type>;
  using dense_matrix_type = Teuchos::SerialDenseMatrix<int, SC>;

  //! Update the vector to convert a Monomial basis to a Newton basis.
  static void
  updateNewtonV (const int iter,
		 MV& V,
		 const complex_type& theta)
  {
    const SC one = STS::one ();
    auto Z = V.getVectorNonConst (iter);
    auto W = V.getVectorNonConst (iter+1);
      
    W->update (-theta, *Z, one);
  }

  /// \brief Update the change-of-basis matrix to convert a Monomial
  ///   basis to a Newton basis.
  static void
  updateNewtonH (const int iter,
		 dense_matrix_type& H,
		 const complex_type& theta)
  {
    H(iter, iter) += theta;
  }

  static SC
  updateNewtonH (const int i,
		 const int j,
		 dense_matrix_type& R,
		 const complex_type& theta)
  {
    return theta * R(i, j);
  }
};

template<class SC, class MV>
struct UpdateNewton<SC, MV, false> {
  using STS = Teuchos::ScalarTraits<SC>;
  using real_type = typename STS::magnitudeType;
  using STM = Teuchos::ScalarTraits<real_type>;    
  using complex_type = std::complex<real_type>;
  using dense_matrix_type = Teuchos::SerialDenseMatrix<int, SC>;
    
  static void
  updateNewtonV (const int iter,
		 MV& V,
		 const complex_type& theta)
  {
    const SC one = STS::one ();
    auto Z = V.getVectorNonConst (iter);
    auto W = V.getVectorNonConst (iter+1);

    W->update (-theta.real(), *Z, one);
    if (theta.imag() != STM::zero()) {
      if (iter > 0) {
	auto Zp = V.getVectorNonConst (iter-1);
	W->update (theta.imag()*theta.imag(), *Zp, one);
      }
    }
  }

  static void
  updateNewtonH (const int iter,
		 dense_matrix_type& H,
		 const complex_type& theta)
  {
    H(iter, iter) += theta.real();
    if (theta.imag() != 0.0) {
      if (iter > 0) {
	H(iter, iter-1) -= theta.imag()*theta.imag();
      }
    }
  }

  static SC
  updateNewtonH (const int i,
		 const int j,
		 dense_matrix_type& R,
		 const complex_type& theta)
  {
    real_type ret = theta.real() * R(i, j);
    if (theta.imag() != 0.0 && j > 0) {
      ret -= (theta.imag()*theta.imag()) * R(i, j-1);
    }
    return ret;
  }
};

} // namespace Impl
} // namespace BelosTpetra
  
#endif // BELOS_TPETRA_UPDATENEWTON_HPP
