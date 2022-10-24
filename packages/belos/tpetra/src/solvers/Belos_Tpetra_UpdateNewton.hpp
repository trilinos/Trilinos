//@HEADER
// ************************************************************************
//
//                 Belos: Block Linear Solvers Package
//                  Copyright 2004 Sandia Corporation
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
// ************************************************************************
//@HEADER

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
