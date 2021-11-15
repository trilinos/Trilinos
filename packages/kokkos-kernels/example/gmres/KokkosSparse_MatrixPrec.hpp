/*
//@HEADER
// ************************************************************************
//
//                        Kokkos v. 3.0
//       Copyright (2020) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
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
// THIS SOFTWARE IS PROVIDED BY NTESS "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL NTESS OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Jennifer Loe (jloe@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
/// @file KokkosKernels_MatrixPrec.hpp

#ifndef KK_MATRIX_PREC_HPP
#define KK_MATRIX_PREC_HPP

#include<KokkosSparse_Preconditioner.hpp>
#include<Kokkos_Core.hpp>
#include<KokkosBlas.hpp>
#include<KokkosSparse_spmv.hpp>

namespace KokkosSparse{

namespace Experimental{

/// \class MatrixPrec
/// \brief  This is a simple class to use if one 
///         already has a matrix representation of their 
///         preconditioner M.  The class applies an
///         SpMV with M as the preconditioning step. 
/// \tparam ScalarType Type of the matrix's entries
/// \tparam Layout Kokkos layout of vectors X and Y to which 
///         the preconditioner is applied
/// \tparam EXSP Execution space for the preconditioner apply
/// \tparam Ordinal Type of the matrix's indices; 
///
/// Preconditioner provides the following methods
///   - initialize() Does nothing; Matrix initialized upon object construction. 
///   - isInitialized() returns true 
///   - compute() Does nothing; Matrix initialized upon object construction.
///   - isComputed() returns true 
///
template< class ScalarType, class Layout, class EXSP, class OrdinalType = int > 
class MatrixPrec : virtual public KokkosSparse::Experimental::Preconditioner<ScalarType, Layout, EXSP, OrdinalType>
{ 
private:
    using crsMat_t = KokkosSparse::CrsMatrix<ScalarType, OrdinalType, EXSP>;
    crsMat_t A;

    bool isInitialized_ = true;
    bool isComputed_ = true;

public:
  //! Constructor:
    MatrixPrec<ScalarType, Layout, EXSP, OrdinalType> (const KokkosSparse::CrsMatrix<ScalarType, OrdinalType, EXSP> &mat) 
     : A(mat) {}

  //! Destructor.
  virtual ~MatrixPrec(){}

  ///// \brief Apply the preconditioner to X, putting the result in Y.
  /////
  ///// \tparam XViewType Input vector, as a 1-D Kokkos::View
  ///// \tparam YViewType Output vector, as a nonconst 1-D Kokkos::View
  /////
  ///// \param transM [in] "N" for non-transpose, "T" for transpose, "C"
  /////   for conjugate transpose.  All characters after the first are
  /////   ignored.  This works just like the BLAS routines.
  ///// \param alpha [in] Input coefficient of M*x
  ///// \param beta [in] Input coefficient of Y
  /////
  ///// If the result of applying this preconditioner to a vector X is
  ///// \f$M \cdot X\f$, then this method computes \f$Y = \beta Y + \alpha M \cdot X\f$.
  ///// The typical case is \f$\beta = 0\f$ and \f$\alpha = 1\f$.
  //
  void apply (const Kokkos::View<ScalarType*, Layout, EXSP> &X, 
      Kokkos::View<ScalarType*, Layout, EXSP> &Y, 
      const char transM[] = "N",
      ScalarType alpha = Kokkos::Details::ArithTraits<ScalarType>::one(),
      ScalarType beta = Kokkos::Details::ArithTraits<ScalarType>::zero()) const 
  {
    KokkosSparse::spmv(transM, alpha, A, X, beta, Y);
  };
  //@}

  //! Set this preconditioner's parameters.
  void setParameters () {}

  void initialize() { }

  //! True if the preconditioner has been successfully initialized, else false.
  bool isInitialized() const { return isInitialized_;}

  void compute(){ }

  //! True if the preconditioner has been successfully computed, else false.
  bool isComputed() const {return isComputed_;}

  //! True if the preconditioner implements a transpose operator apply. 
  bool hasTransposeApply() const { return true; } 

};
} //End Experimental
} //End namespace KokkosSparse

#endif
