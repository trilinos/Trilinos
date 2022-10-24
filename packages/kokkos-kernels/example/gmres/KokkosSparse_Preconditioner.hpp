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
/// @file KokkosKernels_Preconditioner.hpp
//
#ifndef KK_PREC_HPP
#define KK_PREC_HPP

#include<Kokkos_Core.hpp>
#include<KokkosKernels_Controls.hpp>
#include<Kokkos_ArithTraits.hpp>

namespace KokkosSparse{

namespace Experimental{

/// \class Preconditioner
/// \brief Interface for KokkosKernels preconditioners
/// \tparam ScalarType Type of the matrix's entries
/// \tparam Layout Kokkos layout of vectors X and Y to which 
///         the preconditioner is applied
/// \tparam EXSP Execution space for the preconditioner apply
/// \tparam Ordinal Type of the matrix's indices; 
///
/// Preconditioner provides the following methods
///   - initialize() performs all operations based on the graph of the
///     matrix (without considering the numerical values)
///   - isInitialized() returns true if the preconditioner has been
///     successfully initialized
///   - compute() computes everything required to apply the
///     preconditioner, using the matrix's values (and assuming that the
///     graph structure of the matrix has not changed)
///   - isComputed() returns true if the preconditioner has been
///     successfully computed, false otherwise.
///
/// Implementations of compute() must internally call initialize() if
/// isInitialized() returns false. The preconditioner is applied by
/// apply(). 
/// Every time that initialize() is called, the object destroys all the previously
/// allocated information, and reinitializes the preconditioner. Every
/// time compute() is called, the object recomputes the actual values of
/// the preconditioner.
template< class ScalarType, class Layout, class EXSP, class OrdinalType = int > 
class Preconditioner{ 
  
public:
  //! Constructor:
  Preconditioner(){}

  //! Destructor.
  virtual ~Preconditioner(){}

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
  virtual void
  apply (const Kokkos::View<ScalarType*, Layout, EXSP> &X, 
         Kokkos::View<ScalarType*, Layout, EXSP> &Y, 
         const char transM[] = "N",
         ScalarType alpha = Kokkos::Details::ArithTraits<ScalarType>::one(),
         ScalarType beta = Kokkos::Details::ArithTraits<ScalarType>::zero()) const = 0;
  //@}

  //! Set this preconditioner's parameters.
  virtual void setParameters () = 0;

  /// @brief Set up the graph structure of this preconditioner.
  ///
  /// If the graph structure of the constructor's input matrix has
  /// changed, or if you have not yet called initialize(), you must
  /// call initialize() before you may call compute() or apply().
  ///
  /// Thus, initialize() corresponds to the "symbolic factorization"
  /// step of a sparse factorization, whether or not the specific
  /// preconditioner actually does a sparse factorization.
  virtual void initialize() = 0;

  //! True if the preconditioner has been successfully initialized, else false.
  virtual bool isInitialized() const = 0;

  /// @brief Set up the numerical values in this preconditioner.
  ///
  /// If the values of the constructor's input matrix have changed, or
  /// if you have not yet called compute(), you must call compute()
  /// before you may call apply().
  ///
  /// Thus, compute() corresponds to the "numeric factorization"
  /// step of a sparse factorization, whether or not the specific
  /// preconditioner actually does a sparse factorization.
  virtual void compute() = 0;

  //! True if the preconditioner has been successfully computed, else false.
  virtual bool isComputed() const = 0;

  //! True if the preconditioner implements a transpose operator apply. 
  virtual bool hasTransposeApply() const { return false; } 

};

} // End Experimental
} //End namespace KokkosSparse

#endif
