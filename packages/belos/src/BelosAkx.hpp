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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
//
// ************************************************************************
//@HEADER

#ifndef __Belos_Akx_hpp
#define __Belos_Akx_hpp

#include <BelosMultiVecTraits.hpp>
#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_SerialDenseVector.hpp>

namespace Belos {
  /// \class Akx
  /// \brief Matrix powers kernel interface
  /// \author Mark Hoemmen
  ///
  template<class Scalar, class MV>
  class Akx {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    typedef MV multivector_type;

    /// \brief The maximum length of the candidate basis.
    ///
    /// Some implementations of the matrix powers kernel may set an
    /// upper bound on the candidate basis length.  This may be
    /// because initial set up of the matrix powers kernel requires
    /// expensive reorganization of the sparse matrix data structures.
    virtual int maxCandidateBasisLength () const = 0;

    /// \brief Recommended initial candidate basis length.
    ///
    /// Some implementations of the matrix powers kernel may recommend
    /// a particular candidate basis length, at least initially, for
    /// reasons of performance, numerical stability, or both.  Solvers
    /// should revise this dynamically for numerical stability, since
    /// the matrix powers kernel implementation does not promise to do
    /// the extensive analysis necessary to guarantee numerical
    /// stability at this candidate basis length.
    virtual int recommendedCandidateBasisLength () const = 0;

    /// \brief Compute flexible matrix powers kernel.
    ///
    /// For Flexible GMRES only.  Given the first vector q = q_last:
    /// Fill in the first s columns of V_cur with a candidate basis
    /// for the span of the vectors \f$A M^{-1} q, (A M^{-1})^2 q,
    /// \dots, (A M^{-1})^s q\f$, and fill in the first s columns of
    /// Z_cur with a candidate basis for the span of the vectors
    /// \f$M^{-1} q, M^{-1} (A M^{-1}) q, \dots, M^{-1} (A
    /// M^{-1})^{s-1} q\f$.  "Candidate basis" means that the vectors
    /// will be linearly independent in exact arithmetic, as long as
    /// the corresponding Krylov subspace has at least that dimension.
    /// Any remaining columns of V_cur or Z_cur will not be modified.
    ///
    /// We call V_cur the "right basis" because the preconditioner is
    /// applied on the right side of the operator, and Z_cur the "left
    /// basis" because the preconditioner is applied to the left side.
    /// If the preconditioner is the identity matrix (a "trivial
    /// preconditioner"), then the first s-1 columns of the left basis
    /// coincide with the last s-1 columns of the right basis.
    ///
    /// \param q_last [in] First vector of the matrix powers kernel.
    /// \param Z_cur [out] The first s columns will be filled in with
    ///   the right basis.
    /// \param V_cur [out] The first s columns will be filled in with
    ///   the left basis.
    /// \param s [in] Number of vectors to compute.  s >= 0.  If s ==
    ///   0, no vectors will be computed.
    virtual void 
    computeFlexibleBasis (const MV& q_last, 
			  MV& Z_cur, 
			  MV& V_cur, 
			  const int s) = 0;

    /// \brief Compute matrix powers kernel with first vector q_last.
    ///
    /// Fill in the first s columns of V_cur with the results of the
    /// matrix powers kernel (not including the first vector).  Any
    /// remaining columns of V_cur will not be modified.
    ///
    /// \param q_last [in] First vector of the matrix powers kernel.
    /// \param V_cur [out] The first s columns will be filled in with
    ///   the results of the matrix powers kernel.
    /// \param s [in] Number of vectors to compute.  s >= 0.  If s ==
    ///   0, no vectors will be computed.
    virtual void computeBasis (const MV& q_last, MV& V_cur, const int s) = 0;

    /// \brief Improve the basis by updating eigenvalue information.
    ///
    /// Some matrix powers kernel implementations may benefit (in
    /// terms of numerical stability) from (improvements to)
    /// approximations to the eigenvalues of the operator.  Many
    /// Krylov subspace methods (both for linear systems and for
    /// eigenvalue problems) naturally produce such approximations,
    /// which may improve as the iteration progresses.
    ///
    /// Use the given numValues eigenvalue approximations possibly to
    /// improve the numerical properties (e.g., condition number) of
    /// the matrix powers kernel basis.  This call assumes that
    /// eigenvalues may be complex.  Real parts are stored in the
    /// first array, and their corresponding imaginary parts in the
    /// second array.
    ///
    ///
    /// \param realParts [in] Real parts of the eigenvalue
    ///   approximations.  realParts[k] + i*imagParts[k] is the k-th
    ///   approximate eigenvalue, of multiplicity multiplicities[k].
    ///   Only the first numValues entries are read.
    ///
    /// \param imagParts [in] Imaginary parts of the eigenvalue
    ///   approximations.  realParts[k] + i*imagParts[k] is the k-th
    ///   approximate eigenvalue, of multiplicity multiplicities[k].
    ///   If no eigenvalues are complex, all of the entries of
    ///   imagParts must be zero.  Only the first numValues entries
    ///   are read.
    ///
    /// \param multiplicities [in] Multiplicity of the k-th eigenvalue
    ///   approximation.  Only the first numValues entries are read.
    ///
    /// \param numValues [in] The number of eigenvalue approximations,
    ///   not counting multiplicities.
    virtual void
    updateBasis (const Teuchos::ArrayView<const magnitude_type>& realParts,
		 const Teuchos::ArrayView<const magnitude_type>& imagParts,
		 const Teuchos::ArrayView<const int>& multiplicities,
		 const int numValues);

    /// \brief Return the s+1 by s change-of-basis matrix.
    /// 
    /// For details on the meaning and entries of this matrix, see
    /// M. Hoemmen, "Communication-avoiding Krylov subspace methods"
    /// (PhD thesis, UC Berkeley EECS), 2010.  If s is bigger than the
    /// max candidate basis length, this method will throw an
    /// exception.
    ///
    /// \note Implementations of the Newton basis should by default
    /// set all the shifts to zero, so the change-of-basis matrix is
    /// always valid even if the implementation has not yet received
    /// approximate eigenvalue information from the iterative method.
    /// The change-of-basis matrix should also be valid for any s up
    /// to the maximum candidate basis length, even if updateBasis() 
    /// has not yet been called with numValues >= that s value.
    ///
    /// \param s [in] The change-of-basis matrix returned by this
    ///   method will be s+1 by s.
    ///
    /// \return s+1 by s change-of-basis matrix.
    virtual Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,Scalar> > 
    changeOfBasisMatrix (const int s) = 0;
  };

} // namespace Belos

#endif // __Belos_Akx_hpp
