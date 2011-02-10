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
#include <BelosOperatorTraits.hpp>
#include <Teuchos_BLAS.hpp>
#include <Teuchos_OrdinalTraits.hpp>
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

  /// \class OpAkx
  /// \brief Base for Akx implementations using only abstract operators
  /// \author Mark Hoemmen
  ///
  template<class Scalar, class MV, class OP>
  class OpAkx : public Akx<Scalar, MV> {
  public:
    /// \brief Constructor.
    ///
    /// \param A [in] The matrix of interest.
    /// \param M_left [in] If not null, the left preconditioner,
    ///   or split preconditioner if M_right is also not null.
    /// \param M_right [in] If not null, the right preconditioner,
    ///   or split preconditioner if M_left is also not null.
    OpAkx (const Teuchos::RCP<const OP>& A, 
	   const Teuchos::RCP<const OP>& M_left, 
	   const Teuchos::RCP<const OP>& M_right) :
      A_ (A), M_left_ (M_left), M_right_ (M_right)
    {}
	
    int maxCandidateBasisLength () const { 
      // We're not using a special matrix powers kernel implementation
      // here, just whatever operator(s) we were given, so we can
      // apply them as many times as we want.  Of course, this doesn't
      // account for memory or time constraints.
      return Teuchos::OrdinalTraits<int>::max(); 
    }

  protected:
    /// \brief Apply the "operator" M_left*A*M_right.
    ///
    /// The operator is suitable for (CA-)GMRES and for
    /// unpreconditioned iterative methods (where M_left and M_right
    /// are each the identity operator).  This method computes v_cur
    /// := M_left * (A * (M_right * v_prv)), where "*" indicates the
    /// "apply" operation which Belos::OperatorTraits::Apply()
    /// performs.  
    ///
    /// This method has been made available so that subclasses don't
    /// have to access A, M_left_, or M_right_ directly, and also to 
    /// avoid duplicated code.
    ///
    /// \param v_prv [in] Vector to which to apply the operator.
    /// \param v_cur [out] Result of applying the operator to z_cur.
    void 
    applyOp (const MV& v_prv, MV& v_cur) 
    {
      // We use v_scratch_ as scratch space, since
      // OperatorTraits::Apply() doesn't allow aliasing of the input
      // and output vectors when applying an operator.  We only
      // allocate v_scratch_ if necessary.
      if (! M_right_.is_null())
	{
	  if (v_scratch_.is_null())
	    v_scratch_ = MVT::Clone (v_cur);
	  OPT::Apply (*M_right_, v_prv, *v_scratch_);
	  OPT::Apply (*A_, *v_scratch_, v_cur);
	  if (! M_left_.is_null())
	    {
	      MVT::Assign (v_cur, *v_scratch_);
	      OPT::Apply (*M_left_, *v_scratch_, v_cur);
	    }
	}
      else if (! M_left_.is_null())
	{
	  if (v_scratch_.is_null())
	    v_scratch_ = MVT::Clone (v_cur);
	  OPT::Apply (*A_, v_prv, *v_scratch_);
	  OPT::Apply (*M_left_, *v_scratch_, v_cur);
	}
      else
	OPT::Apply (*A_, v_prv, v_cur);
    }

    /// \brief Apply the "flexible operator."
    ///
    /// The flexible operator is suitable for (CA-) Flexible GMRES.
    /// This method first applies the right preconditioner to v_prv,
    /// storing the result in z_cur.  If the right preconditioner is
    /// the identity operator, v_prv is just copied into z_cur.  Then,
    /// this method applies the matrix A to z_cur, storing the result
    /// in v_cur.
    ///
    /// This method has been made available so that subclasses don't
    /// have to access A, M_left_, or M_right_ directly, and also to 
    /// avoid duplicated code.
    ///
    /// \param v_prv [in] Vector to which to apply the right
    ///   preconditioner.  
    /// \param z_cur [out] Result of applying the right preconditioner
    ///   to v_cur.
    /// \param v_cur [out] Result of applying the matrix A to z_cur.
    void
    applyFlexibleOp (const MV& v_prv, MV& z_cur, MV& v_cur)
    {
      TEST_FOR_EXCEPTION(! M_left_.is_null(), std::logic_error,
			 "Flexible GMRES only works with a right "
			 "preconditioner, not a left or split "
			 "preconditioner.");
      if (! M_right_.is_null())
	{ 
	  OPT::Apply (*M_right_, v_prv, z_cur);
	  OPT::Apply (*A_, z_cur, v_cur);
	}
      else
	{ // In this case, the right preconditioner is the identity
	  // operator, so we only have to copy v_prv into z_cur.
	  MVT::Assign (v_prv, z_cur);
	  OPT::Apply (*A_, z_cur, v_cur);
	}
    }

  private:
    //! The matrix A of interest 
    Teuchos::RCP<const OP> A_;
    //! The left preconditioner (null means the identity operator)
    Teuchos::RCP<const OP> M_left_;
    //! The right preconditioner (null means the identity operator)
    Teuchos::RCP<const OP> M_right_;
    /// \brief Scratch space for applying preconditioners.
    ///
    /// This vector is allocated only when necessary, and the
    /// allocation is kept for later use.  
    ///
    /// FIXME (mfh 09 Feb 2011) Since Belos::MultiVecTraits doesn't
    /// have a way to compare different multivectors' Maps, caching
    /// the allocation here is subject to possible bugs if the
    /// operators are changed such that their inputs or outputs have
    /// different Maps.
    Teuchos::RCP<const MV> v_scratch_;
  };

  /// \class MonomialOpAkx
  /// \brief Monomial-basis Akx implementation using abstract operators.
  /// \author Mark Hoemmen
  ///
  template<class Scalar, class MV, class OP>
  class MonomialOpAkx : public OpAkx<Scalar, MV, OP> {
  public:
    MonomialOpAkx (const Teuchos::RCP<const OP>& A, 
		   const Teuchos::RCP<const OP>& M_left,
		   const Teuchos::RCP<const OP>& M_right) :
      OpAkx (A, M_left, M_right) {}

    void 
    computeBasis (const MV& q_last, MV& V_cur, const int s) 
    {
      using Teuchos::Range1D;
      using Teuchos::RCP;
      using Teuchos::rcp_const_cast;
      typedef MultiVecTraits<Scalar, MV> MVT;

      TEST_FOR_EXCEPTION(s < 0, std::invalid_argument, 
			 "Number of basis vectors to compute (s = " 
			 << s << ") is invalid; it must be positive.");
      if (s == 0)
	return; // Nothing to do
      TEST_FOR_EXCEPTION(MVT::GetNumberVecs(V_cur) < s, std::invalid_argument,
			 "You're asking to compute s = " << s << " basis "
			 "vectors, but only have room for " 
			 << MVT::GetNumberVecs(V_cur) << ".");
      RCP<const MV> v_prv = Teuchos::rcpFromRef (q_last);
      for (int k = 0; k < s; ++k, v_prv)
	{
	  RCP<MV> v_cur = MVT::CloneViewNonConst (V_cur, Range1D(k,k));
	  applyOp (v_prv, v_cur);
	  v_prv = rcp_const_cast<const MV> (v_cur);
	}
    }

    void 
    computeFlexibleBasis (const MV& q_last, 
			  MV& Z_cur, 
			  MV& V_cur, 
			  const int s)
    {
      using Teuchos::Range1D;
      using Teuchos::RCP;
      using Teuchos::rcp_const_cast;
      typedef MultiVecTraits<Scalar, MV> MVT;

      TEST_FOR_EXCEPTION(s < 0, std::invalid_argument, 
			 "Number of basis vectors to compute (s = " 
			 << s << ") is invalid; it must be positive.");
      if (s == 0)
	return; // Nothing to do
      TEST_FOR_EXCEPTION(MVT::GetNumberVecs(V_cur) < s, std::invalid_argument,
			 "You're asking to compute s = " << s << " basis "
			 "vectors, but only have room for " 
			 << MVT::GetNumberVecs(V_cur) << " in V_cur.");
      TEST_FOR_EXCEPTION(MVT::GetNumberVecs(Z_cur) < s, std::invalid_argument,
			 "You're asking to compute s = " << s << " basis "
			 "vectors, but only have room for " 
			 << MVT::GetNumberVecs(Z_cur) << " in Z_cur.");
      RCP<const MV> v_prv = Teuchos::rcpFromRef (q_last);
      for (int k = 0; k < s; ++k, v_prv)
	{
	  RCP<MV> z_cur = MVT::CloneViewNonConst (Z_cur, Range1D(k,k));
	  RCP<MV> v_cur = MVT::CloneViewNonConst (V_cur, Range1D(k,k));
	  applyFlexibleOp (v_prv, z_cur, v_cur);
	  v_prv = rcp_const_cast<const MV> (v_cur);
	}
    }

    Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,Scalar> > 
    changeOfBasisMatrix (const int s) 
    {
      using Teuchos::rcp;
      typedef Teuchos::SerialDenseMatrix<int,Scalar> mat_type;
      typedef Teuchos::ScalarTraits<Scalar> STS;
      
      if (B_.is_null())
	B_ = makeChangeOfBasisMatrix (s);
      else if (B_->numRows() != s+1 || B_->numCols() != s)
	B_ = makeChangeOfBasisMatrix (s);
      return B_;
    }

    void
    updateBasis (const Teuchos::ArrayView<const magnitude_type>& realParts,
		 const Teuchos::ArrayView<const magnitude_type>& imagParts,
		 const Teuchos::ArrayView<const int>& multiplicities,
		 const int numValues)
    { // The monomial basis doesn't do anything with approximate
      // eigenvalue information.  Quiet any compiler warnings about
      // unused values.
      (void) realParts;
      (void) imagParts;
      (void) multiplicities;
      (void) numValues;
    }

  private:
    //! Construct and return a new change-of-basis matrix.
    static Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,Scalar> >
    makeChangeOfBasisMatrix (const int s)
    {
      using Teuchos::RCP;
      typedef Teuchos::SerialDenseMatrix<int,Scalar> mat_type;
      typedef Teuchos::ScalarTraits<Scalar> STS;

      RCP<mat_type> B (new mat_type (s+1, s));
      for (int j = 0; j < s; ++j)
	(*B)(j+1, j) = STS::one();
      return B;
    }

    /// \brief Monomial change-of-basis matrix.
    ///
    /// The matrix is cached to avoid unnecessary allocations.
    Teuchos::RCP<const Teuchos::SerialDenseMatrix<int,Scalar> > B_;
  };

} // namespace Belos

#endif // __Belos_Akx_hpp
