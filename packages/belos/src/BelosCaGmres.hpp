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

#ifndef __Belos_CaGmres_hpp
#define __Belos_CaGmres_hpp

/// \file BelosCaGmres.hpp
/// \brief Implementation of Communication-Avoiding GMRES (CA-GMRES)
/// \author Mark Hoemmen

#include "BelosConfigDefs.hpp"
#include "BelosTypes.hpp"
#include "BelosGmresBase.hpp"

namespace Belos {

  //! @name CaGmres Exceptions
  //@{
  
  /// Raised when CaGmres::updateUpperHessenberg() fails.  There's no
  /// good way to recover from such a failure: it means that even
  /// though the previous upper Hessenberg matrix _and_ the new basis
  /// are both full rank, somehow the update resulted in a singular
  /// matrix.
  class CaGmresUpdateFailure : public BelosError {
  public:
    CaGmresUpdateFailure(const std::string& what_arg) : 
      BelosError(what_arg) {}
  };

  
  /// CaGmresOrthoFailure is thrown when the CaGmres object is unable
  /// to compute at least two independent vectors in the advance()
  /// routine.
  class CaGmresOrthoFailure : public BelosError {
  public:
    CaGmresOrthoFailure (const int theCurOuterIter,
			 const int theNumBasisVectors,
			 const int theCandidateBasisLength,
			 const int theRank) :
      curOuterIter_ (theCurOuterIter),
      numBasisVectors_ (theNumBasisVectors)
      candidateBasisLength_ (theCandidateBasisLength),
      rank_ (theRank)
    {}

    //! Current outer iteration
    int curOuterIter() const { return curOuterIter_; }

    //! Total number of basis vectors successfully constructed
    int numBasisVectors() const { return numBasisVectors_; }

    //! The offending candidate basis' length
    int candidateBasisLength() const { return candidateBasisLength_; }

    //! The offending candidate basis' rank
    int rank() const { return rank_; }

    //! Return a human-readable, informative error message
    const char* const 
    what()
    {
      // Lazy initialization via an std::string avoids memory leaks,
      // if we decide to catch the exception and handle it (rather
      // than pass it on up the call chain).  It's likely that we may
      // want to handle the exception directly, perhaps by decreasing
      // the basis length and retrying the outer iteration, so
      // avoiding memory leaks matters in this case.
      if (what_ == "")
	{
	  // This may throw an exception if there isn't enough memory
	  // with which to allocate a string, but then the program is
	  // likely in deep trouble anyway.
	  std::ostringstream os;
	  os << "CA-GMRES: After constructing " << numBasisVectors_ 
	     << " successfully, outer iteration" << (curOuterIter_+1) 
	     << " failed: the candidate basis of length " 
	     << candidateBasisLength_ << " only has rank " 
	     << rank_ << ".";
	  what_ = os.str();
	}
      return what_.c_str();
    }

  private:
    int curOuterIter_, numBasisVectors_, candidateBasisLength_, rank_;
    std::string what_;
  };
  
  //@}

  /// \class CaGmres
  /// \brief Communication-Avoiding GMRES implementation
  /// \ingroup belos_solver_framework
  /// \author Mark Hoemmen
  ///
  /// Communication-Avoiding GMRES (CA-GMRES) (Hoemmen 2010) is an
  /// iterative method for solving nonsymmetric linear systems.  It is
  /// equivalent to GMRES (Saad and Schultz 1986) in exact arithmetic.
  ///
  /// References:
  /// - M. Hoemmen, "Communication-avoiding Krylov subspace methods,"
  ///   PhD dissertation, University of California Berkeley EECS 
  ///   Department, 2010.
  /// - Saad and Schultz, "GMRES: A generalized minimal residual
  ///   algorithm for solving nonsymmetric linear systems", SISSC,
  ///   vol. 7, pp. 856-869, 1986.
  /// - Saad, "A flexible inner-outer preconditioned GMRES algorithm",
  ///   SISC, vol. 14, pp. 461-469, 1993.
  ///
  template<class Scalar, class MV, class OP>
  class CaGmres : public GmresBase<Scalar, MV, OP> {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;

    /// Constructor
    ///
    /// \param problem [in/out] Linear problem.  On input, we use the
    ///   starting guess (x0) and the initial residual (r0) to
    ///   initialize the iteration.  The iteration may call
    ///   updateSolution().  On output, if the solution has been
    ///   updated, the vector returned by getLHS() will be modified.
    /// \param ortho [in] Orthogonalization manager
    /// \param maxIterCount [in] Maximum number of iterations before
    ///   restart.  The number of vectors' worth of storage this
    ///   constructor allocates is proportional to this, so choose
    ///   carefully.
    /// \param flexible [in] Whether or not to run the Flexible
    ///   variant of GMRES (FGMRES).  This requires twice as much
    ///   vector storage, but lets the preconditioner change in every
    ///   iteration.  This only works with right preconditioning, not
    ///   left or split preconditioning.
    /// \param akx [in/out] Callers have the option to provide a
    ///   previously initialized matrix powers kernel.  We offer this
    ///   option both because the LinearProblem doesn't have a slot
    ///   for a matrix powers kernel, and because constructing an
    ///   optimized matrix powers kernel is possibly expensive but
    ///   need only be done once if the matrix and preconditioner
    ///   don't change.  It's allowed for akx to be null, in which
    ///   case we construct a matrix powers kernel implementation
    ///   ourselves.  The object is nonconst because the iteration
    ///   itself is used to update eigenvalue approximations, which
    ///   the matrix powers kernel uses to improve numerical 
    ///   stability.
    /// \param params [in] Options.  May be null, in which case
    ///   default options are used.
    ///
    /// \note The Flexible GMRES (FGMRES) option forces the matrix
    ///   powers kernel to use a conventional implementation, rather
    ///   than an optimized matrix powers kernel implementation.
    ///   Besides this, CA-GMRES works as usual.  However, if the
    ///   varying preconditioner tends to produce a numerically
    ///   rank-deficient candidate basis, this will restrict the
    ///   candidate basis length.
    CaGmres (const Teuchos::RCP<LinearProblem<Scalar,MV,OP> >& lp,
	     const Teuchos::RCP<const OrthoManager<Scalar, MV> >& ortho,
	     const int maxNumIters,
	     const bool flexible,
	     const Teuchos::RCP<Akx<Scalar, MV> >& akx,
	     const Teuchos::RCP<const Teuchos::ParameterList>& params)
      : GmresBase (lp, ortho, maxNumIters, flexible),
	akx_ (akx),
	numOuterIters_ (0)
    {
      initialize (akx, maxNumIters, flexible, params);
    }

    void 
    initAkx (const Teuchos::RCP<Akx<Scalar, MV> >& akx, 
	     const Teuchos::RCP<Teuchos::ParameterList>& params,
	     const bool flexible)
    {
      if (! is_null(akx))
	akx_ = akx;
      else
	{
	  typedef AkxFactory<Scalar, MV, OP> factory_type;
	  Teuchos::RCP<Teuchos::ParameterList> akxParams;
	  if (! is_null(params))
	    {
	      const char* candidateNames[] = {"Akx", "Matrix Powers Kernel"};
	      const int numCandidateNames = 2;
	      for (k = 0; k < numCandidateNames && is_null(akxParams); ++k)
		{
		  try {
		    *akxParams = params->sublist (candidateNames[k]);
		  } catch (???) {
		    akxParams = null;
		  }
		}
	    }
	  akx_ = factory_type::create (lp, flexible, akxParams);
	}
    }

    void 
    initialize (const Teuchos::RCP<Akx<Scalar, MV> >& akx, 
		const Teuchos::RCP<Teuchos::ParameterList>& params,
		const int maxNumIters,
		const bool flexible)
    {
      using Teuchos::is_null;
      using Teuchos::null;

      initAkx (akx, params, flexible);
      // FIXME (mfh 31 Dec 2010) Should we get the candidate basis
      // length from the Akx implementation alone?
      candidateBasisLength_ = akx_->candidateBasisLength();
      // FIXME (mfh 31 Dec 2010) Is this right?
      maxNumOuterIters_ = maxNumIters / candidateBasisLength_ + 1;
    }

    virtual bool canExtendBasis() const { 
      return getNumIters() < maxNumIters(); // FIXME (mfh 17 Jan 2011) ???
    }

    virtual void
    extendBasis (Teuchos::RCP<MV>& V_cur, 
		 Teuchos::RCP<MV>& Z_cur)
    {
      // This does not count the initial basis vector.
      const int m = getNumIters(); 
      TEST_FOR_EXCEPTION(m >= maxNumIters(), GmresCantExtendBasis,
			 "Maximum number of iterations " << getNumIters() 
			 << " reached.");
      // Reduce candidate basis length as necessary, if we are running
      // out of iterations to perform.  However, if this is the case,
      // don't change the saved candidate basis length; we want to
      // keep that in case we restart.  
      //
      // Note that the candidate basis length is not a property of the
      // matrix powers kernel implementation, which in turn is
      // required to support different basis lengths, up to a maximum
      // length with which it was initially configured.
      const int s = (m + candidateBasisLength_ >= maxNumIters()) ? 
	maxNumIters() - m : candidateBasisLength_;

      // The most recently orthogonalized and accepted basis vector in
      // the Q basis.
      RCP<const MV> q_last = MVT::CloneView (V_, Range1D(m, m));
      // Space for new candidate basis vectors in the Q basis.
      V_cur = MVT::CloneViewNonConst (V_, Range1D(m+1, m+s));
      if (flexible_)
	{
	  // Space for new candidate basis vectors in the Z basis.
	  Z_cur = MVT::CloneViewNonConst (Z_, Range1D(m, m+s-1));
	  akx_->computeFlexibleBasis (q_last, Z_cur, V_cur, s);
	}
      else
	{
	  if (! Z_cur.is_null())
	    Z_cur = Teuchos::null;
	  akx_->computeBasis (q_last, V_cur, s);
	}
      // In case we need to retry with a shorter candidate basis
      // length, remember what basis length we used this time.
      lastCandidateBasisLength_ = s;
    }

    virtual void
    acceptCandidateBasis ()
    {
      const int m = getNumIters(); 
      const int s = lastCandidateBasisRank_;

      curNumIters_ += s;
      numInnerIters_.push_back (s);
      numOuterIters_++;

      // FIXME (mfh 31 Dec 2010, 17 Jan 2011) Figure out how best to
      // express updating the approximate eigenvalue information in
      // the matrix powers kernel.  What should we pass in?
      akx_->updateBasis (m, s, H_);
    }

    virtual void
    rejectCandidateBasis ()
    {
      // Whether we should invoke advance() recursively.  The
      // alternative is to throw an exception signalling the inability
      // of CA-GMRES to make progress, given the linear system to
      // solve and the matrix powers kernel implementation.
      const bool shouldRecurse = 
	lastCandidateBasisRank_ >= 2 && candidateBasisLength_ > 2;

      // It's not worth running CA-GMRES with a candidate basis length
      // of 1, so give up if that would happen on recursion.
      if (! shouldRecurse)
	{
	  const std::ostringstream os;
	  os << "You should either use standard GMRES for solving this "
	    "problem, or try setting up the matrix powers kernel differently.";
	  // Rank of the current candidate basis is < 2.
	  TEST_FOR_EXCEPTION(lastCandidateBasisRank_ < 2, 
			     GmresRejectsCandidateBasis,
			     "CA-GMRES cannot produce a candidate basis of rank at "
			     "least 2.  " << os.str());
	  // Candidate basis length (not the "last candidate basis
	  // length," which may be different) is <= 2.
	  TEST_FOR_EXCEPTION(candidateBasisLength_ <= 2, 
			     GmresRejectsCandidateBasis,
			     "CA-GMRES: the current candidate basis length " 
			     << candidateBasisLength_ << " is <= 2.  This means"
			     " we cannot attempt recovery by recursing with a "
			     "shorter candidate basis length." << os.str());
	  TEST_FOR_EXCEPTION(true, std::logic_error, "Should never get here!");
	}

      // Decrease the candidate basis length and try again.  It could
      // have been that lastCandidateBasisLength_ was quite a bit less
      // than candidateBasisLength_, if we were at the end of a
      // restart cycle and candidateBasisLength_ does not divide the
      // restart length evenly.  However, since the last candidate
      // basis was still not full rank, we should still use
      // lastCandidateBasisLength_ as a guide.
      //
      // The max is to ensure the second term of the min is at least 2,
      // for example if candidateBasisLength == 3.
      const int nextCandidateBasisLength = 
	std::min (lastCandidateBasisRank_, 
		  std::max (2, candidateBasisLength_ / 2));
      // Throw an exception if a programming bug would otherwise have
      // caused an infinite loop.
      TEST_FOR_EXCEPTION(nextCandidateBasisLength < 2, std::logic_error, 
			 "CA-GMRES should never get here!");
      candidateBasisLength_ = nextCandidateBasisLength;

      // Invoke advance() recursively.  The cost of recursion should
      // be minimal, since typical candidate basis lengths are short,
      // and advance() (by itself, not counting subclass
      // implementations of hooks, etc.) doesn't allocate any dense
      // vector or sparse matrix data.
      advance();
    }

    virtual void 
    updateUpperHessenbergMatrix (const Teuchos::RCP<mat_type>& C_Q,
				 const Teuchos::RCP<mat_type>& B_Q,
				 const Teuchos::RCP<mat_type>& C_Z,
				 const Teuchos::RCP<mat_type>& B_Z)
    {
      using Teuchos::Copy, Teuchos::View, Teuchos::NO_TRANS, Teuchos::RCP;
      const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();

      // Sanity checks
      TEST_FOR_EXCEPTION(C_Q.is_null(), std::logic_error, "C_Q is null");
      TEST_FOR_EXCEPTION(B_Q.is_null(), std::logic_error, "B_Q is null");
      TEST_FOR_EXCEPTION(flexible_ && C_Z.is_null(), std::logic_error, 
			 "C_Z is null; not allowed when Flexible GMRES "
			 "is being used");
      TEST_FOR_EXCEPTION(flexible_ && B_Z.is_null(), std::logic_error, 
			 "B_Z is null; not allowed when Flexible GMRES "
			 "is being used");
      // Even on the first outer iteration, C_Q should be non-null; it
      // should contain exactly one row.
      const int m = C_Q->numRows();
      const int s = C_Q->numCols();
      if (flexible_)
	{
	  const int m_z = C_Z->numRows();
	  const int s_z = C_Q->numCols();
	  TEST_FOR_EXCEPTION(m_z != m || s_z != s-1, std::logic_error, 
			     "C_Z has the wrong dimensions: C_Q is " << m 
			     << " x " << s << " and so C_Z should be " << m 
			     << " x " << (s-1) << ", but is " << m_z << " x " 
			     << s_z << " instead.");
	}
      // Helper matrices make the algebra easier, though they take up
      // a modest amount of extra space.  Cleverness could make them
      // go away, but we are more interested in clarity at this point.
      mat_type F_Q (m-1, s+1);
      {
	if (m > 1)
	  {
	    mat_type F_Q_left (View, F_Q, m-1, 1, 0, 0);
	    F_Q_left.putScalar (zero);
	    mat_type F_Q_right (View, F_Q, m-1, s, 0, 1);
	    F_Q_right = mat_type (View, *C_Q, m-1, s);
	  }
      }
      mat_type F_Z (m-1, s);
      {
	if (m > 1)
	  {
	    mat_type F_Z_left (View, F_Z, m-1, 1, 0, 0);
	    F_Z_left.putScalar (zero);
	    mat_type F_Z_right (View, F_Z, m-1, s-1, 0, 1);
	    if (flexible_)
	      F_Z_right = mat_type (View, *C_Z, m-1, s-1);
	    else
	      F_Z_right = mat_type (View, *C_Q, m-1, s-1);
	  }
      }
      mat_type R_Q (s+1, s+1);
      {
	R_Q(0,0) = one;
	mat_type R_Q_top (View, R_Q, 1, s, 0, 1);
	R_Q_top = mat_type (View, *C_Q, 1, s, m-1, 1);
	mat_type R_Q_bot (View, R_Q, s, s, 1, 1);
	R_Q_bot = *B_Q;
      }
      mat_type R_Z (s, s);
      {
	R_Z(0,0) = one;
	mat_type R_Z_top (View, R_Z, 1, s-1, 0, 1);
	if (flexible_)
	  R_Z_top = mat_type (View, *C_Z, 1, s-1, m-1, 1);
	else
	  R_Z_top = mat_type (View, *C_Q, 1, s-1, m-1, 1);
	mat_type R_Z_bot (View, R_Z, s-1, s-1, 1, 1);
	R_Z_bot = *B_Z;
      }
      Teuchos::BLAS<int,Scalar> blas;
      
      // E := B * inv(R_Z), an s+1 by s matrix, is a handy
      // intermediate quantity.  In the previous sentence, B is the
      // s+1 by s change-of-basis matrix (yes, naming it B is
      // confusing).
      mat_type E (s+1, s);
      {
	using Teuchos::RIGHT_SIDE, Teuchos::UPPER_TRI, Teuchos::NON_UNIT_DIAG;
	// The change-of-basis matrix is s+1 by s.  We confusingly also
	// name it B.  It is the same whether or not we are doing
	// Flexible GMRES.
	RCP<const mat_type> B = akx_->changeOfBasisMatrix (s);
	// BLAS' TRSM overwrites its second matrix argument, so copy
	// the right-hand side into E before solving.
	E = *B; 
	// E := B * inv(R_Z)
	blas.TRSM (RIGHT_SIDE, UPPER_TRI, NON_TRANS, NON_UNIT_DIAG, s+1, s,
		   one, R_Z.values(), R_Z.stride(), E.values(), E.stride());
      }
      //
      // Compute H_kk, the lower right s+1 by s part of the upper
      // Hessenberg matrix.  Note that m-1 >= 0, even on the first
      // outer iteration (where m==1, since we still have to project
      // against the initial basis vector).
      //
      mat_type H_kk (View, *H_, s+1, s, m-1, m-1);
      // H_kk := R_Q * E (result is s+1 by s)
      blas.GEMM (NO_TRANS, NO_TRANS, s+1, s, s+1,
		 one, R_Q.values(), R.stride(), 
		 E.values(), E.stride(),
		 zero, H_kk.values(), H_kk.stride());
      // If this is the first outer iteration, then we're done.
      if (getNumIters() == 0)
	return;
      // H_kk(1, 2:s) -= beta_km1 F_Z(m-1, 2:s)
      const Scalar beta_km1 = (*H_)(m-1, m-2);
      {
	mat_type H_kk_top (View, H_kk, 1, s-1, 0, 1);
	const mat_type F_Z_bot (View, F_Z, 1, s-1, m-2, 1);
	H_kk_top = F_Z_bot;
	H_kk_top *= -beta_km1; // elementwise multiplication
      }
      //
      // Compute H_km1_k
      //
      mat_type H_km1_k = mat_type (View, *H_, m-1, s, 0, m-1);
      // H_km1_k := -H_km1 * F_Z (result is m-1 by s)
      const mat_type H_km1 = mat_type (View, *H_, m-1, m-1);
      blas.GEMM (NO_TRANS, NO_TRANS, m-1, s, m-1, 
		 -one, H_km1.values(), H_km1.stride(), 
		 F_Z.values(), F_Z.stride(),
		 zero, H_km1_k.values(), H_km1_k.stride());
      // H_km1_k += F_Q * E (result is m-1 by s)
      blas.GEMM (NO_TRANS, NO_TRANS, m-1, s, s+1,
		 one, F_Q.values(), F.stride(), 
		 E.values(), E.stride(),
		 one, H_km1_k.values(), H_km1_k.stride());
    } 


  private:
    /// \brief Candidate basis length for next outer iteration
    ///
    /// Number of "inner" iterations to attempt to execute in the next
    /// "outer" iteration of CA-GMRES; the "s" candidate basis length
    /// parameter.  CA-GMRES may change (either increase or decrease)
    /// this adaptively for numerical stability (decrease) or
    /// performance (usually increase, though not necessarily).
    int candidateBasisLength_;

    //! Length of the most recently generated candidate basis
    int lastCandidateBasisLength_;

    //! Rank of the most recently generated candidate basis
    int lastCandidateBasisRank_;

    /// \brief Number of "outer" iterations completed thus far
    ///
    /// Inclusive range: [0, maxNumOuterIters_].  Incremented at the
    /// end of each successful outer iteration, so a value of
    /// maxNumOuterIters means that all maxNumOuterIters_ outer
    /// iterations were successful.
    int numOuterIters_;

    //! Maximum number of "outer" iterations before restart
    int maxNumOuterIters_;

    /// \brief The number of inner iterations per outer iteration
    ///
    /// The number of inner iterations in each completed outer
    /// iteration.  numInnerIters_[k] more basis vectors were
    /// generated by successful outer iteration k.  This may differ
    /// from candidateBasisLength_, since CA-GMRES is allowed to
    /// change (either increase or decrease) the latter adaptively.
    /// Does not include the initial basis vector (which is the scaled
    /// initial residual vector).
    std::vector<int> numInnerIters_;

    /// \brief Matrix powers kernel implementation
    ///
    /// Pointer to the matrix powers kernel implementation.  This is
    /// nonconst, because the matrix powers kernel may make use of
    /// eigenvalue approximations collected during the solve.  As the
    /// solve progresses, these approximations may get better and
    /// better.
    ///
    /// \note Ask this object for the current change-of-basis matrix.
    ///
    /// \note OrthoManager normalizes "in place," whereas
    /// TsqrOrthoManager uses internal scratch space.  Thus, we
    /// compute the matrix powers kernel in place as well, using the
    /// storage V_ inherited from GmresBase.
    ///
    /// \note The Flexible GMRES (FGMRES) option forces the matrix
    /// powers kernel to use a conventional implementation, rather
    /// than an optimized matrix powers kernel implementation.
    /// Besides this, CA-GMRES works as usual.  If the varying
    /// preconditioner tends to produce a numerically rank-deficient
    /// candidate basis, this will restrict the candidate basis
    /// length, however.
    Teuchos::RCP<Akx<Scalar, MV> > akx_;
  };
} // namespace Belos

#endif //  __Belos_CaGmres_hpp
