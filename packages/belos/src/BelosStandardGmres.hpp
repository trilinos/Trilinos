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

#ifndef __Belos_StandardGmres_hpp
#define __Belos_StandardGmres_hpp

/// \file BelosStandardGmres.hpp
/// \brief Implementation of standard (and Flexible) GMRES
/// \author Mark Hoemmen

#include <BelosGmresBase.hpp>

namespace Belos {

  /// \class StandardGmres
  /// \brief Implementation of (Standard) GMRES
  /// \author Mark Hoemmen
  ///
  /// Implementation of the Generalized Minimal Residual (GMRES)
  /// method of Saad and Schultz (1986) for solving nonsymmetric
  /// linear systems.
  ///
  /// References:
  /// - Saad and Schultz, "GMRES: A generalized minimal residual
  ///   algorithm for solving nonsymmetric linear systems", SISSC,
  ///   vol. 7, pp. 856-869, 1986.
  template<class Scalar, class MV, class OP>
  class StandardGmres : public GmresBase<Scalar,MV,OP> {
  public:
    typedef Scalar scalar_type;
    typedef typename Teuchos::ScalarTraits<Scalar>::magnitudeType magnitude_type;
    typedef MV multivector_type;
    typedef OP operator_type;

  private:
    /// \typedef base_type
    /// \brief Base class typedef
    typedef GmresBase<Scalar, MV, OP> base_type;
    typedef Teuchos::ScalarTraits<Scalar> STS;
    typedef MultiVecTraits<Scalar, MV> MVT;
    typedef Teuchos::SerialDenseMatrix<int, Scalar> mat_type;
    
  public:
    /// Constructor
    ///
    /// \param problem [in/out] Linear problem.  On input, we use the
    ///   starting guess (x0) and the initial residual (r0) to
    ///   initialize the iteration.  The iteration may call
    ///   updateSolution().  On output, if the solution has been
    ///   updated, the vector returned by getLHS() will be modified.
    /// \param ortho [in] Orthogonalization manager
    /// \param outMan [in/out] Output manager
    /// \param maxIterCount [in] Maximum number of iterations before
    ///   restart.  The number of vectors' worth of storage this
    ///   constructor allocates is proportional to this, so choose
    ///   carefully.
    /// \param flexible [in] Whether or not to run the Flexible
    ///   variant of GMRES (FGMRES).  This requires twice as much
    ///   vector storage, but lets the preconditioner change in every
    ///   iteration.  This only works with right preconditioning, not
    ///   left or split preconditioning.
    StandardGmres (const Teuchos::RCP<LinearProblem<Scalar,MV,OP> >& lp,
		   const Teuchos::RCP<const OrthoManager<Scalar, MV> >& ortho,
		   const Teuchos::RCP<OutputManager<Scalar> >& outMan,
		   const int maxIterCount,
		   const bool flexible) :
      GmresBase<Scalar, MV, OP> (lp, ortho, outMan, maxIterCount, flexible) {}

    virtual bool canExtendBasis() const {
      return this->getNumIters() < this->maxNumIters();
    }

    virtual void
    extendBasis (Teuchos::RCP<MV>& V_cur, 
		 Teuchos::RCP<MV>& Z_cur)
    { 
      using Teuchos::Range1D;
      using std::endl;
      const bool verboseDebug = false;
      //
      // mfh 16 Feb 2011: The use of "this->..." here and elsewhere to
      // refer to GmresBase member data is obligatory, since we are
      // inheriting from a templated class.  See the C++ FAQ:
      //
      // http://www.parashift.com/c++-faq-lite/templates.html#faq-35.19
      // 
      RCP<LinearProblem<Scalar, MV, OP> > lp = this->lp_;
      RCP<MV> V = this->V_;
      RCP<MV> Z = this->Z_;
      // This does not count the initial basis vector.
      const int k = this->getNumIters(); 
      const int m = this->maxNumIters();
      TEST_FOR_EXCEPTION(k >= m, GmresCantExtendBasis,
			 "Belos::StandardGmres::extendBasis: "
			 "Maximum number of iterations " << m << "reached; "
			 "cannot extend basis further.");
      std::ostream& dbg = this->outMan_->stream(Debug);
      if (verboseDebug)
	{
	  dbg << "---- StandardGmres::extendBasis: "
	    "V_prv = V[" << k << "," << k << "], "
	    "V_cur = V[" << k+1 << "," << k+1 << "]";
	  if (this->flexible_)
	    dbg << ", Z_cur = Z[" << k << "," << k << "]" << endl;
	  else
	    dbg << endl;
	}
      RCP<const MV> V_prv = MVT::CloneView(*V, Range1D(k, k));
      V_cur = MVT::CloneViewNonConst(*V, Range1D(k+1, k+1));
      if (this->flexible_)
	{
	  Z_cur = MVT::CloneViewNonConst(*Z, Range1D(k, k));
	  lp->applyOp (*V_prv, *Z_cur);
	  lp->applyRightPrec (*Z_cur, *V_cur);
	}
      else
	{
	  lp->apply (*V_prv, *V_cur);
	  if (! Z_cur.is_null())
	    Z_cur = Teuchos::null;
	}
    }

    virtual void
    orthogonalize (const Teuchos::RCP<MV>& V_cur,
		   const Teuchos::RCP<const MV>& V_prv,
		   Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& C_V,
		   Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& B_V,
		   const Teuchos::RCP<MV>& Z_cur,
		   const Teuchos::RCP<const MV>& Z_prv,
		   Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& C_Z,
		   Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& B_Z)
    {
      using Teuchos::null;
      using Teuchos::rcp;
      using Teuchos::tuple;
      const bool verboseDebug = true;

      // Standard GMRES doesn't generate projection and normalization
      // coefficients for the Z basis, regardless of whether we are
      // performing Flexible GMRES.
      if (! C_Z.is_null())
	C_Z = null;
      if (! B_Z.is_null())
	B_Z = null;
      // Help C++'s type inference by referring to the "this" pointer.
      mat_type& H = *(this->H_);
      mat_type& z = *(this->z_);
      const int k = this->getNumIters();

      // Flexible standard GMRES only needs to orthogonalize the new V
      // basis vector; it doesn't need to orthogonalize the Z basis
      // vector.
      //
      // We don't need to do anything with the "rank" output of
      // projectAndNormalize(), since the single value in B will tell
      // us (or rather, the caller) what we (rather, they) need to
      // know.
      C_V = rcp (new mat_type (Teuchos::View, H, k+1, 1, 0, k));
      B_V = rcp (new mat_type (Teuchos::View, H, 1, 1, k+1, k));
      (void) this->ortho_->projectAndNormalize (*V_cur, tuple(C_V), 
						B_V, tuple(V_prv));
      if (verboseDebug)
	{
	  using std::endl;
	  std::ostream& dbg = this->outMan_->stream(Debug);
	  dbg << "---- Current upper Hessenberg matrix (" 
	      << k+2 << " x " << k+1 << "): "
	      << mat_type(Teuchos::View, H, k+2, k+1) << endl;
	  dbg << "---- Current projected right-hand side ("
	      << k+2 << " x " << 1 << "): "
	      << mat_type(Teuchos::View, z, k+2, 1) << endl;
	}
    }

    virtual bool 
    acceptedCandidateBasis () const 
    {
      mat_type& H = *(this->H_);

      // OK, because we haven't accepted the candidate basis yet,
      // so the number of iterations hasn't yet been incremented.
      const int k = this->getNumIters();
      const Scalar H_kp1k = H(k+1, k);
      
      // NOTE (mfh {15,16} Jan 2011) This test should perhaps be more
      // sophisticated.  However, perhaps the right place for such a
      // test is the status check, rather than here.  Certainly if
      // H_(k+1,k) is zero, the basis cannot be extended, even in
      // exact arithmetic.
      return STS::magnitude (H_kp1k) > STS::zero();
    }

    virtual void 
    updateUpperHessenbergMatrix (const Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& C_V,
				 const Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& B_V,
				 const Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& C_Z,
				 const Teuchos::RCP<Teuchos::SerialDenseMatrix<int,Scalar> >& B_Z)
    {
      // Standard GMRES just writes to the upper Hessenberg matrix in
      // place in its implementation of orthogonalize(), so we don't
      // need to do anything here.  Just include the usual boilerplate
      // to avoid compiler warnings for unused inputs.
      (void) C_V;
      (void) B_V;
      (void) C_Z;
      (void) B_Z;
    }

  };

} // namespace Belos

#endif // __Belos_StandardGmres_hpp
