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

#ifndef __Belos_OpAkx_hpp
#define __Belos_OpAkx_hpp

#include <BelosAkx.hpp>
#include <BelosOperatorTraits.hpp>
#include <Teuchos_OrdinalTraits.hpp>
#include <Teuchos_SerialDenseVector.hpp>

namespace Belos {

  /// \class OpAkx
  /// \brief Base for Akx implementations using only abstract operators.
  /// \author Mark Hoemmen
  ///
  /// It's not always beneficial or desirable to use a matrix powers
  /// kernel implementation that expensively reorganizes the sparse
  /// matrix and preconditioners.  Also, in some cases, we might only
  /// have access to the matrix as an abstract operator, or the
  /// operator might not even be represented as a matrix.  In these
  /// cases, subclasses of OpAkx can be useful.  They implement the
  /// matrix powers kernel straightforwardly, using s applications of
  /// the operator (including left and/or right preconditioner(s), if
  /// provided) to generate s basis vectors.  Subclasses implement
  /// specific basis types, like monomial or Newton.
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
	
    //! Maximum candidate basis length.
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

} // namespace Belos

#endif // __Belos_OpAkx_hpp
