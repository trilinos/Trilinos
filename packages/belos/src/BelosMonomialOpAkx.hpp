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

#ifndef __Belos_MonomialOpAkx_hpp
#define __Belos_MonomialOpAkx_hpp

#include <BelosOpAkx.hpp>

namespace Belos {

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

#endif // __Belos_MonomialOpAkx_hpp
