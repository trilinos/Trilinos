// @HEADER
// ***********************************************************************
// 
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//                Copyright (2006) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
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
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_LINEAR_COMBINATION_DECL_HPP
#define RTOPPACK_TOP_LINEAR_COMBINATION_DECL_HPP

#include "RTOpPack_RTOpTHelpers.hpp"
#include "Teuchos_Workspace.hpp"


namespace RTOpPack {


/** \brief Linear combination transformation operator: <tt>z0[i] = beta*z0[i]
 * + sum( alpha[k]*v[k][i], k=0...num_vecs-1 ), i=0...n-1</tt>.
 *
 * This transformation operator only accepts <tt>num_targ_vec==1</tt>
 * but accepts any <tt>num_vecs > 0</tt>.
 *
 * Warning! this class can only be used in SPMD mode and not
 * client/server or master/slave.  You know what needs to happen for
 * this to work!
 */
template<class Scalar>
class TOpLinearCombination : public RTOpT<Scalar> {
public:

  /** \brief . */
  TOpLinearCombination(
    const ArrayView<const Scalar> &alpha_in = Teuchos::null,
    const Scalar &beta = Teuchos::ScalarTraits<Scalar>::zero()
    );

  /** \brief . */
  void alpha( const ArrayView<const Scalar> &alpha_in );

  /** \brief . */
  const ArrayView<const Scalar> alpha() const;

  /** \brief . */
  void beta( const Scalar& beta_in );

  /** \brief . */
  Scalar beta() const;

  /** \brief . */
  int num_vecs() const;

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj_inout
    ) const;

  //@}

private:

  Scalar beta_;
  Array<Scalar> alpha_;

};


} // namespace RTOpPack


#endif // RTOPPACK_TOP_LINEAR_COMBINATION_DECL_HPP
