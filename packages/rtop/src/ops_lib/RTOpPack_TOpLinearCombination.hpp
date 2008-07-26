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
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_LINEAR_COMBINATION_HPP
#define RTOPPACK_TOP_LINEAR_COMBINATION_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


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


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANIATION
#  include "RTOpPack_TOpLinearCombination_def.hpp"
#endif


#endif // RTOPPACK_TOP_LINEAR_COMBINATION_HPP 
