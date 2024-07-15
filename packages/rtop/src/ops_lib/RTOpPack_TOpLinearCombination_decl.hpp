// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
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
