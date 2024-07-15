// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_SET_SUB_VECTOR_HPP
#define RTOPPACK_TOP_SET_SUB_VECTOR_HPP

#include "RTOpPack_RTOpTHelpers.hpp"
#include "RTOpPack_SparseSubVectorT.hpp"


namespace RTOpPack {


/** \brief Advanced transformation operator that assigns elements from a
 * sparse explicit vector.
 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class TOpSetSubVector : public RTOpT<Scalar> {
public:
  
  /** \brief . */
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;

  /** \name Constructors/initializers. */
  //@{

  /** \brief . */
  TOpSetSubVector();

  /** \brief . */
  TOpSetSubVector( const SparseSubVectorT<Scalar> &sub_vec );

  /** \brief . */
  void set_sub_vec( const SparseSubVectorT<Scalar> &sub_vec );

  //@}

protected:

  /** \name Overridden protected functions from RTOpT. */
  //@{

  /** \brief . */
  bool coord_invariant_impl() const;

  /** \brief . */
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj
    ) const;

  //@}

private:

  enum { num_sub_vec_members = 6 };

  SparseSubVectorT<Scalar> sub_vec_;

}; // class TOpSetSubVector


} // namespace RTOpPack


#ifndef HAVE_TEUCHOS_EXPLICIT_INSTANIATION
#  include "RTOpPack_TOpSetSubVector_def.hpp"
#endif


#endif // RTOPPACK_TOP_SET_SUB_VECTOR_HPP
