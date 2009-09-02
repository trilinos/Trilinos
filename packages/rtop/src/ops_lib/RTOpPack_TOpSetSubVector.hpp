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
