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

#ifndef RTOPPACK_UNARY_FUNC_PTR_DECL_HPP
#define RTOPPACK_UNARY_FUNC_PTR_DECL_HPP

#include "RTOpPack_RTOpT.hpp"

namespace RTOpPack {

/** \brief <tt>RTOpT</tt> subclass for unary transformation functions using a
 * function pointer.
 *
 * ToDo: Finish documentation!
 *
 * \ingroup RTOpPack_RTOpTHelpers_grp
 */
template<class Scalar>
class TOpUnaryFuncPtr : public RTOpT<Scalar> {
public:

  /** \brief . */
  typedef void (*unary_func_ptr_t) ( const Scalar x[], int x_dim, Scalar out[] );

  /// Construct to uninitialized
  TOpUnaryFuncPtr();

  /// Calls <tt>initialize()</tt>
  TOpUnaryFuncPtr(
    unary_func_ptr_t        unary_func_ptr
    ,const std::string      &op_name
    );

  /** \brief Initialize.
   *
   * @param  unary_func_ptr
   *               [in] Pointer to function that actually performs the unary operation.
   * @param  op_name
   *               [in] Name of the operation (for debugging mostly by clients)
   *
   * Preconditions:<ul>
   * <li> <tt>unary_func_ptr != NULL</tt> (throw <tt>std::invalid_argument</tt>)
   * </ul>
   */
  void initialize(
    unary_func_ptr_t        unary_func_ptr
    ,const std::string      &op_name
    );

  /** \brief Set uninitialized.
   *
   * @param  unary_func_ptr
   *               [out] If <tt>unary_func_ptr!=NULL</tt> then <tt>*unary_func_ptr</tt>
   *               is set to pointer to function that was passed in to <tt>initialize()</tt>.
   * @param  op_name
   *               [out] If <tt>op_name!=NULL</tt> then <tt>*op_name</tt>
   *               is set to the operation name that was passed in to <tt>initialize()</tt>.
   */
  void set_initialized(
    unary_func_ptr_t    *unary_func_ptr  = NULL
    ,std::string        *op_name         = NULL
    );

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  const char* op_name() const;
  /** \brief . */
  void apply_op(
    const int   num_vecs,       const ConstSubVectorView<Scalar>         sub_vecs[]
    ,const int  num_targ_vecs,  const SubVectorView<Scalar>  targ_sub_vecs[]
    ,ReductTarget *reduct_obj
    ) const;

  //@}

private:
  
  std::string        op_name_;
  unary_func_ptr_t   unary_func_ptr_;

  // Not defined and not to be called
  TOpUnaryFuncPtr(const TOpUnaryFuncPtr&);
  TOpUnaryFuncPtr& operator=(const TOpUnaryFuncPtr&);

};

} // end namespace RTOpPack

#endif // RTOPPACK_UNARY_FUNC_PTR_DECL_HPP
