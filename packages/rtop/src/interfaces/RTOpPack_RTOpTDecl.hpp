// @HEADER
// ***********************************************************************
// 
//      Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms 
//                 Copyright (2004) Sandia Corporation
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
// Questions? Contact Michael A. Heroux (maherou@sandia.gov)
// 
// ***********************************************************************
// @HEADER

#ifndef RTOPPACK_RTOP_NEW_T_DECL_HPP
#define RTOPPACK_RTOP_NEW_T_DECL_HPP

#include "RTOpPack_Types.hpp"
#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_PrimitiveTypeTraits.hpp"
#include "Teuchos_Describable.hpp"

namespace RTOpPack {

/** \defgroup RTOpCpp_grp Templated interfaces for generalized vector reduction/transformation operators in C++.
 */
//@{

/** \brief Abstract base class for all reduction objects.
 */
class ReductTarget : public Teuchos::Describable
{};

/** \brief Templated interface to vector reduction/transformation operators {abstract}.
 *
 * The purpose of this base class is to allow users to specify
 * arbitrary reduction/transformation operations on vectors without
 * requiring the vectors to reveal their implementation details.  The
 * design is motivated partly by the "Visitor" patter (Gamma, 1995).
 *
 * This interface is designed to allow implementation of a distributed
 * parallel application without the explicit knowledge by the
 * application.
 *
 * In the following discussion, <tt>v[k]</tt>, <tt>x</tt>, <tt>y</tt>
 * and <tt>z</tt> are some abstract vector objects of dimension
 * <tt>n</tt>.  Users can define operators to perform reduction and/or
 * transformation operations.  Reduction operations applied over all
 * of the elements of a vector require communication between nodes in
 * a parallel environment but do not change any of the vectors
 * involved.  Transformation operations don't require communication
 * between nodes in a parallel environment.  The targets of a
 * transformation operation is a set of one or more vectors which are
 * mutated in some way.
 *
 * The idea is that the user may want to perform reduction
 * operations of the form:
 *
 \verbatim

 op(v[0]...v[*],z[0]...z[*]) -> reduct_obj
 \endverbatim
 *
 * where <tt>reduct_obj</tt> is a single object based on a reduction over
 * all the elements of the vector arguments, or transformation
 * operations of the form:
 *
 \verbatim

 op(v[0](i)...v[*](i),z[0](i)...z[*](i)) -> z[0](i)...z[*](i), for i = 0...n-1
 \endverbatim
 *
 * Operators can also be defined that perform reduction and
 * transformation operations.
 *
 * The tricky part though, is that the <tt>reduct_obj</tt> object of
 * the reduction operation may be more complex than a single scalar
 * value.  For instance, it could be a <tt>double</tt> and an
 * <tt>int</tt> pair such as in the reduction operation:
 *
 \verbatim

 min{ |x(i)|, i = 0...n-1 } -> [ x(j_min), j_min ]
 \endverbatim
 *
 * or it could perform several reductions at once and store
 * several scalar values such as in:
 *
 \verbatim

 min_max_sum{ x(i), i = 0...n-1 } -> [ x(j_min), j_min, x(j_max), j_max, x_sum ]
 \endverbatim
 *
 * Transformation operations are much simpler to think about and to
 * deal with.  Some off-the-wall examples of transformation operations
 * that this design will support are:
 *
 \verbatim

 max{ |x(i)|, |y(i)| } + |z(i)| -> z(i), for i = 0...n-1
 
 alpha * |z(i)| / x(i) -> z(i), for i = 0...n-1
 
 alpha * x(i) * y(i) + beta * z(i) -> z(i), for i = 0...n-1
 \endverbatim
 *
 * Reduction operations present the more difficult technical
 * challenge since they require information gathered from all of the
 * elements to arrive at the final result.  This design allows
 * operator classes to be defined that can simultaneously perform
 * reduction and transformation operations:
 *
 \verbatim

   op(v[0](i)...v[*](i),z[0](i)...z[*](i)) -> z[0](i)...z[*](i),reduct_obj
      , for i = 0...n-1
 \endverbatim
 *
 * This design is based on a few assumptions about the reduction and
 * transformation operations and vector implementations themselves.
 * First, we will assume that vectors are stored and manipulated as
 * chunks of sub-vectors (of dimension <tt>subDim</tt>) where each
 * sub-vector is sufficiently large.  This design supports dense
 * strided sub-vectors (see <tt>ConstSubVectorView</tt> and <tt>SubVectorView</tt>)
 * and is relatively flexible.
 *
 * It is strictly the responsibilities of the vector implementations
 * to determine how these operators are applied.  For instance, if we
 * are performing a transformation operation of the form:
 *
 \verbatim

 op( x(i), y(i), z(i) ) -> z(i), for i = 0...n-1
 \endverbatim
 *
 * where <tt>x</tt>, <tt>y</tt>, and <tt>z</tt> are distributed
 * parallel vectors, then we would assume that the elements would be
 * partitioned onto the various processors with the same local
 * elements stored on each processor so as not to require any
 * communication between processors.
 */
template<class Scalar>
class RTOpT : public Teuchos::Describable {
public:

  /** @name public types */
  //@{

  /** \brief . */
  typedef typename Teuchos::PrimitiveTypeTraits<Scalar>::primitiveType  primitive_value_type;

  //@}

  /// Constructor that creates an operator name appended with the type.
  RTOpT( const std::string &op_name_base );

  /** @name Reduction object functions */
  //@{

  /** \brief Get the number of entries of each basic data type in the externalized state for
   * a reduction object for this operator.
   *
   * Note that a specific reduction object is not passed in as an
   * argument.  This is because the structure of a reduction object is
   * completely determined by its associated operator object and this
   * structure can not change as a result of a reduction operation
   * (this is needed to simplify global communication code when used *
   * with MPI).
   *
   * The default implementation returns zeros for
   * <tt>*num_values</tt>, <tt>*num_indexes</tt> and
   * <tt>*num_chars</tt> (i.e. by default there is no reduction
   * operation performed).
   */
  virtual void get_reduct_type_num_entries(
    int*   num_values
    ,int*  num_indexes
    ,int*  num_chars
    ) const;
  /** \brief Creates a new reduction target object initialized and ready to be used in
   * a reduction.
   *
   * To delete this object simply let the returned
   * <tt>RefCountPtr<></tt> object go out of scope.
   *
   *
   * The default implementation returns <tt>return.get()==NULL</tt>
   * (i.e. by default there is no reduction operation performed).
   */
  virtual Teuchos::RefCountPtr<ReductTarget> reduct_obj_create() const;
  /** \brief Reduce intermediate reduction target objects.
   *
   * The default implementation does not do anything (i.e. by default
   * there is no reduction operation performed).
   */
  virtual void reduce_reduct_objs(
    const ReductTarget& in_reduct_obj, ReductTarget* inout_reduct_obj
    ) const;
  /** \brief Reinitialize an already created reduction object.
   *
   * The default implementation does nothing (i.e. by default there is
   * no reduction operation performed).
   *
   * @param  reduct_obj  [in/out] Reduction object is reinitialized on output.
   */
  virtual void reduct_obj_reinit( ReductTarget* reduct_obj ) const;
  /** \brief Extract the state of an already created reduction object.
   *
   * This method allows the state of a reduction target object to be
   * externalized so that it can be passed over a heterogeneous
   * network of computers.
   *
   * The default implementation does nothing (i.e. by default there is
   * no reduction operation performed).
   */
  virtual void extract_reduct_obj_state(
    const ReductTarget        &reduct_obj
    ,int                      num_values
    ,primitive_value_type     value_data[]
    ,int                      num_indexes
    ,index_type               index_data[]
    ,int                      num_chars
    ,char_type                char_data[]
    ) const;
  /** \brief Load the state of an already created reduction object given
   * arrays of primitive objects.
   *
   * The default implementation does nothing.
   */
  virtual void load_reduct_obj_state(
    int                            num_values
    ,const primitive_value_type    value_data[]
    ,int                           num_indexes
    ,const index_type              index_data[]
    ,int                           num_chars
    ,const char_type               char_data[]
    ,ReductTarget                 *reduct_obj
    ) const;

  //@}

  /** @name Operator functions */
  //@{

  /** \brief Return the name (as a null-terminated C-style string) of the operator.
   *
   * This name is used to differentiate an operator subclass from all
   * other operator subclasses.  This is an important property needed
   * for a client/server or master/slave runtime configuration.
   *
   * The default implementation uses the value created in the
   * constructor <tt>RTOpT()</tt>.
   */
  virtual const char* op_name() const;
  /** \brief Copy the state of another operator object into this one.
   *
   * This default version uses <tt>op->ref extract_op_state()</tt> and
   * <tt>this->\ref load_op_state()</tt> to perform the copy.  No
   * override should be needed by subclasses (unless a slightly more
   * efficient implementation is desired).
   */
  virtual RTOpT<Scalar>& operator=(const RTOpT<Scalar>& op);
  /** \brief Return the number of entries of each type of basic data type in
   * the externalized state for the operator object.
   *
   * The default implementation returns zeros for
   * <tt>*num_values</tt>, <tt>*num_indexes</tt> and
   * <tt>*num_chars</tt> (i.e. the default reduction/transformation
   * operator has no state).
   */
  virtual void get_op_type_num_entries(
    int*  num_values
    ,int* num_indexes
    ,int* num_chars
    ) const;
  /** \brief Extract the state of the object in a portable format.
   *
   * This method allows the state of a reduction/transformation
   * operator to be externalized so that it can be passed over a
   * heterogeneous network of computers.
   *
   * The default implementation does nothing (i.e. the default
   * reduction/transformation operator has no state).
   *
   * @param  num_values
   *              [in] Value returned from <tt>this->get_op_type_num_entries()</tt>.
   * @param  value_data
   *              [out] Array (size <tt>num_values</tt>) of floating point data.
   * @param  num_indexes
   *              [in] Value returned from <tt>this->get_op_type_num_entries()</tt>.
   * @param  index_data
   *              [out] Array (size <tt>num_indexes</tt>) of integral data.
   * @param  num_chars
   *              [in] Value returned from <tt>this->get_op_type_num_entries()</tt>.
   * @param  char_data
   *              [out] Array (size <tt>num_chars</tt>) of character data.
   */
  virtual void extract_op_state(
    int                             num_values
    ,primitive_value_type           value_data[]
    ,int                            num_indexes
    ,index_type                     index_data[]
    ,int                            num_chars
    ,char_type                      char_data[]
    ) const;
  /** \brief Load the state of the object from a portable format.
   *
   * This method allows the state of the operator object to be set given an the externalized
   * state as extracted using <tt>\ref extract_op_state "extract_op_state(...)"</tt>
   * called on a compatible operator object (possibly on a different heterogeneous computer).
   *
   * The default implementation does nothing (i.e. the default reduction/transformation
   * operator has no state).
   *
   * @param  num_values
   *              [in] Value returned from <tt>this->get_op_type_num_entries()</tt>.
   * @param  value_data
   *              [out] Array (size <tt>num_values</tt>) of floating point data.
   * @param  num_indexes
   *              [in] Value returned from <tt>this->get_op_type_num_entries()</tt>.
   * @param  index_data
   *              [out] Array (size <tt>num_indexes</tt>) of integral data.
   * @param  num_chars
   *              [in] Value returned from <tt>this->get_op_type_num_entries()</tt>.
   * @param  char_data
   *              [out] Array (size <tt>num_chars</tt>) of character data.
   */
  virtual void load_op_state(
    int                           num_values
    ,const primitive_value_type   value_data[]
    ,int                          num_indexes
    ,const index_type             index_data[]
    ,int                          num_chars
    ,const char_type              char_data[]
    );
  /** \brief Returns <tt>true</tt> if this operator is coordinate invariant.
   *
   * The default implementation returns <tt>true</tt> as most vector
   * operators are coordinate invariant.
   */
  virtual bool coord_invariant() const;
  /** \brief Apply the reduction/transformation operator to a set of
   * sub-vectors.
   *
   * <tt>op(sub_vecs[],targ_sub_vecs[]),reduct_obj) -> targ_sub_vecs[],reduct_obj</tt>.
   *
   * This is the bread and butter of the whole design.  Through this
   * method, a vector implementation applies a
   * reduction/transformation operator to a set of sub-vectors.
   *
   *	@param	num_vecs
   *             [in] Number of non-mutable sub-vectors in <tt>sub_vec[]</tt>.
   *	@param	sub_vecs
   *             [in] Array (length <tt>num_vecs</tt>) of non-mutable vectors to apply the
   *             operator over.  The ordering of these sub-vectors
   *             <tt>sub_vecs[k], for k = 0...num_vecs-1</tt>, is significant to the this operator object.
   *             If <tt>num_vecs==0</tt> then <tt>sub_vecs</tt> can be <tt>NULL</tt>.
   *	@param	num_targ_vecs
   *             [in] Number of mutable sub-vectors in <tt>targ_sub_vec[]</tt>.
   *	@param	targ_sub_vecs
   *             [in/out] Array (length <tt>num_targ_vecs</tt>) of mutable vectors to apply the
   *             operator over and be mutated.  The ordering of these sub-vectors
   *             <tt>targ_sub_vecs[k], for k = 0...num_targ_vecs-1</tt>, is significant to this
   *             operator object.  If <tt>num_targ_vecs==0</tt> then <tt>targ_sub_vecs</tt> can be <tt>NULL</tt>.
   *	@param	reduct_obj
   *             [in/out] This reduction object must have been created by
   *              a <tt>this->reduct_obj_create()</tt> call and it may have
   *              already passed through one or more other
   *              reduction operations (accumulating the reductions along the way).
   *              If <tt>this->get_reduct_type_num_entries()</tt>
   *              returns <tt>num_values==0</tt>, <tt>num_indexes==0</tt> and <tt>num_chars==0</tt>,
   *              then <tt>reduct_obj</tt> should be set to <tt>NULL</tt>
   *              as no reduction will be performed.
   *
   * Preconditions:<ul>
   *	<li> <tt>num_vecs > 0 || num_targ_vecs > 0</tt>
   *	<li> <tt>[num_vecs > 0] sub_vecs != NULL</tt>
   *	<li> <tt>[num_targ_vecs > 0] targ_sub_vecs != NULL</tt>
   *	<li> <tt>[num_vecs > 0] globalOffset==sub_vecs[k].globalOffset</tt>
   *        , for <tt>k = 0,...,num_vecs-1</tt>
   *	<li> <tt>[num_targ_vecs > 0] globalOffset==targ_sub_vecs[k].globalOffset</tt>
   *        , for <tt>k = 0,...,num_targ_vecs-1</tt>
   *	<li> <tt>[num_vecs > 0] sub_dim==sub_vecs[k].subDim()</tt>
   *       , for <tt>k = 0,...,num_vecs-1</tt>
   *	<li> <tt>[num_targ_vecs > 0] sub_dim==targ_sub_vecs[k].subDim()</tt>
   *       , for <tt>k = 0,...,num_targ_vecs-1</tt>
   *	</ul>
   *
   * If <tt>reduct_obj != NULL</tt> then the reduction operation will
   * be accumulated as:
   *
   \verbatim
     op(op(sub_vecs[],targ_sub_vecs[]),reduct_obj) -> reduct_obj
   \endverbatim
   *
   * By allowing an in/out <tt>reduct_obj</tt> and an accumulation of
   * the reduction, the maximum reuse of memory is achieved.  If
   * <tt>this->reduct_obj_create()</tt> or
   * <tt>this->reduct_obj_reinit()</tt> (passing in
   * <tt>reduct_obj</tt>) was called immediately before this function,
   * then on return, <tt>reduct_obj</tt> will contain only the
   * reduction from this function call.
   *
   * If <tt>num_vecs</tt> is incompatible with the underlying operator object then
   * <tt>InvalidNumVecs</tt> is thrown.  If the sub-vectors are not compatible
   * (i.e. <tt>global_offset</tt> and/or <tt>sub_dim</tt> not the same) then
   * IncompatibleVecs is thrown.
   */
  virtual void apply_op(
    const int   num_vecs,       const ConstSubVectorView<Scalar>  sub_vecs[]
    ,const int  num_targ_vecs,  const SubVectorView<Scalar>       targ_sub_vecs[]
    ,ReductTarget *reduct_obj
    ) const = 0;

  //@}

private:

  std::string op_name_;

}; // end class RTOpT

} // end namespace RTOpPack

#endif // RTOPPACK_RTOP_NEW_T_DECL_HPP
