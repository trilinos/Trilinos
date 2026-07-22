// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER


#ifndef RTOPPACK_RTOP_NEW_T_DECL_HPP
#define RTOPPACK_RTOP_NEW_T_DECL_HPP


#include "RTOpPack_Types.hpp"
#include "Teuchos_Describable.hpp"


namespace RTOpPack {


/** \brief Abstract base class for all reduction objects.
 */
class ReductTarget : public Teuchos::Describable
{};


/** \brief Templated interface to vector reduction/transformation operators
 * {abstract}.
 *
 * The purpose of this base class is to allow users to specify
 * arbitrary reduction/transformation operations on vectors without
 * requiring the vectors to reveal their implementation details.  The
 * design is motivated partly by the "Visitor" patter (Gamma, 1995).
 *
 * This interface is designed to allow implementation of a distributed
 * parallel abstract numerical algorithm without the explicit knowledge of the
 * algorithm.
 *
 * In the following discussion, <tt>v[k]</tt>, <tt>x</tt>, <tt>y</tt> and
 * <tt>z</tt> are some abstract vector objects of dimension <tt>n</tt>.  Users
 * can define operators to perform reduction and/or transformation operations.
 * Reduction operations applied over all of the elements of a vector require
 * communication between nodes in a parallel environment but do not change any
 * of the vectors involved.  Transformation operations don't require
 * communication between nodes in a parallel environment.  The targets of a
 * transformation operation is a set of one or more vectors which are changed
 * in some way.
 *
 * The idea is that the user may want to perform reduction operations of the
 * form:
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
 * Operators can also be defined that perform reduction and transformation
 * operations on the same vectors that that should only be done for efficiency
 * reasons.
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
 * Transformation operations are much simpler to think about and to deal with.
 * Some off-the-wall examples of transformation operations that this design
 * will support are:
 *
 \verbatim

 max{ |x(i)|, |y(i)| } + |z(i)| -> z(i), for i = 0...n-1
 
 alpha * |z(i)| / x(i) -> z(i), for i = 0...n-1
 
 alpha * x(i) * y(i) + beta * z(i) -> z(i), for i = 0...n-1
 \endverbatim
 *
 * Reduction operations present the more difficult technical challenge since
 * they require information gathered from all of the elements to arrive at the
 * final result.  This design allows operator classes to be defined that can
 * simultaneously perform reduction and transformation operations:
 *
 \verbatim

   op(v[0](i)...v[*](i),z[0](i)...z[*](i)) -> z[0](i)...z[*](i),reduct_obj,
      for i = 0...n-1
 \endverbatim
 *
 * This design is based on a few assumptions about the reduction and
 * transformation operations and vector implementations themselves.  First, we
 * will assume that vectors are stored and manipulated as chunks of
 * sub-vectors (of dimension <tt>subDim</tt>) where each sub-vector is
 * sufficiently large to overcome the inherent overhead of this design.  This
 * design supports dense strided sub-vectors (see <tt>ConstSubVectorView</tt>
 * and <tt>SubVectorView</tt>) but is relatively flexible.
 *
 * It is strictly the responsibility of the vector implementations to
 * determine how these operators are applied.  For instance, if we are
 * performing a transformation operation of the form:
 *
 \verbatim

 op( x(i), y(i), z(i) ) -> z(i), for i = 0...n-1
 \endverbatim
 *
 * where <tt>x</tt>, <tt>y</tt>, and <tt>z</tt> are distributed parallel
 * vectors, then we would assume that the elements would be partitioned onto
 * the various processors with the same local elements stored on each
 * processor so as not to require any communication between processors.
 */
template<class Scalar>
class RTOpT : public Teuchos::Describable {
public:

  /** @name public types */
  //@{

  /** \brief . */
  typedef typename PrimitiveTypeTraits<Scalar,Scalar>::primitiveType
  primitive_value_type;

  //@}

  /** @name Reduction object functions (NVI) */
  //@{

  /** \brief Get the number of entries of each basic data type in the externalized state for
   * a reduction object for this operator.
   *
   * Note that a specific reduction object is not passed in as an
   * argument. This is because the structure of a reduction object is
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
  void get_reduct_type_num_entries(
    const Ptr<int> &num_values,
    const Ptr<int> &num_indexes,
    const Ptr<int> &num_chars
    ) const
    {
      get_reduct_type_num_entries_impl(num_values, num_indexes, num_chars);
    }

  /** \brief Creates a new reduction target object initialized and ready to be used in
   * a reduction.
   *
   * The default implementation returns <tt>returnVal.get()==NULL</tt>
   * (i.e. by default there is no reduction operation performed).
   */
  Teuchos::RCP<ReductTarget> reduct_obj_create() const
    {
      return reduct_obj_create_impl();
    }

  /** \brief Reduce intermediate reduction target objects.
   *
   * The default implementation does not do anything (i.e. by default
   * there is no reduction operation performed).
   */
  void reduce_reduct_objs(
    const ReductTarget& in_reduct_obj, const Ptr<ReductTarget>& inout_reduct_obj
    ) const
    {
      reduce_reduct_objs_impl( in_reduct_obj, inout_reduct_obj );
    }

  /** \brief Reinitialize an already created reduction object.
   *
   * The default implementation does nothing (i.e. by default there is
   * no reduction operation performed).
   *
   * \param reduct_obj [in/out] Reduction object is reinitialized on output.
   */
  void reduct_obj_reinit( const Ptr<ReductTarget> &reduct_obj ) const
    {
      reduct_obj_reinit_impl(reduct_obj);
    }

  /** \brief Extract the state of an already created reduction object.
   *
   * This method allows the state of a reduction target object to be
   * externalized so that it can be passed over a heterogeneous
   * network of computers.
   *
   * The default implementation does nothing (i.e. by default there is
   * no reduction operation performed).
   */
  void extract_reduct_obj_state(
    const ReductTarget &reduct_obj,
    const ArrayView<primitive_value_type> &value_data,
    const ArrayView<index_type> &index_data,
    const ArrayView<char_type> &char_data
    ) const
    {
      extract_reduct_obj_state_impl( reduct_obj,
        value_data, index_data, char_data );
    }

  /** \brief Load the state of an already created reduction object given
   * arrays of primitive objects.
   *
   * The default implementation does nothing.
   */
  void load_reduct_obj_state(
    const ArrayView<const primitive_value_type> &value_data,
    const ArrayView<const index_type> &index_data,
    const ArrayView<const char_type> &char_data,
    const Ptr<ReductTarget> &reduct_obj
    ) const
    {
      load_reduct_obj_state_impl( value_data, index_data, char_data, reduct_obj );
    }

  //@}

  /** @name Operator functions (NIV) */
  //@{

  /** \brief Return the name (as a null-terminated C-style string) of the operator.
   *
   * This name is used to differentiate an operator subclass from all
   * other operator subclasses. This is an important property needed
   * for a client/server or master/slave runtime configuration.
   *
   * The default implementation uses the value created in the
   * constructor <tt>RTOpT()</tt>.
   */
  std::string op_name() const
    {
      return op_name_impl();
    }

  // 2007/11/14: rabartl: ToDo: Above: change to return std::string.  Don't
  // bother deprecating the old function since you can't really do it very
  // well.

  /** \brief Returns <tt>true</tt> if this operator is coordinate invariant.
   *
   * The default implementation returns <tt>true</tt> as most vector
   * operators are coordinate invariant.
   */
  bool coord_invariant() const
    {
      return coord_invariant_impl();
    }

  /** \brief Returns the continuous range of elements that this operator is
   * defined over.
   *
   * Vector client implementations are free to ignore this but they can use
   * this information to optimize rare operators that only interact with a
   * subset of elements.
   *
   * The default implementation return <tt>Range1D()</tt> which means all of
   * the elements.
   */
  Range1D range() const
    {
      return range_impl();
    }

  /** \brief Apply the reduction/transformation operator to a set of
   * sub-vectors.
   *
   * <tt>op(sub_vecs[], targ_sub_vecs[], reduct_obj) -> targ_sub_vecs[], reduct_obj</tt>.
   *
   * This is the bread and butter of the whole design. Through this method, a
   * vector implementation applies a reduction/transformation operator to a
   * set of sub-vectors.
   *
   * \param sub_vecs [in] Array (length <tt>num_vecs</tt>) of non-mutable
   * vectors to apply the operator over. The ordering of these sub-vectors
   * <tt>sub_vecs[k], for k = 0...num_vecs-1</tt>, is significant to the this
   * operator object. If <tt>num_vecs==0</tt> then <tt>sub_vecs</tt> can be
   * <tt>NULL</tt>.
   *
   * \param targ_sub_vecs [in/out] Array (length <tt>num_targ_vecs</tt>) of
   * mutable vectors to apply the operator over and be mutated. The ordering
   * of these sub-vectors <tt>targ_sub_vecs[k], for k =
   * 0...num_targ_vecs-1</tt>, is significant to this operator object. If
   * <tt>num_targ_vecs==0</tt> then <tt>targ_sub_vecs</tt> can be
   * <tt>NULL</tt>.
   *
   * \param reduct_obj [in/out] This reduction object must have been created
   * by a <tt>this->reduct_obj_create()</tt> call and it may have * already
   * passed through one or more other reduction operations (accumulating the
   * reductions along the way). If
   * <tt>this->get_reduct_type_num_entries()</tt> returns
   * <tt>num_values==0</tt>, <tt>num_indexes==0</tt> and
   * <tt>num_chars==0</tt>, then <tt>reduct_obj</tt> should be set to
   * <tt>NULL</tt> as no reduction will be performed.
   *
   * Preconditions:<ul>
   *
   * <li> <tt>globalOffset==sub_vecs[k].globalOffset</tt> , for <tt>k =
   * 0,...,sub_vecs.subDim()-1</tt>
   *
   * <li> <tt>globalOffset==targ_sub_vecs[k].globalOffset</tt> , for <tt>k =
   * 0,...,targ_vecs.subDim()-1</tt>
   *
   * <li> <tt>subDim==sub_vecs[k].subDim()</tt> , for <tt>k =
   * 0,...,sub_vecs.subDim()-1</tt>
   *
   * <li> <tt>subDim==targ_sub_vecs[k].subDim()</tt> , for <tt>k =
   * 0,...,targ_vecs.subDim()-1</tt>
   *
   * </ul>
   *
   * If <tt>nonnull(reduct_obj)==true</tt> then the reduction operation will
   * be accumulated as:
   *
   \verbatim
   op(sub_vecs[], targ_sub_vecs[], reduct_obj) -> reduct_obj
   \endverbatim
   *
   * By allowing an in/out <tt>reduct_obj</tt> and an accumulation of
   * the reduction, the maximum reuse of memory is achieved. If
   * <tt>this->reduct_obj_create()</tt> or
   * <tt>this->reduct_obj_reinit()</tt> (passing in
   * <tt>reduct_obj</tt>) was called immediately before this function,
   * then on return, <tt>reduct_obj</tt> will contain only the
   * reduction from this function call.
   *
   * If the sizes of <tt>sub_vecs</tt> and <tt>targ_sub_vecs</tt> is
   * incompatible with the underlying operator object then
   * <tt>InvalidNumVecs</tt> is thrown. If the sub-vectors are not compatible
   * (i.e. <tt>globalOffset</tt> and/or <tt>subDim</tt> not the same) then
   * <tt>IncompatibleVecs</tt> is thrown.
   */
  void apply_op(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj
    ) const
    {
      apply_op_impl(sub_vecs, targ_sub_vecs, reduct_obj);
    }

  //@}

protected:

  /** \name Protected virtual functions to be overridden by subclasses. */
  //@{

  /** \brief . */
  virtual void get_reduct_type_num_entries_impl(
    const Ptr<int> &num_values,
    const Ptr<int> &num_indexes,
    const Ptr<int> &num_chars
    ) const;

  /** \brief . */
  virtual Teuchos::RCP<ReductTarget> reduct_obj_create_impl() const;

  /** \brief . */
  virtual void reduce_reduct_objs_impl(
    const ReductTarget& in_reduct_obj, const Ptr<ReductTarget>& inout_reduct_obj
    ) const;

  /** \brief . */
  virtual void reduct_obj_reinit_impl( const Ptr<ReductTarget> &reduct_obj ) const;

  /** \brief . */
  virtual void extract_reduct_obj_state_impl(
    const ReductTarget &reduct_obj,
    const ArrayView<primitive_value_type> &value_data,
    const ArrayView<index_type> &index_data,
    const ArrayView<char_type> &char_data
    ) const;

  /** \brief . */
  virtual void load_reduct_obj_state_impl(
    const ArrayView<const primitive_value_type> &value_data,
    const ArrayView<const index_type> &index_data,
    const ArrayView<const char_type> &char_data,
    const Ptr<ReductTarget> &reduct_obj
    ) const;

  /** \brief . */
  virtual std::string op_name_impl() const;

  /** \brief . */
  virtual bool coord_invariant_impl() const;

  /** \brief . */
  virtual Range1D range_impl() const;

  /** \brief . */
  virtual void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj
    ) const = 0;

  //@}

  /** \name Nonvirtual protected functions. */
  //@{

  /** \brief Constructor that creates an operator name appended with the
   * type. */
  RTOpT( const std::string &op_name_base = "" );

  /** \brief Just set the operator name. */
  void setOpNameBase( const std::string &op_name_base );

  //@}

public:


private:

  std::string op_name_;

  void throwNoReductError() const;

}; // end class RTOpT


// 2007/11/14: rabartl: ToDo: Break off an RTOpDefaultBase interface and put
// all default implementation functions in there.


} // end namespace RTOpPack


#endif // RTOPPACK_RTOP_NEW_T_DECL_HPP
