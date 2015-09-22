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

#ifndef RTOPPACK_RTOP_SUB_RANGE_DECORATOR_DECL_HPP
#define RTOPPACK_RTOP_SUB_RANGE_DECORATOR_DECL_HPP


#include "RTOpPack_RTOpT.hpp"
#include "Teuchos_ConstNonconstObjectContainer.hpp"


namespace RTOpPack {


/** \brief Decorator subclass that restricts the range of elements to apply
 * the underlying RTOpT object to.
 *
 * This standard decorator class can wrap any <tt>RTOpT</tt> object and
 * restrict what elements in the vector that <tt>%RTOpT</tt> will be applied.
 * The ordinals <tt>first_ele_offset</tt> and <tt>sub_dim</tt> determine which
 * of the global elements the logical vector <tt>v</tt> that the RTOp will be
 * applied.  The ordinal <tt>global_offset_shift</tt> is used to shift adjust
 * the vector element indicies before passing it to the underlying
 * <tt>%RTOpT</tt> object.
 *
 * Therefore, the logical vector <tt>x</tt> derived from the input local
 * vector <tt>v</tt> passed to <tt>getOp()->apply_op(...)</tt> in subvector
 * checks is:

 \verbatim

  z(k + global_offset_shift) = v(k),
    for k = first_ele_offset ... first_ele_offset+sub_dim-1
 \endverbatim

 * For example, if <tt>first_ele_offset=10, sub_dim=5,
 * global_offset_shift=0</tt>, the subvectors elements with indices <tt>k= 10,
 * 11, 12, 13, 14</tt> would be pass into the underling <tt>%RTOpT</tt> object
 * <tt>*this->getOp()</tt> to be processed and would be given the global
 * indices 

 *
 * ToDo: Finish documentation!
 */
template<class Scalar>
class RTOpSubRangeDecorator : public RTOpT<Scalar> {
public:

  /** \name Public types. */
  //@{

  /** \brief . */
  typedef typename RTOpT<Scalar>::primitive_value_type primitive_value_type;
  
  //@}

  /** \name Constructors, accessors. */
  //@{

  /** \brief . */
  RTOpSubRangeDecorator();

  /** \brief . */
  RTOpSubRangeDecorator(
    const RCP<RTOpT<Scalar> > &op,
    const Ordinal first_ele_offset = 0,
    const Ordinal sub_dim = -1
    );

  /** \brief . */
  RTOpSubRangeDecorator(
    const RCP<const RTOpT<Scalar> > &op,
    const Ordinal first_ele_offset = 0,
    const Ordinal sub_dim = -1
    );

  /** \brief . */
  void nonconstInitialize(
    const RCP<RTOpT<Scalar> > &op,
    const Ordinal first_ele_offset = 0,
    const Ordinal sub_dim = -1
    );

  /** \brief . */
  void initialize(
    const RCP<const RTOpT<Scalar> > &op,
    const Ordinal first_ele_offset = 0,
    const Ordinal sub_dim = -1
    );

  /** \brief . */
  RCP<RTOpT<Scalar> > getNonconstOp();

  /** \brief . */
  RCP<const RTOpT<Scalar> > getOp() const;

  //@}

  /** @name Overridden from RTOpT */
  //@{

  /** \brief . */
  void get_reduct_type_num_entries_impl(
    const Ptr<int> &num_values,
    const Ptr<int> &num_indexes,
    const Ptr<int> &num_chars
    ) const;
  /** \brief . */
  Teuchos::RCP<ReductTarget> reduct_obj_create_impl() const;
  /** \brief . */
  void reduce_reduct_objs_impl( const ReductTarget &in_reduct_obj,
    const Ptr<ReductTarget> &inout_reduct_obj) const;
  /** \brief . */
  void reduct_obj_reinit_impl( const Ptr<ReductTarget> &reduct_obj ) const;
  /** \brief . */
  void extract_reduct_obj_state_impl(
    const ReductTarget &reduct_obj,
    const ArrayView<primitive_value_type> &value_data,
    const ArrayView<index_type> &index_data,
    const ArrayView<char_type> &char_data
    ) const;
  /** \brief . */
  void load_reduct_obj_state_impl(
    const ArrayView<const primitive_value_type> &value_data,
    const ArrayView<const index_type> &index_data,
    const ArrayView<const char_type> &char_data,
    const Ptr<ReductTarget> &reduct_obj
    ) const;
  /** \brief . */
  std::string op_name_impl() const;
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

  Teuchos::ConstNonconstObjectContainer<RTOpT<Scalar> > op_;
  Ordinal first_ele_offset_;
  Ordinal sub_dim_;

};


} // namespace RTOpPack


#endif // RTOPPACK_RTOP_SUB_RANGE_DECORATOR_DECL_HPP
