// @HEADER
// *****************************************************************************
// RTOp: Interfaces and Support Software for Vector Reduction Transformation
//       Operations
//
// Copyright 2006 NTESS and the RTOp contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#ifndef RTOPPACK_TOP_RANDOMIZE_HPP
#define RTOPPACK_TOP_RANDOMIZE_HPP

#include "RTOpPack_RTOpTHelpers.hpp"


namespace RTOpPack {


/** \brief Generate a random vector in the range [l,u]: <tt>z0[i] =
 * 0.5*((u-l)*Teuchos::ScalarTraits<Scalar>::random()+(u+l)), i=0...n-1</tt>.
 *
 * The seed for the random number generator can be set by
 * <tt>TOpRandomize<Scalar>::set_seed(s)</tt> where <tt>s</tt> is some
 * unsigned integer.  Note that this class generates random numbers based on
 * the initial seed and the global element ID so this should produce the same
 * pseudo-random elements independent of the number of processors being used.
 *
 * The seed changes every time a new object is created in order to improve the
 * randomness to some degree.
 */
template<class Scalar>
class TOpRandomize : public RTOpT<Scalar> {
public:
  using RTOpT<Scalar>::apply_op;
  /** \brief . */
  static void set_static_seed( const unsigned int static_seed )
    { static_seed_ = static_seed; }
  /** \brief . */
  static unsigned int get_static_seed() { return static_seed_; }
  /** \brief . */
  TOpRandomize(
    const Scalar& l = -ScalarTraits<Scalar>::one(),
    const Scalar& u = +ScalarTraits<Scalar>::one()
    )
    {
      this->setOpNameBase("TOpRandomize");
      set_bounds(l, u);
      set_seed(static_seed_);
      ++static_seed_; // By default we will just increment the seed!
    }
  /** \brief . */
  void set_bounds( const Scalar& l, const Scalar& u )
    { l_ = l; u_ = u; }
  /** \brief . */
  void set_seed( const unsigned int seed ) { seed_ = seed; }
  /** \brief . */
  unsigned int get_seed() const { return seed_; }
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void apply_op_impl(
    const ArrayView<const ConstSubVectorView<Scalar> > &sub_vecs,
    const ArrayView<const SubVectorView<Scalar> > &targ_sub_vecs,
    const Ptr<ReductTarget> &reduct_obj_inout
    ) const
    {
      typedef typename Teuchos::ArrayRCP<Scalar>::iterator iter_t;

#ifdef TEUCHOS_DEBUG
      validate_apply_op<Scalar>(*this, 0, 1, false,
        sub_vecs, targ_sub_vecs, reduct_obj_inout.getConst());
#else
      (void)sub_vecs;
      (void)reduct_obj_inout;
#endif
      
      const index_type subDim = targ_sub_vecs[0].subDim();
      const index_type globalOffset =  targ_sub_vecs[0].globalOffset();

      iter_t z0_val = targ_sub_vecs[0].values().begin();
      const ptrdiff_t z0_s = targ_sub_vecs[0].stride();

      // Linear coefficients for translating from [-1,+1] to [l,b]
      const Scalar a = Scalar(0.5)*(u_ - l_);
      const Scalar b = Scalar(0.5)*(u_ + l_);
      for( index_type i = 0; i < subDim; ++i, z0_val += z0_s )
      {
        Teuchos::ScalarTraits<Scalar>::seedrandom(seed_+globalOffset+i);
        *z0_val = a * Teuchos::ScalarTraits<Scalar>::random() + b;
        // Above should be in the range [l,b]
      }
    }
  //@}
private:
  static unsigned int static_seed_;
  unsigned int seed_;
  Scalar l_;
  Scalar u_;
};


template<class Scalar>
unsigned int TOpRandomize<Scalar>::static_seed_ = 0;


} // namespace RTOpPack


#endif // RTOPPACK_TOP_RANDOMIZE_HPP
