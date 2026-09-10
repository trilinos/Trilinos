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

namespace random_impl {

inline std::uint64_t mix64(std::uint64_t x)
{
  x = (x ^ (x >> 30)) * 0xbf58476d1ce4e5b9ULL;
  x = (x ^ (x >> 27)) * 0x94d049bb133111ebULL;
  return x ^ (x >> 31);
}

inline double unit_uniform(
    const std::uint64_t seed,
    const std::uint64_t counter,
    const std::uint64_t global_index)
{
  std::uint64_t x = mix64(seed+0x9e3779b97f4a7c15ULL);
  x = mix64(x + mix64(counter));
  x = mix64(x + mix64(global_index));

  // Shift right by 11 to keep the top 53 high-quality bits
  // and multiply by 2^-53 to produce a double in [0,1).
  return static_cast<double>(x >> 11) * 0x1.0p-53;
}
}


/** \brief Generate a uniform random vector in the range [l,u].
 *
 * The seed for the random number generator can be set by
 * <tt>TOpRandomize<Scalar>::set_static_seed(s)</tt> where <tt>s</tt> is some
 * unsigned integer. 
 * A counter is incremented every time a new object is created. 
 * This class generates random numbers based on the initial seed, the counter, and the global element ID so this produces the same
 * pseudo-random elements independent of the number of processors being used.
 * The counter is set to zero when the static seed is set via <tt>TOpRandomize<Scalar>::set_static_seed(s)</tt>.
 * The implementation is based on the Stateless SplitMix64 finalizer.
 * See:
 * G. L. Steele Jr., D. Lea, C. H. Flood,
 * "Fast Splittable Pseudorandom Number Generators",
 * OOPSLA 2014, DOI: 10.1145/2660193.2660195.
 *
 * When available it might be worth adopting the philox_engine counter-based random generator that should be part of C++26.
 */
template<class Scalar>
class TOpRandomize : public RTOpT<Scalar> {
public:
  using RTOpT<Scalar>::apply_op;
  /** \brief . */
  static void set_static_seed( const unsigned int static_seed )
  { 
    static_seed_ = static_seed;
    static_counter_ = 0;
  }

  /** \brief . */
  static unsigned int get_static_seed() { return static_seed_; }

  /** \brief . */
  TOpRandomize(
    const Scalar& l = -ScalarTraits<Scalar>::one(),
    const Scalar& u = +ScalarTraits<Scalar>::one()
    ) : seed_(static_seed_), counter_(static_counter_++), l_(l), u_(u)
    {
      this->setOpNameBase("TOpRandomize");
    }

  /** \brief . */
  void set_bounds( const Scalar& l, const Scalar& u ) { l_ = l; u_ = u; }

  /** \brief . */
  void set_seed( const unsigned int seed ) { seed_ = seed; }

  /** \brief . */
  void set_counter( const unsigned int counter ) { counter_ = counter; }

  /** \brief . */
  unsigned int get_seed() const { return seed_; }

  /** \brief . */
  unsigned int get_counter() const { return counter_; }

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
        auto rand = random_impl::unit_uniform(
                     static_cast<std::uint64_t>(seed_),
                     static_cast<std::uint64_t>(counter_),
                     static_cast<std::uint64_t>(globalOffset+i)
                    );
        *z0_val = a * static_cast<Scalar>(rand) + b;
        // Above should be in the range [l,b]
      }
    }
  //@}
private:
  static unsigned int static_seed_;
  static unsigned int static_counter_;
  unsigned int seed_;
  unsigned int counter_;
  Scalar l_;
  Scalar u_;
};


template<class Scalar>
unsigned int TOpRandomize<Scalar>::static_seed_ = 0;

template<class Scalar>
unsigned int TOpRandomize<Scalar>::static_counter_ = 0;


} // namespace RTOpPack


#endif // RTOPPACK_TOP_RANDOMIZE_HPP
