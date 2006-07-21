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
 * the initial seed and the global element ID so this should produce the save
 * elements independent of the number of processors being used.
 *
 * The seed changes every time a new object is created
 */
template<class Scalar>
class TOpRandomize : public ROpScalarScalarTransformationBase<Scalar> {
public:
  /** \brief . */
  static void set_static_seed( const unsigned int static_seed ) { static_seed_ = static_seed; }
  /** \brief . */
  static unsigned int get_static_seed() { return static_seed_; }
  /** \brief . */
  TOpRandomize(
    const Scalar& l   = -Teuchos::ScalarTraits<Scalar>::one()
    ,const Scalar& u  = +Teuchos::ScalarTraits<Scalar>::one()
    )
    :RTOpT<Scalar>("TOpRandomize"), ROpScalarScalarTransformationBase<Scalar>(l,u)
    {
      seed_ = static_seed_;
      ++static_seed_; // By default we will just increment the seed!
    }
  /** \brief . */
  void set_bounds( const Scalar& l, const Scalar& u ) { this->scalarData1(l); this->scalarData2(u); }
  /** \brief . */
  void set_seed( const unsigned int seed ) { seed_ = seed; }
  /** \brief . */
  unsigned int get_seed() const { return seed_; }
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void apply_op(
    const int   num_vecs,       const ConstSubVectorView<Scalar>         sub_vecs[]
    ,const int  num_targ_vecs,  const SubVectorView<Scalar>              targ_sub_vecs[]
    ,ReductTarget *reduct_obj
    ) const
    {
      const Scalar l = this->scalarData1(), u = this->scalarData2();
      const Scalar a = Scalar(0.5)*(u-l), b = Scalar(0.5)*(u+l) ; // Linear coefficients for translating from [-1,+1] to [l,b]
      RTOP_APPLY_OP_0_1(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      for( Teuchos_Index i = 0; i < subDim; ++i, z0_val += z0_s ) {
        Teuchos::ScalarTraits<Scalar>::seedrandom(seed_+globalOffset+i);
        *z0_val = a * Teuchos::ScalarTraits<Scalar>::random() + b; // Should be in the range [l,b]
      }
    }
  //@}
private:
  static unsigned int     static_seed_;
  unsigned int            seed_;
}; // class TOpRandomize

template<class Scalar>
unsigned int TOpRandomize<Scalar>::static_seed_ = 0;

} // namespace RTOpPack

#endif // RTOPPACK_TOP_RANDOMIZE_HPP
