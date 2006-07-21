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

#ifndef RTOPPACK_TOP_ABS_HPP
#define RTOPPACK_TOP_ABS_HPP

#include "RTOpPack_RTOpTHelpers.hpp"

namespace RTOpPack {

/** \brief Transformation operator that takes absolute values of elements: <tt>z0[i] = abs(v0[i]), i=0...n-1</tt>.
 */
template<class Scalar>
class TOpAbs : public RTOpT<Scalar> {
public:
  /** \brief . */
  TOpAbs() : RTOpT<Scalar>("TOpAbs") {}
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void apply_op(
    const int   num_vecs,       const ConstSubVectorView<Scalar>         sub_vecs[]
    ,const int  num_targ_vecs,  const SubVectorView<Scalar>  targ_sub_vecs[]
    ,ReductTarget *reduct_obj
    ) const
    {
      RTOP_APPLY_OP_1_1(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      for( Teuchos_Index i = 0; i < subDim; ++i, v0_val += v0_s, z0_val += z0_s ) {
        *z0_val = Teuchos::ScalarTraits<Scalar>::magnitude(*v0_val);
      }
    }
  //@}
}; // class TOpAbs

} // namespace RTOpPack

#endif // RTOPPACK_TOP_ABS_HPP
