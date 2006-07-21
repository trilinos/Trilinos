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

#ifndef RTOPPACK_ROP_DOT_PROD_HPP
#define RTOPPACK_ROP_DOT_PROD_HPP

#include "RTOpPack_RTOpTHelpers.hpp"

namespace RTOpPack {

/** \brief Dot product reduction operator: <tt>result = sum( conj(v0[i])*v1[i], i=0...n-1 )</tt>.
 */
template<class Scalar>
class ROpDotProd : public ROpScalarReductionBase<Scalar> {
public:
  /** \brief . */
  ROpDotProd() : RTOpT<Scalar>("ROpDotProd") {}
  /** \brief . */
  Scalar operator()(const ReductTarget& reduct_obj ) const { return this->getRawVal(reduct_obj); }
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void apply_op(
    const int   num_vecs,       const ConstSubVectorView<Scalar>         sub_vecs[]
    ,const int  num_targ_vecs,  const SubVectorView<Scalar>  targ_sub_vecs[]
    ,ReductTarget *_reduct_obj
    ) const
    {
      using Teuchos::dyn_cast;
      typedef Teuchos::ScalarTraits<Scalar> ST;
      ReductTargetScalar<Scalar> &reduct_obj = dyn_cast<ReductTargetScalar<Scalar> >(*_reduct_obj); 
      RTOP_APPLY_OP_2_0(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      Scalar dot_prod = reduct_obj.get();
      if( v0_s == 1 && v1_s == 1 ) {
        for( Teuchos_Index i = 0; i < subDim; ++i )
          dot_prod += ST::conjugate(*v0_val++)*(*v1_val++);
      }
      else {
        for( Teuchos_Index i = 0; i < subDim; ++i, v0_val += v0_s, v1_val += v1_s )
          dot_prod += ST::conjugate(*v0_val)*(*v1_val);
      }
      reduct_obj.set(dot_prod);
    }
  //@}
}; // class ROpDotProd

} // namespace RTOpPack

#endif // RTOPPACK_ROP_DOT_PROD_HPP
