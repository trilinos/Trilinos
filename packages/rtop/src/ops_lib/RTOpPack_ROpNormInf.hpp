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

#ifndef RTOPPACK_ROP_NORMINF_HPP
#define RTOPPACK_ROP_NORMINF_HPP

#include "RTOpPack_RTOpTHelpers.hpp"

namespace RTOpPack {

/** \brief Infinity norm reduction operator: <tt>result = sum( |v0[i]|, i=0...n-1 )</tt>.
 */
template<class Scalar>
class ROpNormInf
  : public ROpScalarReductionBase<Scalar,typename Teuchos::ScalarTraits<Scalar>::magnitudeType>
{
public:
  /** \brief . */
  typedef   typename Teuchos::ScalarTraits<Scalar>::magnitudeType  ScalarMag;
  /** \brief . */
  ROpNormInf() : RTOpT<Scalar>("ROpNormInf") {}
  /** \brief . */
  ScalarMag operator()(const ReductTarget& reduct_obj ) const
    { return this->getRawVal(reduct_obj); }
  /** @name Overridden from RTOpT */
  //@{
  /** \brief . */
  void reduce_reduct_objs(
    const ReductTarget& _in_reduct_obj, ReductTarget* _inout_reduct_obj
    ) const
    {
      using Teuchos::dyn_cast;
      const ReductTargetScalar<ScalarMag>
        &in_reduct_obj = dyn_cast<const ReductTargetScalar<ScalarMag> >(_in_reduct_obj);
      ReductTargetScalar<ScalarMag>
        &inout_reduct_obj = dyn_cast<ReductTargetScalar<ScalarMag> >(*_inout_reduct_obj);
      if( in_reduct_obj.get() > inout_reduct_obj.get() ) inout_reduct_obj.set(in_reduct_obj.get());
    }
  /** \brief . */
  void apply_op(
    const int   num_vecs,       const ConstSubVectorView<Scalar>    sub_vecs[]
    ,const int  num_targ_vecs,  const SubVectorView<Scalar>         targ_sub_vecs[]
    ,ReductTarget *_reduct_obj
    ) const
    {
      using Teuchos::dyn_cast;
      ReductTargetScalar<ScalarMag>
        &reduct_obj = dyn_cast<ReductTargetScalar<ScalarMag> >(*_reduct_obj); 
      RTOP_APPLY_OP_1_0(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      ScalarMag norm_inf = reduct_obj.get();
      if( v0_s == 1 ) {
        for( Teuchos_Index i = 0; i < subDim; ++i ) {
          const ScalarMag
            mag = Teuchos::ScalarTraits<Scalar>::magnitude(*v0_val++);
          norm_inf = mag > norm_inf ? mag : norm_inf;
        }
      }
      else {
        for( Teuchos_Index i = 0; i < subDim; ++i, v0_val += v0_s ) {
          const ScalarMag
            mag = Teuchos::ScalarTraits<Scalar>::magnitude(*v0_val);
          norm_inf = mag > norm_inf ? mag : norm_inf;
        }
      }
      reduct_obj.set(norm_inf);
    }
  //@}
}; // class ROpNormInf

} // namespace RTOpPack

#endif // RTOPPACK_ROP_NORMINF_HPP
