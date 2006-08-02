// @HEADER
// ***********************************************************************
// 
//      Thyra: Interfaces and Support Code for the Interoperability 
//                  of Abstract Numerical Algorithms 
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

#ifndef THYRA_SUNDIALS_OPS_HPP
#define THYRA_SUNDIALS_OPS_HPP



#include "Thyra_VectorBase.hpp"
#include "RTOpPack_SUNDIALS_Ops.hpp"

namespace Thyra
{ 
  /** \brief element-wise product */
  template<class Scalar> inline
  void VProd(const VectorBase<Scalar>& x, 
             const VectorBase<Scalar>& y,
             VectorBase<Scalar>* z)
  {
    RTOpPack::SUNDIALS_VProd<Scalar> op;
    const VectorBase<Scalar>* vecs[]      = { &x, &y };
    VectorBase<Scalar>*       targ_vecs[] = { z };
    applyOp<Scalar>(op,2,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
  }  



  
  /** \brief element-wise divide */
  template<class Scalar> inline
  void VDiv(const VectorBase<Scalar>& x, 
            const VectorBase<Scalar>& y,
            VectorBase<Scalar>* z)
  {
    RTOpPack::SUNDIALS_VDiv<Scalar> op;
    const VectorBase<Scalar>* vecs[]      = { &x, &y };
    VectorBase<Scalar>*       targ_vecs[] = { z };
    applyOp<Scalar>(op,2,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
  }



  /** \brief scale by a constant */
  template<class Scalar> inline
  void VScale(const Scalar& c,
              const VectorBase<Scalar>& x, 
              VectorBase<Scalar>* z)
  {
    RTOpPack::SUNDIALS_VScale<Scalar> op(c);
    const VectorBase<Scalar>* vecs[]      = { &x };
    VectorBase<Scalar>*       targ_vecs[] = { z };
    applyOp<Scalar>(op,1,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
  }



  /** \brief add a constant */
  template<class Scalar> inline
  void VAddConst(const Scalar& c,
                 const VectorBase<Scalar>& x, 
                 VectorBase<Scalar>* z)
  {
    RTOpPack::SUNDIALS_VAddConst<Scalar> op(c);
    const VectorBase<Scalar>* vecs[]      = { &x };
    VectorBase<Scalar>*       targ_vecs[] = { z };
    applyOp<Scalar>(op,1,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
  }

  
  /** \brief weighted L2 norm */
  template<class Scalar> inline
  Scalar VWL2Norm(const VectorBase<Scalar>& x, 
                  const VectorBase<Scalar>& w)
  {
    RTOpPack::SUNDIALS_VWL2Norm<Scalar> op;
    Teuchos::RefCountPtr<RTOpPack::ReductTarget> red_targ 
      = op.reduct_obj_create();

    const VectorBase<Scalar>* vecs[]      = { &x, &w };
    
    applyOp<Scalar>(op,2,vecs,0,(VectorBase<Scalar>**)NULL,&*red_targ);
    return sqrt(op(*red_targ));
  }


  /** \brief weighted rms norm */
  template<class Scalar> inline
  Scalar VWrmsNorm(const VectorBase<Scalar>& x, 
                   const VectorBase<Scalar>& w)
  {
    return VWL2Norm(x,w)/sqrt(x.space()->dim());
  }



  /** \brief weighted L2 norm with mask */
  template<class Scalar> inline
  Scalar VWrmsMaskNorm(const VectorBase<Scalar>& x, 
                       const VectorBase<Scalar>& w, 
                       const VectorBase<Scalar>& id)
  {
    RTOpPack::SUNDIALS_VWrmsMaskNorm<Scalar> op;
    Teuchos::RefCountPtr<RTOpPack::ReductTarget> red_targ 
      = op.reduct_obj_create();

    const VectorBase<Scalar>* vecs[]      = { &x, &w, &id };
    
    applyOp<Scalar>(op,3,vecs,0,(VectorBase<Scalar>**)NULL,&*red_targ);
    return sqrt(op(*red_targ));
  }
  
  /** \brief minimum quotient */
  template<class Scalar> inline
  Scalar VMinQuotient(const VectorBase<Scalar>& num, 
                      const VectorBase<Scalar>& denom)
  {
    RTOpPack::SUNDIALS_VMinQuotient<Scalar> op;
    Teuchos::RefCountPtr<RTOpPack::ReductTarget> red_targ 
      = op.reduct_obj_create();

    const VectorBase<Scalar>* vecs[]      = { &num, &denom };
    
    applyOp<Scalar>(op,2,vecs,0,(VectorBase<Scalar>**)NULL,&*red_targ);
    return op(*red_targ);
  }


  /** \brief constraint mask */
  template<class Scalar> inline
  bool VConstrMask(const VectorBase<Scalar>& x, 
                   const VectorBase<Scalar>& constraint, 
                   VectorBase<Scalar>* mask)
  {
    RTOpPack::SUNDIALS_VConstrMask<Scalar> op;
    Teuchos::RefCountPtr<RTOpPack::ReductTarget> red_targ 
      = op.reduct_obj_create();

    const VectorBase<Scalar>* vecs[]      = { &x, &constraint };
    VectorBase<Scalar>*       targ_vecs[] = { mask  };
    
    applyOp<Scalar>(op,2,vecs,1,targ_vecs,&*red_targ);
    return op(*red_targ);
  }

  /** \brief invert with nonzero test */
  template<class Scalar> inline
  bool VInvTest(const VectorBase<Scalar>& v_rhs, 
                VectorBase<Scalar>* v_lhs)
  {
    RTOpPack::SUNDIALS_VInvTest<Scalar> op;
    Teuchos::RefCountPtr<RTOpPack::ReductTarget> red_targ 
      = op.reduct_obj_create();

    const VectorBase<Scalar>* vecs[]      = { &v_rhs };
    VectorBase<Scalar>*       targ_vecs[] = { v_lhs  };
    
    applyOp<Scalar>(op,1,vecs,1,targ_vecs,&*red_targ);
    return op(*red_targ);
  }
  
  /** \brief compare to a scalar */
  template<class Scalar> inline
  void VCompare(const Scalar& alpha,
                const VectorBase<Scalar>& x,
                VectorBase<Scalar>* y)
  {
    TEST_FOR_EXCEPTION(y==NULL,std::logic_error,
                       "compareToScalar(...), Error!");
    RTOpPack::SUNDIALS_VCompare<Scalar> op(alpha);
    const VectorBase<Scalar>* vecs[]      = { &x };
    VectorBase<Scalar>* targ_vecs[] = {y };
    applyOp<Scalar>(op,1,vecs,1,targ_vecs,(RTOpPack::ReductTarget*)NULL);
  }
} 

#endif 
