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

#ifndef RTOPPACK_SUNDIALS_OPS_HPP
#define RTOPPACK_SUNDIALS_OPS_HPP

#include "RTOpPack_RTOpTHelpers.hpp"
#include "RTOpPack_RTOpT.hpp"
#include "RTOpPack_Types.hpp"


namespace RTOpPack 
{
  /* ------------------------------------------------------------------------------
   *
   * RTOp implementation of operations required by the SUNDIALS N_Vector interface
   * but not present in the Thyra StdVectorOps set. A few of these methods are
   * minor variations on the standard ops and could be composed from them, but 
   * composing them in that way requires vector copies not required by the operations
   * themselves. Thus for the sake of performance, we've implemented specialized RTOps
   * for all N_Vector functions not supported exactly as standard ops. 
   *
   * The following methods can be implemented directly with standard operations:
   * 
   * N_VLinearSum   -->  Thyra::linear_combination()
   * N_VConst       -->  Thyra::put_scalar()
   * N_VAbs         -->  Thyra::abs()
   * N_VInv         -->  Thyra::reciprocal() 
   * N_VDotProd     -->  Thyra::dot() 
   * N_VL1Norm      -->  Thyra::norm_1() 
   * N_VMaxNorm     -->  Thyra::norm_inf()
   * N_VMin         -->  Thyra::min() 
   *
   * The following, while similar to functions in the standard ops set, must 
   * be implemented here for optimal efficiency. 
   * 
   * N_VProd            - elementwise product
   * N_VDiv             - elementwise divide
   * N_VScale           - scale by constant
   * N_VAddConst        - add a constant to all elements
   *
   *
   * Weighted norms:
   *
   * The Thyra standard ops do include a weighted norm, but N_Vector uses a different
   * definition of the weighting factor. Thyra defines the weighted L2Norm as 
   * Sqrt[Sum[w_i x_i^2]],  but SUNDIALS wants Sqrt[Sum[w_i^2 x_i^2]]. 
   * 
   * N_VWL2Norm         - weighted 2-norm
   * N_VWrmsNorm        - weighted norm, normalized by vector length
   * N_VWrmsNormMask    - weighted rms norm, ignoring certain masked elements
   *
   * Notice that because RTOps work with chunks of data, 
   * the normalization by vector length *cannot* be performed inside
   * the RTOp and must be deferred until the Thyra wrapper functions. 
   *
   * 
   * Finally, there are several methods unlike those in the standard ops, which
   * are of course implemented anew:
   *
   * N_VConstrMask
   * N_VInvTest
   * N_VCompare
   * N_VMinQuotient
   * 
   * We do not implement the data access methods
   *
   * N_VGetArrayPointer
   * N_VSetArrayPointer
   *
   * which means that the SUNDIALS dense and banded solvers cannot be used. 
   *
   * \author Kevin Long
   * \date 10 December 2005
   *
   * ------------------------------------------------------------------------------*/



  /**
   * Performs the operation z[i] = x[i]*y[i] for i = 0, 1, ..., N-1
   */
  template<class Scalar>
  class SUNDIALS_VProd : public ROpScalarTransformationBase<Scalar> 
  {
  public:
    /** \brief . */
    SUNDIALS_VProd()
      : RTOpT<Scalar>("SUNDIALS_VProd"), ROpScalarTransformationBase<Scalar>() 
    {}
    /** @name Overridden from RTOpT */
    //@{
    /** \brief . */
    void apply_op(const int num_vecs, 
                  const ConstSubVectorView<Scalar> sub_vecs[],
                  const int  num_targ_vecs,  
                  const SubVectorView<Scalar> targ_sub_vecs[],
                  ReductTarget *reduct_obj) const
    {
      RTOP_APPLY_OP_2_1(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      for( index_type i = 0; i < subDim; 
           ++i,  v0_val += v0_s,  v1_val += v1_s, z0_val += z0_s ) 
        {
          *z0_val = (*v0_val) * (*v1_val);
        }
    }
    //@}
  }; 

  /**
   * Performs the operation z[i] = x[i]/y[i] for i = 0, 1, ..., N-1
   */
  template<class Scalar>
  class SUNDIALS_VDiv : public ROpScalarTransformationBase<Scalar> 
  {
  public:
    /** \brief . */
    SUNDIALS_VDiv()
      : RTOpT<Scalar>("SUNDIALS_VDiv"), ROpScalarTransformationBase<Scalar>() 
    {}
    /** @name Overridden from RTOpT */
    //@{
    /** \brief . */
    void apply_op(const int num_vecs, 
                  const ConstSubVectorView<Scalar> sub_vecs[],
                  const int  num_targ_vecs,  
                  const SubVectorView<Scalar> targ_sub_vecs[],
                  ReductTarget *reduct_obj) const
    {
      RTOP_APPLY_OP_2_1(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      for( index_type i = 0; i < subDim; 
           ++i,  v0_val += v0_s,  v1_val += v1_s, z0_val += z0_s ) 
        {
          *z0_val = (*v0_val) / (*v1_val);
        }
    }
    //@}
  }; 


  /**
   * Performs the operation z = c*x
   */
  template<class Scalar>
  class SUNDIALS_VScale : public ROpScalarTransformationBase<Scalar> 
  {
  public:
    /** \brief . */
    void alpha( const Scalar& alpha ) { this->scalarData(alpha); }
    /** \brief . */
    Scalar alpha() const { return this->scalarData(); }
    /** \brief . */
    SUNDIALS_VScale(const Scalar& alpha)
      : RTOpT<Scalar>("SUNDIALS_VScale"), ROpScalarTransformationBase<Scalar>(alpha) 
    {}
    /** @name Overridden from RTOpT */
    //@{
    /** \brief . */
    void apply_op(const int num_vecs, 
                  const ConstSubVectorView<Scalar> sub_vecs[],
                  const int  num_targ_vecs,  
                  const SubVectorView<Scalar> targ_sub_vecs[],
                  ReductTarget *reduct_obj) const
    {
      RTOP_APPLY_OP_1_1(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      Scalar a = alpha();
      for( index_type i = 0; i < subDim; 
           ++i,  v0_val += v0_s,  z0_val += z0_s ) 
        {
          *z0_val = a * (*v0_val);
        }
    }
    //@}
  }; 



  /**
   * Performs the operation z[i] = x[i] + b   for i = 0, 1, ..., N-1
   */
  template<class Scalar>
  class SUNDIALS_VAddConst : public ROpScalarTransformationBase<Scalar> 
  {
  public:
    /** \brief . */
    void alpha( const Scalar& alpha ) { this->scalarData(alpha); }
    /** \brief . */
    Scalar alpha() const { return this->scalarData(); }
    /** \brief . */
    SUNDIALS_VAddConst(const Scalar& alpha)
      : RTOpT<Scalar>("SUNDIALS_VAddConst"), ROpScalarTransformationBase<Scalar>(alpha) 
    {}
    /** @name Overridden from RTOpT */
    //@{
    /** \brief . */
    void apply_op(const int num_vecs, 
                  const ConstSubVectorView<Scalar> sub_vecs[],
                  const int  num_targ_vecs,  
                  const SubVectorView<Scalar> targ_sub_vecs[],
                  ReductTarget *reduct_obj) const
    {
      RTOP_APPLY_OP_1_1(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      Scalar a = alpha();
      for( index_type i = 0; i < subDim; 
           ++i,  v0_val += v0_s,  z0_val += z0_s ) 
        {
          *z0_val = (*v0_val) + a;
        }
    }
    //@}
  }; 




  /** 
   * Returns the weighted root mean square norm of x with weight
   * vector w:
   * \code
   *         sqrt [(sum (i = 0 to N-1) {(x[i]*w[i])^2})/N]
   * \endcode
   * Note that this is a different definition of the weighted norm than
   * used in ROp_WeightedNorm2(). 
   *
   * \author K. Long
   */
  
  template <class Scalar>
  class SUNDIALS_VWL2Norm : public ROpScalarReductionBase<Scalar>
  {
  public:
    /** \brief . */
    SUNDIALS_VWL2Norm()
      : RTOpT<Scalar>(""),
        ROpScalarReductionBase<Scalar>(0.0)
    {;}


    
    Scalar operator()(const ReductTarget& reduct_obj ) const 
    { return this->getRawVal(reduct_obj); }
    /** @name Overridden from RTOpT */
    //@{
      
    /** \brief . */
    void apply_op(
                  const int num_vecs,    
                  const ConstSubVectorView<Scalar> sub_vecs[],
                  const int num_targ_vecs, 
                  const SubVectorView<Scalar> targ_sub_vecs[],
                  ReductTarget *reduct_obj) const
    {
      ReductTargetScalar<Scalar> &_reduct_obj 
        = Teuchos::dyn_cast<ReductTargetScalar<Scalar> >(*reduct_obj); 
      RTOP_APPLY_OP_2_0(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      Scalar sum = this->getRawVal(*reduct_obj);
      for( index_type i = 0; i < subDim; ++i,  
             v0_val += v0_s,  
             v1_val += v1_s) 
        {
          double x = (*v0_val) * (*v1_val);
          sum += x*x;
        }
      _reduct_obj.set(sum);
      
    }
    //@}
  }; 


  /** 
   * Returns the weighted root mean square norm of x with weight
   * vector w:
   * \code
   *         sqrt [(sum (i = 0 to N-1) {(x[i]*w[i])^2})/N]
   * \endcode
   * Note that this is a *different* definition of the weighted norm than
   * used in ROp_WeightedNorm2(). 
   *
   * \author K. Long
   */
  
  template <class Scalar>
  class SUNDIALS_VWrmsMaskNorm : public ROpScalarReductionBase<Scalar>
  {
  public:
    
    /** \brief . */
    SUNDIALS_VWrmsMaskNorm()
      : RTOpT<Scalar>(""),
        ROpScalarReductionBase<Scalar>(0.0)
    {;}


    
    Scalar operator()(const ReductTarget& reduct_obj ) const 
    { return this->getRawVal(reduct_obj); }
    /** @name Overridden from RTOpT */
    //@{
      
    /** \brief . */
    void apply_op(
                  const int num_vecs,    
                  const ConstSubVectorView<Scalar> sub_vecs[],
                  const int num_targ_vecs, 
                  const SubVectorView<Scalar> targ_sub_vecs[],
                  ReductTarget *reduct_obj) const
    {
      ReductTargetScalar<Scalar> &_reduct_obj 
        = Teuchos::dyn_cast<ReductTargetScalar<Scalar> >(*reduct_obj); 
      RTOP_APPLY_OP_3_0(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      Scalar sum = this->getRawVal(*reduct_obj);
      const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
      for( index_type i = 0; i < subDim; ++i,  
             v0_val += v0_s,  
             v1_val += v1_s,  
             v2_val += v2_s) 
        {
          if (*v2_val > zero)
            {
              double x = (*v0_val) * (*v1_val);
              sum += x*x;
            }
        }
      _reduct_obj.set(sum);
      
    }
    //@}
  }; 


  /** 
   *   Performs the operation :
   * \code
   *       m[i] = 1.0 if constraint test fails for x[i]
   *       m[i] = 0.0 if constraint test passes for x[i]
   * \endcode
   *   where the constraint tests are as follows:
   * \code
   *      If c[i] = +2.0, then x[i] must be >  0.0.
   *      If c[i] = +1.0, then x[i] must be >= 0.0.
   *      If c[i] = -1.0, then x[i] must be <= 0.0.
   *      If c[i] = -2.0, then x[i] must be <  0.0.
   * \endcode
   *   This routine returns a boolean FALSE if any element failed
   *   the constraint test, TRUE if all passed. It also sets a
   *   mask vector m, with elements equal to 1.0 where the
   *   corresponding constraint test failed, and equal to 0.0
   *   where the constraint test passed.
   *   This routine is specialized in that it is used only for
   *   constraint checking.
   *
   * \author K. Long
   */
  
  template<class Scalar>
  class SUNDIALS_VConstrMask : public RTOpBoolReduceAndTransform<Scalar>
  {
  public:
    
    /** \brief . */
    SUNDIALS_VConstrMask()
      : RTOpT<Scalar>(""), 
        RTOpBoolReduceAndTransform<Scalar>() 
    {}
    
    /** @name Overridden from RTOpT */
    //@{
    /** \brief . */
    void apply_op(
                  const int num_vecs,    
                  const ConstSubVectorView<Scalar> sub_vecs[],
                  const int num_targ_vecs, 
                  const SubVectorView<Scalar> targ_sub_vecs[],
                  ReductTarget *reduct_obj) const
    {
      ReductTargetScalar<index_type>& _reduct_obj 
        = Teuchos::dyn_cast<ReductTargetScalar<index_type> >(*reduct_obj); 
      RTOP_APPLY_OP_2_1(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
      bool feasible = _reduct_obj.get();
      for( index_type i = 0; i < subDim; ++i,  
             v0_val += v0_s,  
             v1_val += v1_s,  
             z0_val += z0_s ) 
        {
          Scalar x_i = *v0_val;
          Scalar c_i = *v1_val;
          bool ok = true;
          if (c_i > 1.5 && c_i < 2.5)
            {
              ok = x_i > zero;
            }
          else if (c_i > 0.5 && c_i < 1.5)
            {
              ok = x_i >= zero;
            }
          else if (c_i < -0.5 && c_i > -1.5)
            {
              ok = x_i <= zero;
            }
          else if (c_i < -1.5)
            {
              ok = x_i < zero;
            }
          else
            {
              TEST_FOR_EXCEPTION(true, runtime_error,
                                 "illegal constraint flag = " << c_i
                                 << ". Allowed values are {-2.0, -1.0, 1.0, 2.0}");
            }

          /* Set the result element to 1 if the point is infeasible. */
          if (ok)
            {
              *z0_val = 0.0;
            }
          else
            {
              feasible = false;
              *z0_val = 1.0;
            }
        }
      _reduct_obj.set(feasible);
    }
    //@}
  }; 



  /** 
   * Performs the operation :
   * \code
   *       minq  = min ( num[i]/denom[i]) over all i such that
   *       denom[i] != 0.
   *   This routine returns the minimum of the quotients obtained
   *   by term-wise dividing num[i] by denom[i]. A zero element
   *   in denom will be skipped. If no such quotients are found,
   *   then the large value BIG_REAL is returned.
   * \endcode
   *
   * \author K. Long
   */
  
  template <class Scalar>
  class SUNDIALS_VMinQuotient : public ROpScalarReductionBase<Scalar>
  {
  public:
    
    /** \brief . */
    SUNDIALS_VMinQuotient()
      : RTOpT<Scalar>(""),
        ROpScalarReductionBase<Scalar>(Teuchos::ScalarTraits<Scalar>::rmax())
    {;}


    
    Scalar operator()(const ReductTarget& reduct_obj ) const 
    { return this->getRawVal(reduct_obj); }
    /** @name Overridden from RTOpT */
    //@{
    /** \brief . */
    void reduce_reduct_objs(const ReductTarget& in_reduct_obj, 
                            ReductTarget* inout_reduct_obj) const
    {
      const Scalar in_min_ele    = this->getRawVal(in_reduct_obj);
      const Scalar inout_min_ele = this->getRawVal(*inout_reduct_obj);
      this->setRawVal( in_min_ele < inout_min_ele ? in_min_ele : inout_min_ele, inout_reduct_obj );
    }
      
    /** \brief . */
    void apply_op(
                  const int num_vecs,    
                  const ConstSubVectorView<Scalar> sub_vecs[],
                  const int num_targ_vecs, 
                  const SubVectorView<Scalar> targ_sub_vecs[],
                  ReductTarget *reduct_obj) const
    {
      ReductTargetScalar<Scalar> &_reduct_obj 
        = Teuchos::dyn_cast<ReductTargetScalar<Scalar> >(*reduct_obj); 
      RTOP_APPLY_OP_2_0(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      Scalar min_ele = this->getRawVal(*reduct_obj);
      const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
      for( index_type i = 0; i < subDim; ++i,  
             v0_val += v0_s,  
             v1_val += v1_s) 
        {
          if (*v1_val != zero)
            {
              Scalar x = (*v0_val)/(*v1_val);
              if (x < min_ele) min_ele = x;
            }
        }
      _reduct_obj.set(min_ele);
      
    }
    //@}
  }; 




  /** 
   *
   * Performs the operation z[i] = 1/x[i] with a test for
   * x[i] == 0.0 before inverting x[i]. If z[i] == 0.0,
   * skip the operation. 
   *
   * \author K. Long
   */
  
  template<class Scalar>
  class SUNDIALS_VInvTest : public RTOpBoolReduceAndTransform<Scalar>
  {
  public:
    
    /** \brief . */
    SUNDIALS_VInvTest()
      : RTOpT<Scalar>(""), 
        RTOpBoolReduceAndTransform<Scalar>() 
    {}
    
    /** @name Overridden from RTOpT */
    //@{
    /** \brief . */
    void apply_op(
                  const int num_vecs,    
                  const ConstSubVectorView<Scalar> sub_vecs[],
                  const int num_targ_vecs, 
                  const SubVectorView<Scalar> targ_sub_vecs[],
                  ReductTarget *reduct_obj) const
    {
      ReductTargetScalar<index_type>& _reduct_obj 
        = Teuchos::dyn_cast<ReductTargetScalar<index_type> >(*reduct_obj); 
      RTOP_APPLY_OP_1_1(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
      const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
      bool pass =  _reduct_obj.get();
      for( index_type i = 0; i < subDim; ++i,  
             v0_val += v0_s,  
             z0_val += z0_s ) 
        {
          if (*v0_val != zero) 
            {
              *z0_val = one/(*v0_val);
            }
          else 
            {
              *z0_val = 0.0;
              pass = false;
            }
        }
      _reduct_obj.set(pass);
    }
    //@}
  }; 



  /** 
   *  Performs the operation
   *     z[i] = 1.0 if ABS(x[i]) >= alpha
   *            0.0 otherwise
   *
   * \author K. Long
   */
  
  template<class Scalar>
  class SUNDIALS_VCompare : public ROpScalarTransformationBase<Scalar> {
  public:
    /** \brief . */
    void alpha( const Scalar& alpha ) { this->scalarData(alpha); }
    /** \brief . */
    Scalar alpha() const { return this->scalarData(); }
    /** \brief . */
    SUNDIALS_VCompare( const Scalar &alpha)
      : RTOpT<Scalar>("SUNDIALS_Compare"), 
        ROpScalarTransformationBase<Scalar>(alpha) 
    {}
    /** @name Overridden from RTOpT */
    //@{
    /** \brief . */
    void apply_op(
                  const int num_vecs,    
                  const ConstSubVectorView<Scalar> sub_vecs[],
                  const int num_targ_vecs, 
                  const SubVectorView<Scalar> targ_sub_vecs[],
                  ReductTarget *reduct_obj) const
    {
      RTOP_APPLY_OP_1_1(num_vecs,sub_vecs,num_targ_vecs,targ_sub_vecs);
      const Scalar one = Teuchos::ScalarTraits<Scalar>::one();
      const Scalar zero = Teuchos::ScalarTraits<Scalar>::zero();
      for( index_type i = 0; i < subDim; ++i,  v0_val += v0_s,  z0_val += z0_s ) 
        {
          if (fabs(*v0_val) >= alpha()) *z0_val = one;
          else *z0_val = zero;
        }
    }
    //@}
  };

}


#endif 
