// @HEADER
// ***********************************************************************
// 
//    Thyra: Interfaces and Support for Abstract Numerical Algorithms
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

#ifndef THYRA_LINEARCOMBINATION_IMPL_HPP
#define THYRA_LINEARCOMBINATION_IMPL_HPP


#include "Thyra_ConfigDefs.hpp"
#include "Thyra_LinearCombinationDecl.hpp"
#include "Thyra_LinearOperatorImpl.hpp"
#include "Thyra_VectorHandleOpsImpl.hpp"

namespace Thyra
{
  /* -------- methods of OpTimesLC ------ */

  template <class Scalar, class Node> inline
  OpTimesLC<Scalar, Node>::OpTimesLC(const Scalar& alpha, 
                               const Node& x)
    : alpha_(alpha), op_(), x_(x) 
  {;}

  template <class Scalar, class Node> inline
  OpTimesLC<Scalar, Node>
  ::OpTimesLC(const Scalar& alpha,
              const Thyra::ConstLinearOperator<Scalar>& op, 
              const Node& x)
    : alpha_(alpha), op_(op), x_(x) 
  {;}
  
  
  template <class Scalar, class Node> inline
  void OpTimesLC<Scalar, Node>::evalInto(Thyra::Vector<Scalar>& result) const
  {
    if (op_.constPtr().get() != 0)
      {
        op_.apply(x_.evalToConst(), result, alpha_);
      }
    else
      {
        x_.evalInto(result);
        if (alpha_ != one()) scale(result, alpha_);
      }
  }

  template <class Scalar, class Node> inline
  void OpTimesLC<Scalar, Node>::addInto(Thyra::Vector<Scalar>& result,
                            LCSign sign) const
  {
    Scalar s = sign;
    if (op_.constPtr().get() != 0)
      {
        Thyra::Vector<Scalar> tmp;
        op_.apply(x_.evalToConst(), tmp);
        axpy(s*alpha_, tmp, result);
      }
    else
      {
        axpy(s*alpha_, x_.evalToConst(), result);
      }
  } 

  template <class Scalar, class Node> inline
  Thyra::Vector<Scalar> OpTimesLC<Scalar, Node>::formVector() const 
  {
    Thyra::Vector<Scalar> result;
    if (op_.constPtr().get() != 0)
      {
        result = op_.range().createMember();
        cout << "result = " << result.constPtr().get() << endl;
        op_.apply(x_.evalToConst(), result, alpha_);
        cout << "result = " << result.constPtr().get() << endl;
      }
    else
      {
        result = Thyra::formVector(x_);
        if (alpha_ != one()) scale(result, alpha_);    
      }
    cout << "result = " << result.constPtr().get() << endl;
    return result;
  }


  template <class Scalar, class Node> inline
  bool OpTimesLC<Scalar, Node>::containsVector(const Thyra::VectorBase<Scalar>* vec) const 
  {return x_.containsVector(vec);}


 


  
  /* ------------------------ methods of LC2 --------------------------- */
  
  template <class Scalar, class Node1, class Node2> inline
  LC2<Scalar, Node1, Node2>::LC2(const Node1& x1, const Node2& x2, LCSign sign)
    : x1_(x1), x2_(x2), sign_(sign) 
  {;}

  template <class Scalar, class Node1, class Node2> inline
  bool LC2<Scalar, Node1, Node2>::containsVector(const Thyra::VectorBase<Scalar>* vec) const
  {return x1_.containsVector(vec) || x2_.containsVector(vec);}

  template <class Scalar, class Node1, class Node2> inline
  void LC2<Scalar, Node1, Node2>::evalInto(Thyra::Vector<Scalar>& result) const
  {
    x1_.evalInto(result);
    x2_.addInto(result, sign_);
  } 

  template <class Scalar, class Node1, class Node2> inline
  void LC2<Scalar, Node1, Node2>::addInto(Thyra::Vector<Scalar>& result,
                                          Thyra::LCSign sign) const
  {
    x1_.addInto(result, sign);
    if (sign_*sign < 0) x2_.addInto(result, LCSubtract);
    else x2_.addInto(result, LCAdd);
  }

  template <class Scalar, class Node1, class Node2> inline
  Thyra::Vector<Scalar> LC2<Scalar, Node1, Node2>::formVector() const
  {
    Thyra::Vector<Scalar> result = Thyra::formVector(x1_);
    x2_.addInto(result, sign_);
    return result;
  }

}

namespace Thyra
{

  /* ------------------------ global methods ----------------------- */


  /*======================================================================
   *
   *    scalar times vector
   *
   *======================================================================*/

  /* scalar * vec */
  template <class Scalar> inline
  OpTimesLC<Scalar, Thyra::ConstVector<Scalar> > 
  operator*(const Scalar& alpha, 
            const Thyra::ConstVector<Scalar>& x)
  {
    return OpTimesLC<Scalar, Thyra::ConstVector<Scalar> >(alpha, x);
  } 

  /* vec * scalar */
  template <class Scalar> inline
  OpTimesLC<Scalar, Thyra::ConstVector<Scalar> > 
  operator*(const Thyra::ConstVector<Scalar>& x, 
            const Scalar& alpha)
  {
    return OpTimesLC<Scalar, Thyra::ConstVector<Scalar> >(alpha, x);
  }


  /*======================================================================
   *
   *    scalar times OpTimesLC
   *
   *======================================================================*/

  /* scalar * OpTimesLC */
  template <class Scalar, class Node> inline
  OpTimesLC<Scalar, Node> 
  operator*(const Scalar& alpha, 
            const OpTimesLC<Scalar, Node>& x)
  {
    return OpTimesLC<Scalar, Node>(alpha * x.alpha(), x.op(), x.node());
  }

  /* OpTimesLC * scalar */
  template <class Scalar, class Node> inline
  OpTimesLC<Scalar, Node> 
  operator*(const OpTimesLC<Scalar, Node>& x, const Scalar& alpha)
  {
    return alpha * x;
  }


  /*======================================================================
   *
   *    scalar times LC2
   *
   *======================================================================*/

  /* scalar * LC2 */
  template <class Scalar, class Node1, class Node2> inline
  OpTimesLC<Scalar, LC2<Scalar, Node1, Node2> > 
  operator*(const Scalar& alpha, 
            const LC2<Scalar, Node1, Node2>& x)
  {
    return OpTimesLC<Scalar, LC2<Scalar, Node1, Node2> >(alpha, x);
  }

  /* LC2 * scalar */
  template <class Scalar, class Node1, class Node2> inline
  OpTimesLC<Scalar, LC2<Scalar, Node1, Node2> > 
  operator*(const LC2<Scalar, Node1, Node2>& x, const Scalar& alpha)
  {
    return alpha * x;
  }
  


  /*======================================================================
   *
   *    operator times [vectors, OpTimesLC, LC2]
   *
   *======================================================================*/

  /* op * vec */
  template <class Scalar> inline
  OpTimesLC<Scalar, Thyra::ConstVector<Scalar> > 
  operator*(const ConstLinearOperator<Scalar>& op, 
            const Thyra::ConstVector<Scalar>& x)
  {
    return OpTimesLC<Scalar, Thyra::ConstVector<Scalar> >(Teuchos::ScalarTraits<Scalar>::one(), op, x);
  }


  /* op * OpTimesLC */
  template <class Scalar, class Node> inline
  OpTimesLC<Scalar, Node> 
  operator*(const ConstLinearOperator<Scalar>& op, 
            const OpTimesLC<Scalar, Node>& x)
  {
    TEST_FOR_EXCEPTION(op.constPtr().get()==0, runtime_error,
                       "null operator in LinearOperator * ( OpTimesLC )");
    if (x.op().constPtr().get()==0)
      {
        return OpTimesLC<Scalar, Node>(x.alpha(), op, x.node());
      }
    else
      {
        return OpTimesLC<Scalar, Node>(x.alpha(), op * x.op(), x.node());
      }
  }


  /* op * LC2 */
  template <class Scalar, class Node1, class Node2> inline
  OpTimesLC<Scalar, LC2<Scalar, Node1, Node2> > 
  operator*(const ConstLinearOperator<Scalar>& op, 
            const LC2<Scalar, Node1, Node2>& x)
  {
    return OpTimesLC<Scalar, LC2<Scalar, Node1, Node2> >(Teuchos::ScalarTraits<Scalar>::one(), op, x);
  }


  /*======================================================================
   *
   *    add/subtract vector, vector
   *
   *======================================================================*/
  
  /* vec + vec */
  template <class Scalar> inline
  LC2<Scalar, Thyra::ConstVector<Scalar>, Thyra::ConstVector<Scalar> >
  operator+(const Thyra::ConstVector<Scalar>& x1, 
            const Thyra::ConstVector<Scalar>& x2)
  {
    return LC2<Scalar, Thyra::ConstVector<Scalar>, Thyra::ConstVector<Scalar> >(x1, x2);
  }
  
  /* vec - vec */
  template <class Scalar> inline
  LC2<Scalar, Thyra::ConstVector<Scalar>, Thyra::ConstVector<Scalar> >
  operator-(const Thyra::ConstVector<Scalar>& x1, 
            const Thyra::ConstVector<Scalar>& x2)
  {
    return LC2<Scalar, Thyra::ConstVector<Scalar>, Thyra::ConstVector<Scalar> >(x1, x2, LCSubtract);
  }

  /*======================================================================
   *
   *    add/subtract vector, OpTimesLC
   *
   *======================================================================*/

  /* vec + OpTimesLC */
  template <class Scalar, class Node> inline
  LC2<Scalar, Thyra::ConstVector<Scalar>, OpTimesLC<Scalar, Node> >
  operator+(const Thyra::ConstVector<Scalar>& x1, 
            const OpTimesLC<Scalar, Node>& x2)
  {
    return LC2<Scalar, Thyra::ConstVector<Scalar>, OpTimesLC<Scalar, Node> >(x1, x2);
  }

  /* vec - OpTimesLC */
  template <class Scalar, class Node> inline
  LC2<Scalar, Thyra::ConstVector<Scalar>, OpTimesLC<Scalar, Node> >
  operator-(const Thyra::ConstVector<Scalar>& x1, 
            const OpTimesLC<Scalar, Node>& x2)
  {
    return LC2<Scalar, Thyra::ConstVector<Scalar>, OpTimesLC<Scalar, Node> >(x1, x2, 
                                                                 LCSubtract);
  }
  
  /* OpTimesLC + vec */
  template <class Scalar, class Node> inline
  LC2<Scalar, OpTimesLC<Scalar, Node>, Thyra::ConstVector<Scalar> >
  operator+(const OpTimesLC<Scalar, Node>& x1, 
            const Thyra::ConstVector<Scalar>& x2)
  {
    return LC2<Scalar, OpTimesLC<Scalar, Node>, Thyra::ConstVector<Scalar> >(x1, x2);
  }
  
  /* OpTimesLC - vec */
  template <class Scalar, class Node> inline
  LC2<Scalar, OpTimesLC<Scalar, Node>, Thyra::Vector<Scalar> >
  operator-(const OpTimesLC<Scalar, Node>& x1, 
            const Thyra::ConstVector<Scalar>& x2)
  {
    return LC2<Scalar, OpTimesLC<Scalar, Node>, Thyra::Vector<Scalar> >(x1, x2,
                                                                 LCSubtract);
  }

  
  /*======================================================================
   *
   *    add/subtract OpTimesLC, OpTimesLC
   *
   *======================================================================*/
  
  /* OpTimesLC + OpTimesLC */
  template <class Scalar, class Node1, class Node2> inline
  LC2<Scalar, OpTimesLC<Scalar, Node1>, OpTimesLC<Scalar, Node2> >
  operator+(const OpTimesLC<Scalar, Node1>& x1, 
            const OpTimesLC<Scalar, Node2>& x2)
  {
    return LC2<Scalar, OpTimesLC<Scalar, Node1>, 
      OpTimesLC<Scalar, Node2> >(x1, x2);
  }
  
  /* OpTimesLC - OpTimesLC */
  template <class Scalar, class Node1, class Node2> inline
  LC2<Scalar, OpTimesLC<Scalar, Node1>, OpTimesLC<Scalar, Node2> >
  operator-(const OpTimesLC<Scalar, Node1>& x1, 
            const OpTimesLC<Scalar, Node2>& x2)
  {
    return LC2<Scalar, OpTimesLC<Scalar, Node1>, 
      OpTimesLC<Scalar, Node2> >(x1, x2, LCSubtract);
  }
  

  
  /*======================================================================
   *
   *    add/subtract Vector, LC2
   *
   *======================================================================*/

  
  /* vec + LC2 */
  template <class Scalar, class Node1, class Node2> inline
  LC2<Scalar, Thyra::ConstVector<Scalar>, LC2<Scalar, Node1, Node2> >
  operator+(const Thyra::ConstVector<Scalar>& x1, 
            const LC2<Scalar, Node1, Node2>& x2)
  {
    return LC2<Scalar, Thyra::ConstVector<Scalar>, LC2<Scalar, Node1, Node2> >(x1, x2);
  }

  /* vec - LC2 */
  template <class Scalar, class Node1, class Node2> inline
  LC2<Scalar, Thyra::ConstVector<Scalar>, LC2<Scalar, Node1, Node2> >
  operator-(const Thyra::ConstVector<Scalar>& x1, 
            const LC2<Scalar, Node1, Node2>& x2)
  {
    return LC2<Scalar, Thyra::ConstVector<Scalar>, LC2<Scalar, Node1, Node2> >(x1, x2,
                                                                   LCSubtract);
  }


  /* LC2 + vec */
  template <class Scalar, class Node1, class Node2> inline
  LC2<Scalar, LC2<Scalar, Node1, Node2>, Thyra::ConstVector<Scalar> >
  operator+(const LC2<Scalar, Node1, Node2>& x1, 
            const Thyra::ConstVector<Scalar>& x2)
  {
    return LC2<Scalar, LC2<Scalar, Node1, Node2>, Thyra::ConstVector<Scalar> >(x1, x2);
  }

  /* LC2 - vec */
  template <class Scalar, class Node1, class Node2> inline
  LC2<Scalar, LC2<Scalar, Node1, Node2>, Thyra::ConstVector<Scalar> >
  operator-(const LC2<Scalar, Node1, Node2>& x1, 
            const Thyra::ConstVector<Scalar>& x2)
  {
    return LC2<Scalar, LC2<Scalar, Node1, Node2>, Thyra::ConstVector<Scalar> >(x1, x2,
                                                                   LCSubtract);
  }


  /*======================================================================
   *
   *    add/subtract OpTimesLC, LC2
   *
   *======================================================================*/


  /* OpTimesLC + LC2 */
  template <class Scalar, class Node0, class Node1, class Node2> inline
  LC2<Scalar, OpTimesLC<Scalar, Node0>, LC2<Scalar, Node1, Node2> > 
  operator+(const OpTimesLC<Scalar, Node0>& x1, 
            const LC2<Scalar, Node1, Node2>& x2)
  {
    return LC2<Scalar, OpTimesLC<Scalar, Node0>,
      LC2<Scalar, Node1, Node2> >(x1, x2);
  }

  /* OpTimesLC - LC2 */
  template <class Scalar, class Node0, class Node1, class Node2> inline
  LC2<Scalar, OpTimesLC<Scalar, Node0>, LC2<Scalar, Node1, Node2> > 
  operator-(const OpTimesLC<Scalar, Node0>& x1, 
            const LC2<Scalar, Node1, Node2>& x2)
  {
    return LC2<Scalar, OpTimesLC<Scalar, Node0>,
      LC2<Scalar, Node1, Node2> >(x1, x2, LCSubtract);
  }


  /* LC2 + OpTimesLC */
  template <class Scalar, class Node1, class Node2, class Node3> inline
  LC2<Scalar, LC2<Scalar, Node1, Node2>, OpTimesLC<Scalar, Node3> > 
  operator+(const LC2<Scalar, Node1, Node2>& x1, 
            const OpTimesLC<Scalar, Node3>& x2)
  {
    return LC2<Scalar, LC2<Scalar, Node1, Node2>, 
      OpTimesLC<Scalar, Node3> >(x1, x2);
  }

  /* LC2 - OpTimesLC */
  template <class Scalar, class Node1, class Node2, class Node3> inline
  LC2<Scalar, LC2<Scalar, Node1, Node2>, OpTimesLC<Scalar, Node3> > 
  operator-(const LC2<Scalar, Node1, Node2>& x1, 
            const OpTimesLC<Scalar, Node3>& x2)
  {
    return LC2<Scalar, LC2<Scalar, Node1, Node2>, 
      OpTimesLC<Scalar, Node3> >(x1, x2, LCSubtract);
  }


  /*======================================================================
   *
   *    add/subtract LC2, LC2
   *
   *======================================================================*/
  
  /* LC2 + LC2 */
  template <class Scalar, class Node1, class Node2, 
            class Node3, class Node4> inline
  LC2<Scalar, LC2<Scalar, Node1, Node2>, LC2<Scalar, Node3, Node4> >
  operator+(const LC2<Scalar, Node1, Node2>& x1, 
            const LC2<Scalar, Node3, Node4>& x2)
  {
    return LC2<Scalar, LC2<Scalar, Node1, Node2>, 
      LC2<Scalar, Node3, Node4> >(x1, x2);
  }

  /* LC2 - LC2 */
  template <class Scalar, class Node1, class Node2, 
            class Node3, class Node4> inline
  LC2<Scalar, LC2<Scalar, Node1, Node2>, LC2<Scalar, Node3, Node4> >
  operator-(const LC2<Scalar, Node1, Node2>& x1, 
            const LC2<Scalar, Node3, Node4>& x2)
  {
    return LC2<Scalar, LC2<Scalar, Node1, Node2>, 
      LC2<Scalar, Node3, Node4> >(x1, x2, LCSubtract);
  }


  /*======================================================================
   *
   *    assignment of [OpTimesLC, LC2] to vector
   *
   *======================================================================*/
  
  
  /* definition of assignment from 1-term linear combination to a vector */
  template <class Scalar> 
  template <class Node> inline
  Thyra::Vector<Scalar>& Thyra::Vector<Scalar>::operator=(const Thyra::OpTimesLC<Scalar, Node>& x)
  {
    if (this->ptr().get()==0)
      {
        *this = x.formVector();
      }
    else if (x.containsVector(this->ptr().get()))
      {
        Thyra::Vector<Scalar> rtn = x.formVector();
        acceptCopyOf(rtn);
      }
    else
      {
        x.evalInto(*this);
      }
    return *this;
  }

 
  /* definition of assignment from N-term linear combination to a vector */
  template <class Scalar>
  template <class Node1, class Node2> inline
  Thyra::Vector<Scalar>& Thyra::Vector<Scalar>::operator=(const Thyra::LC2<Scalar, Node1, Node2>& x)
  {
    if (this->ptr().get()==0)
      {
        *this = x.formVector();
      }
    else if (x.containsVector(this->ptr().get()))
      {
        Thyra::Vector<Scalar> rtn = x.formVector();
        acceptCopyOf(rtn);
      }
    else
      {
        x.evalInto(*this);
      }
    return *this;
  }

 

  /*======================================================================
   *
   *    construction of vectors from [OpTimesLC, LC2]
   *
   *======================================================================*/
   

  template <class Scalar>
  Thyra::ConstVector<Scalar>::ConstVector(const Thyra::ConvertibleToVector<Scalar>& x)
    : Teuchos::ConstHandle<Thyra::VectorBase<Scalar> >(x.formVector().constPtr())
  {;}
  /*
  template <class Scalar> 
  template <class Node> inline
  Thyra::ConstVector<Scalar>::ConstVector(const Thyra::OpTimesLC<Scalar, Node>& x)
    : Teuchos::ConstHandle<Thyra::VectorBase<Scalar> >(x.formVector().constPtr())
  {;}
  */

  template <class Scalar>
  template <class Node1, class Node2> inline
  Thyra::Vector<Scalar>::Vector(const Thyra::LC2<Scalar, Node1, Node2>& x)
    : Thyra::ConstVector<Scalar>(x.formVector().constPtr())
  {;}

  template <class Scalar> 
  template <class Node> inline
  Thyra::Vector<Scalar>::Vector(const Thyra::OpTimesLC<Scalar, Node>& x)
    : Thyra::ConstVector<Scalar>(x.formVector())
  {
    cout << "ptr().get() = " << this->ptr().get() << endl;
  }


#ifdef EXTENDED_OPS_ARE_READY
  /*======================================================================
   *
   *    formation of scaled operator 
   *
   *======================================================================*/

  template <class Scalar> inline
  LinearOperator<Scalar> operator*(const Scalar& a, const LinearOperator<Scalar>& A)
  {
    return new ScaledOperator<Scalar>(A, a);
  }
  

  /*======================================================================
   *
   *    composition of operators
   *
   *======================================================================*/

  template <class Scalar> inline
  LinearOperator<Scalar> operator*(const LinearOperator<Scalar>& A, 
                                   const LinearOperator<Scalar>& B)
  {
    return new ComposedOperator<Scalar>(A, B);
  }
#endif

}

#endif
