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

#ifndef THYRA_LINEARCOMBINATION_DECL_HPP
#define THYRA_LINEARCOMBINATION_DECL_HPP

#include "Teuchos_Handle.hpp"
#include "Thyra_ConfigDefs.hpp"
#include "Thyra_LinearOperatorDecl.hpp"
#include "Thyra_VectorHandleOpsDecl.hpp"

namespace Thyra
{
  /** 
   *
   */
  template <class Scalar> 
  class ConvertibleToVector : public virtual Thyra::Converter<Scalar, Thyra::ConstVector<Scalar> >
  {
  public:
    /** */
    virtual ~ConvertibleToVector(){;}

    /** */
    virtual Thyra::Vector<Scalar> formVector() const = 0 ;

    /** */
    Thyra::ConstVector<Scalar> convert() const 
    {
      Thyra::ConstVector<Scalar> rtn = this->formVector();
      return this->formVector();
    }
  };


  /** 
   * Class OpTimesLC holds an operator times something convertible to a vector
   */
  template <class Scalar, class Node>
  class OpTimesLC : public ConvertibleToVector<Scalar>
  {
  public:

    /** */
    virtual ~OpTimesLC(){;}

    /** */
    OpTimesLC(const Scalar& alpha, const Node& x);

    /** */
    OpTimesLC(const Scalar& alpha,
              const Thyra::ConstLinearOperator<Scalar>& op, 
              const Node& x);

    /** 
     * Evaluate the term into the argument vector, overwriting 
     * the previous value of the argument. */
    void evalInto(Thyra::Vector<Scalar>& result) const ;

    /** Add the term into the argument vector */
    void addInto(Thyra::Vector<Scalar>& result, 
                 LCSign sign = LCAdd) const ;

    /** Evaluate the term and return its value */
    virtual Thyra::Vector<Scalar> formVector() const ;

    /** Determine whether this term contains the given vector */
    bool containsVector(const Thyra::VectorBase<Scalar>* vec) const ;

    /** */
     const Thyra::ConstLinearOperator<Scalar>& op() const {return op_;}

    /** */
    const Scalar& alpha() const {return alpha_;}

    /** */
    const Node& node() const {return x_;}

  private:
    Scalar alpha_;
    
    Thyra::ConstLinearOperator<Scalar> op_;

    Node x_;

    /** */
    static Scalar one() {return Teuchos::ScalarTraits<Scalar>::one();}

    /** */
    static Scalar zero() {return Teuchos::ScalarTraits<Scalar>::zero();}
  };


  
  /** */
  template <class Scalar, class Node> inline
  Thyra::Vector<Scalar> formVector(const OpTimesLC<Scalar, Node>& lc)
  {
    return lc.formVector();
  }


  /**
   * Class LC2 is a 2-term linear combination
   */
  template <class Scalar, class Node1, class Node2>
  class LC2  : public ConvertibleToVector<Scalar>
  {
  public:
    /** */
    virtual ~LC2(){;}

    /** */
    LC2(const Node1& x1, const Node2& x2, LCSign sign = LCAdd);

    /** */
    void evalInto(Thyra::Vector<Scalar>& result) const ;

    /** */
    void addInto(Thyra::Vector<Scalar>& result, 
                 LCSign sign = LCAdd) const ;

    /** */
    virtual Thyra::Vector<Scalar> formVector() const ;

    /** */
    bool containsVector(const Thyra::VectorBase<Scalar>* vec) const ;

    
  private:
    Node1 x1_;

    Node2 x2_;

    LCSign sign_;

    /** */
    static Scalar one() {return Teuchos::ScalarTraits<Scalar>::one();}

    /** */
    static Scalar zero() {return Teuchos::ScalarTraits<Scalar>::zero();}
  };

  
  /** */
  template <class Scalar, class Node1, class Node2> inline
  Thyra::Vector<Scalar> formVector(const LC2<Scalar, Node1, Node2>& lc)
  {
    return lc.formVector();
  }
  
  
}

#endif
