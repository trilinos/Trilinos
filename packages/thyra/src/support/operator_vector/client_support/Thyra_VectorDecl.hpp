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

#ifndef THYRA_VECTOR_DECL_HPP
#define THYRA_VECTOR_DECL_HPP

#include "Teuchos_Handle.hpp"
#include "Thyra_ConfigDefs.hpp"
#include "Thyra_VectorBase.hpp"
#include "Thyra_VectorStdOps.hpp"
#include "Thyra_VecOpMacros.hpp"
#include "RTOpPack_Types.hpp"
#include "Teuchos_VerboseObject.hpp"

// namespace ThyraOverloadedOps
// {
//   template <class Scalar, class Node1, class Node2> class LC2;
//   template <class Scalar, class Node> class OpTimesLC; 
//   template <class Scalar> class ConvertibleToVector; 
//   /** 
//    * 
//    */
//   enum LCSign {LCAdd = 1, LCSubtract = -1};
// }


namespace Thyra
{
  template <class Scalar> class Vector;

  template <class Scalar, class Node1, class Node2> class LC2;
  template <class Scalar, class Node> class OpTimesLC; 
  template <class Scalar> class ConvertibleToVector; 
  /** 
   * 
   */
  enum LCSign {LCAdd = 1, LCSubtract = -1};


  template <class Scalar, class TargetType> class Converter
  {
  public:
    /** */
    virtual ~Converter(){;}

    /** */
    virtual TargetType convert() const = 0 ;

    /** */
    virtual void evalInto(Vector<Scalar>& acceptor) const = 0 ;

    /** */
    virtual bool containsVector(const Thyra::VectorBase<Scalar>* vec) const = 0 ;

    /** */
    virtual void addInto(Vector<Scalar>& other, Thyra::LCSign sign) const = 0 ;
  };

  /** 
   *
   */
  template <class Scalar>
  class ConstVector : public virtual Teuchos::ConstHandle<VectorBase<Scalar> >,
                      public virtual Converter<Scalar, ConstVector<Scalar> >
  {
  public:
    TEUCHOS_CONST_HANDLE_CTORS(ConstVector<Scalar>, VectorBase<Scalar>);

    /** Construct a vector from the result of an overloaded operator  */
    ConstVector(const Thyra::ConvertibleToVector<Scalar>& x);

    /** */
    ConstVector<Scalar> evalToConst() const {return *this;}

    /** */
    bool containsVector(const Thyra::VectorBase<Scalar>* vec) const ;

    /** */
    void evalInto(Vector<Scalar>& other) const ;

    /** */
    void addInto(Vector<Scalar>& other, Thyra::LCSign sign) const ;

    /** */
    virtual Scalar operator[](Index globalIndex) const ;

    /** */
    ConstVector<Scalar> convert() const {return *this;}
      
    /** get read-only block */
    ConstVector<Scalar> getBlock(Index i) const ;
    
  };

  /** 
   * \relates ConstVector
   * Return the dimension of the vector 
   */
  template <class Scalar> 
  Index dim(const ConstVector<Scalar>& x) ;


  /** */
  template <class Scalar> 
  std::ostream& operator<<(std::ostream& os, const ConstVector<Scalar>& v);

  /** */
  template <class Scalar> inline 
  ConstVector<Scalar> toVector(const Converter<Scalar, ConstVector<Scalar> >& x) 
  {return x.convert();}

  /** 
   *
   */
  template <class Scalar>
  class Vector : public Teuchos::Handle<VectorBase<Scalar> >,
                 public ConstVector<Scalar>
  {
  public:
    TEUCHOS_HANDLE_CTORS(Vector<Scalar>, VectorBase<Scalar>);

    /** */
    class IndexObject
    {
    public:
      IndexObject(const Teuchos::RefCountPtr<VectorBase<Scalar> >& v, Index i)
        : v_(v), count_(new int), i_(i)
      {
        *count_ = 1;
        val_ = valGotten_ = get_ele(*v_,i_);
      }
      
      IndexObject(const IndexObject& other)
        : v_(other.v_), count_(other.count_), 
          valGotten_(other.valGotten_), val_(other.val_), i_(other.i_)
      {
        *Teuchos::VerboseObjectBase::getDefaultOStream()
          << "IO copy ctor" << endl;
        (*count_)++;
      }

      ~IndexObject()
      {
        if (--(*count_)==0) 
          {
            if( val_ != valGotten_ )
              set_ele( i_, val_, &*v_ );
            delete count_;
          }
      }

      operator Scalar () const {return val_;}

      IndexObject& operator=(const double& value)
      {
        val_ = value;
        return *this;
      }
    private:
      Teuchos::RefCountPtr<VectorBase<Scalar> > v_;
      int* count_;
      Scalar valGotten_;
      Scalar val_;
      Index i_;
      // undefined empty ctor
      IndexObject();
      // undefined assignment op
      IndexObject& operator=(const IndexObject& other);
    };

    

    /** Construct a vector from a 2-term LC */
    template<class Node1, class Node2>
    Vector(const Thyra::LC2<Scalar, Node1, Node2>& x);

    /** Construct a vector from an operator times a linear combination */
    template<class Node>
    Vector(const Thyra::OpTimesLC<Scalar, Node>& x);

    /** Assign a linear combination of vectors to this vector */
    template<class Node1, class Node2>
    Vector& operator=(const Thyra::LC2<Scalar, Node1, Node2>& x);

    /** Assign a scaled linear combination to this vector */
    template<class Node>
    Vector& operator=(const Thyra::OpTimesLC<Scalar, Node>& x);


    /** Write the contents of another vector into this vector */
    Vector<Scalar>& acceptCopyOf(const ConstVector<Scalar>& x);

    /** */
    IndexObject operator[](Index globalIndex)
    {
      //      *Teuchos::VerboseObjectBase::getDefaultOStream()
      // << "calling non-const [] " << endl;
      return IndexObject(this->ptr(), globalIndex);
    }

    /** */
    Scalar operator[](Index globalIndex) const ;

    /** \name Product vector operations */
    //@{
    /** set block  */
    void setBlock(int i, const Vector<Scalar>& v);
      
    /** get modifiable block */
    Vector<Scalar> getBlock(int i);
      
    //@}
  };

  /* copy */
  THYRA_UNARY_VECTOR_OP_DECL(copy, copyInto, assign, "copy");

  /** */
  template <class Scalar> inline
  Thyra::Vector<Scalar> formVector(const Thyra::Vector<Scalar>& x) {return x;}
  
  /** */
  template <class Scalar> inline
  Thyra::Vector<Scalar> formVector(const Thyra::ConstVector<Scalar>& x) 
  {return copy(x);}
}


#endif
