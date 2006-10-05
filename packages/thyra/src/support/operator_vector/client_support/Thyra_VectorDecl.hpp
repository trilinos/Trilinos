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


namespace Thyra
{
  
  template <class Scalar> class VectorSpace;

  template <class Scalar> class Vector;

  template <class Scalar, class Node1, class Node2> class LC2;
  template <class Scalar, class Node> class OpTimesLC; 
  template <class Scalar> class ConvertibleToVector; 

  /** \brief LCSign is used to indicate whether a linear combination object
   * represents addition or subtraction.
   *
   * \ingroup Thyra_Op_Vec_ANA_Development_grp
   */
  enum LCSign {LCAdd = 1, LCSubtract = -1};

  /** \brief Converter that defines the interface for objects that can be
   * converted to vectors.
   *
   * Obviously, vectors can be converted to vectors, but so can linear
   * combinations of vectors or operators times vectors.
   *
   * This interface is key to efficient overloaded operators. Operators do not
   * perform vector operations directly; rather, they construct Converter subtypes
   * (such as LC2 or OpTimesLC) that represent the operation to be performed. The
   * actual operations are carried out only upon either (a) assignment to a vector,   
   * or (b) the Converter is used in a context in which its vector value is required,
   * for instance, when an operation such as a norm is to be performed.
   * 
   * Because overloaded operators must always create and return temporary objects,
   * returning constant-size deferred-evaluation Converter subtypes rather than
   * vectors results in constant-time overhead rather than the \f$O(N)\f$ 
   * overhead that would be incurred with vector return values.
   *
   * \ingroup Thyra_Op_Vec_ANA_Development_grp
   */
  template <class Scalar, class TargetType> class Converter
  {
  public:
    /** \brief . */
    virtual ~Converter(){;}

    /** \brief Convert to the specified target type (e.g., Vector or ConstVector). */
    virtual TargetType convert() const = 0 ;

    /** \brief Evaluate this object, writing the results into the acceptor vector. */
    virtual void evalInto(Vector<Scalar>& acceptor) const = 0 ;

    /** \brief Determine whether this object contains the given vector. */
    virtual bool containsVector(const Thyra::VectorBase<Scalar>* vec) const = 0 ;

    /** \brief Evaluate this object, adding the results into the argument vector. 
     * The sign argument indicates whether this operation is an addition
     * or a subtraction. */
    virtual void addInto(Vector<Scalar>& other, Thyra::LCSign sign) const = 0 ;
  };

  /** \brief Read-only handle class for wrapping <tt>Thyra::VectorBase</tt>
   * objects and allowing for operator-overloading linear algebra.
   *
   * \ingroup Thyra_Op_Vec_ANA_Development_grp
   */
  template <class Scalar>
  class ConstVector : public virtual Teuchos::ConstHandle<VectorBase<Scalar> >,
                      public virtual Converter<Scalar, ConstVector<Scalar> >
  {
  public:

    TEUCHOS_CONST_HANDLE_CTORS(ConstVector<Scalar>, VectorBase<Scalar>);

    /** \brief Construct a vector from the result of an overloaded operator.  */
    ConstVector(const Thyra::ConvertibleToVector<Scalar>& x);

    /** \brief . */
    ConstVector<Scalar> convert() const 
    {
      return *this;
    }

    /** \name Implementation of the Converter interface */
    //@{
    /** \brief . */
    bool containsVector(const Thyra::VectorBase<Scalar>* vec) const ;

    /** \brief . */
    void evalInto(Vector<Scalar>& other) const ;

    /** \brief . */
    void addInto(Vector<Scalar>& other, Thyra::LCSign sign) const ;
    //@}

    /** Element access */
    //@{
    /** \brief Read-only access to an element. */
    virtual Scalar operator[](Index globalIndex) const ;
    //@}

    /** \name Block-related functions */
    //@{
    /** \brief Return number the of blocks in this vector. If the vector is not
     * a product vector, this function will return 1. */
    int numBlocks() const ;
      
    /** \brief Read-only access to the \f$i\f$-th block. If the vector is not
     * a product vector, this function will throw an exception if \f$i\ne 0\f$,
     * or otherwise return the whole vector. */
    ConstVector<Scalar> getBlock(Index i) const ;
    //@}
    
  };

  /** \brief Return the dimension of the vector.
   *
   * \relates ConstVector
   */
  template <class Scalar> 
  Index dim(const ConstVector<Scalar>& x) ;
  
  /** \brief Return the vector space for a vector.
   *
   * \relates ConstVector
   */
  template <class Scalar> inline
  VectorSpace<Scalar> space(const ConstVector<Scalar>& x);

  /** \brief Write to a stream.
   *
   * \relates ConstVector
   */
  template <class Scalar> 
  std::ostream& operator<<(std::ostream& os, const ConstVector<Scalar>& v);

  /* \brief Convert to a ConstVector.
   *
   * \relates ConstVector
   */
  template <class Scalar> inline 
  ConstVector<Scalar> toVector(const Converter<Scalar, ConstVector<Scalar> >& x) 
  {return x.convert();}

  /** \brief Read-write handle class for wrapping <tt>Thyra::VectorBase</tt>
   * objects and allowing for operator-overloading linear algebra.
   *
   * \ingroup Thyra_Op_Vec_ANA_Development_grp
   */
  template <class Scalar>
  class Vector : public Teuchos::Handle<VectorBase<Scalar> >,
                 public ConstVector<Scalar>
  {
  public:

    TEUCHOS_HANDLE_CTORS(Vector<Scalar>, VectorBase<Scalar>);

    /** \brief Construct from a vector space . */
    Vector( const VectorSpace<Scalar> &space );

    /** \brief Allows an element to be changed using <tt>operator=()</tt>. */
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
      /** \brief Writes back the value if it changed. */
      ~IndexObject()
      {
        if (--(*count_)==0) 
          {
            if( val_ != valGotten_ )
              set_ele( i_, val_, &*v_ );
            delete count_;
          }
      }
      /** \brief Implicit conversion to the underlying Scalar. */
      operator Scalar () const {return val_;}
      /** \brief Assignment from a scalar. */
      IndexObject& operator=(const Scalar& value)
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

    /** \brief Construct a vector from a 2-term LC */
    template<class Node1, class Node2>
    Vector(const Thyra::LC2<Scalar, Node1, Node2>& x);

    /** \brief  Construct a vector from an operator times a linear combination */
    template<class Node>
    Vector(const Thyra::OpTimesLC<Scalar, Node>& x);

    /**  \brief Assign a linear combination of vectors to this vector */
    template<class Node1, class Node2>
    Vector& operator=(const Thyra::LC2<Scalar, Node1, Node2>& x);

    /**  \brief Add and assign a linear combination of vectors to this vector */
    template<class Node1, class Node2>
    Vector& operator+=(const Thyra::LC2<Scalar, Node1, Node2>& x);

    /**  \brief Assign a scaled linear combination to this vector */
    template<class Node>
    Vector& operator=(const Thyra::OpTimesLC<Scalar, Node>& x);

    /**  \brief Add and assign a scaled linear combination to this vector */
    template<class Node>
    Vector& operator+=(const Thyra::OpTimesLC<Scalar, Node>& x);

    /** Write the contents of another vector into this vector */
    Vector<Scalar>& acceptCopyOf(const ConstVector<Scalar>& x);

    /** \brief . */
    Scalar operator[](Index globalIndex) const
      {
        return ConstVector<Scalar>::operator[](globalIndex);
      }

    /** \brief Index operator that allows changes to the element.
     *
     * Note: The object returned is of type <tt>IndexObject</tt> which allows
     * for the customary operations to be performed.
     */
    IndexObject operator[](Index globalIndex)
    {
      return IndexObject(this->ptr(), globalIndex);
    }

    /** \name Product vector operations */
    //@{
    /**  \brief set block  */
    void setBlock(int i, const ConstVector<Scalar>& v);

    /**  \brief set block  */
    void setBlock(int i, const Vector<Scalar>& v);
      
    /**  \brief get modifiable block */
    Vector<Scalar> getBlock(int i);
    //@}

  };

  /* copy */
  THYRA_UNARY_VECTOR_OP_DECL(copy, copyInto, assign, "copy");

  /** \brief Form a Vector from this object.
   *
   * For Vector, the operation is simply a pass-through.
   *
   * \relates Vector
   */
  template <class Scalar> inline
  Thyra::Vector<Scalar> formVector(const Thyra::Vector<Scalar>& x) {return x;}
  
  /** \brief Form a Vector from this object.
   *
   * \relates ConstVector
   */
  template <class Scalar> inline
  Thyra::Vector<Scalar> formVector(const Thyra::ConstVector<Scalar>& x) 
  {return copy(x);}

} // namespace Thyra

#endif
