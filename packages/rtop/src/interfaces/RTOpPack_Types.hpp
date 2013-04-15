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
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Roscoe A. Bartlett (rabartl@sandia.gov) 
// 
// ***********************************************************************
// @HEADER


#ifndef RTOPPACK_TYPES_HPP
#define RTOPPACK_TYPES_HPP


#include "RTOp_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_RCP.hpp"
#include "Teuchos_ArrayRCP.hpp"
#include "Teuchos_Range1D.hpp"
#include "Teuchos_ScalarTraits.hpp"
#include "Teuchos_TypeNameTraits.hpp"
#include "Teuchos_Assert.hpp"
#include "Teuchos_implicit_cast.hpp"


namespace RTOpPack {


//
// Basic types
//

/** \brief . */
typedef Teuchos_Ordinal Ordinal;
/** \brief . */
using Teuchos::Ptr;
/** \brief . */
using Teuchos::RCP;
/** \brief . */
using Teuchos::ArrayRCP;
/** \brief . */
using Teuchos::ArrayView;
/** \brief . */
using Teuchos::Array;
/** \brief . */
using Teuchos::Range1D;
/** \brief . */
using Teuchos::ScalarTraits;
/** \brief . */
using Teuchos::TypeNameTraits;

/** \brief . */
typedef Teuchos_Ordinal index_type;
/** \brief . */
typedef char  char_type;


//
// Exceptions
//


/** \brief . */
class UnknownError : public std::logic_error
{public: UnknownError(const std::string& what_arg) : std::logic_error(what_arg) {}};
/** \brief . */
class InvalidUsage : public std::logic_error
{public: InvalidUsage(const std::string& what_arg) : std::logic_error(what_arg) {}};
/** \brief . */
class InvalidNumVecs : public std::logic_error
{public: InvalidNumVecs(const std::string& what_arg) : std::logic_error(what_arg) {}};
/** \brief . */
class InvalidNumTargVecs : public std::logic_error
{public: InvalidNumTargVecs(const std::string& what_arg) : std::logic_error(what_arg) {}};
/** \brief . */
class IncompatibleVecs : public std::logic_error
{public: IncompatibleVecs(const std::string& what_arg) : std::logic_error(what_arg) {}};
/** \brief . */
class IncompatibleReductObj : public std::logic_error
{public: IncompatibleReductObj(const std::string& what_arg) : std::logic_error(what_arg) {}};


//
// VectorBase subviews
//


/** \brief Class for a non-changeable sub-vector.
 *
 * For a sub-vector <tt>vec</tt>, the corresponding entries in the global
 * vector <tt>x(j)</tt> (one based) are as follows:

 \verbatim

   x( vec.globalOffset() + k ) = v(k), for k = 0...vec.subDim()-1

 \endverbatim

 * The stride <tt>vec.stride()</tt> may be positive (>0) or negative (<0) but
 * not zero (0).  A negative stride <tt>vec.stride() < 0</tt> allows a reverse
 * traversal of the elements.
 *
 * <b>WARNING!</b> the default copy constructor and assignment operators are
 * allowed which results in only pointer copy, not deep copy!  You have been
 * warned!
 */
template<class Scalar>
class ConstSubVectorView {
public:
  /** \brief . */
  ConstSubVectorView() : globalOffset_(0), subDim_(0), stride_(0) {}
  /** \brief . */
  ConstSubVectorView(const ArrayRCP<const Scalar> &values_in)
    :globalOffset_(0), subDim_(0), stride_(0)
    { initialize(0, values_in.size(), values_in, 1); }
  /** \brief . */
  ConstSubVectorView(Ordinal globalOffset_in, Ordinal subDim_in,
    const ArrayRCP<const Scalar> &values_in, ptrdiff_t stride_in)
    :globalOffset_(0), subDim_(0), stride_(0)
    { initialize(globalOffset_in, subDim_in, values_in, stride_in); }
  /** \brief . */
  ConstSubVectorView( const ConstSubVectorView<Scalar>& sv )
    :globalOffset_(sv.globalOffset()), subDim_(sv.subDim()),
     values_(sv.values()), stride_(sv.stride()) 
    {}
  /** \brief . */
  void initialize(Ordinal globalOffset_in, Ordinal subDim_in,
    const ArrayRCP<const Scalar> &values_in, ptrdiff_t stride_in)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT(globalOffset_in >= 0);
      if (!is_null(values_in)) {
        TEUCHOS_ASSERT(subDim_in >= 0);
        TEUCHOS_ASSERT(stride_in != 0);
        TEUCHOS_ASSERT(
          subDim_in*std::abs(Teuchos::as<int>(stride_in)) - 1 <= values_in.upperOffset());
        TEUCHOS_ASSERT(values_in.lowerOffset() <= 0);
      }
      else {
        TEUCHOS_ASSERT(subDim_in==0);
      }
#endif
      globalOffset_=globalOffset_in;
      subDim_=subDim_in;
      values_=values_in;
      stride_=stride_in;
    }
  /** \brief . */
  void uninitialize()
    { globalOffset_ = 0; subDim_=0; values_ = Teuchos::null; stride_ = 0; }
  /** \brief . */
  void setGlobalOffset(Ordinal globalOffset_in)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT(globalOffset_in >= 0);
#endif
      globalOffset_ = globalOffset_in;
    } 
  /** \brief . */
  Ordinal globalOffset() const { return globalOffset_; }
  /** \brief . */
  Ordinal subDim() const { return subDim_; }
  /** \brief . */
  const ArrayRCP<const Scalar>  values() const { return values_; }
  /** \brief . */
  ptrdiff_t stride() const { return stride_; }
  /** \brief Zero-based indexing (Preconditions: <tt>values()!=NULL && (0 <= i
   * < subDim())</tt>). */
  const Scalar& operator[](Ordinal i) const
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(i, 0, subDim_);
#endif
      return valuesBegin()[stride_*i];
    }
  /** \brief Zero-based indexing (Preconditions: <tt>values()!=NULL && (0 <= i
   * < subDim())</tt>). */
  const Scalar& operator()(Ordinal i) const { return (*this)[i]; }
private:
  Ordinal globalOffset_;
  Ordinal subDim_;
  ArrayRCP<const Scalar> values_;
  ptrdiff_t stride_;
  const typename ArrayRCP<const Scalar>::iterator valuesBegin() const
    {
      if (stride_ > 0)
        return values_.begin();
      return values_.begin() + (subDim_*std::abs(Teuchos::as<int>(stride_)) - 1);
    } 
public:
};


/** \brief Class for a changeable sub-vector.
 *
 * This class derives from <tt>ConstSubVectorView</tt> and adds methods to
 * change the data.  Note, a <tt>const SubVectorView</tt> object allows
 * clients to change the values in the underlying subvector.  The meaning of
 * <tt>const</tt> in this context is that the view of the data can not change.
 *
 * <b>WARNING!</b> the default copy constructor and assignment operators are
 * allowed which results in only pointer copy, not deep copy.  This means this
 * class has shallow copy semantics. You have been warned!
 *
 * NOTE: It is perfectly safe to derive this class from ConstSubVectorView
 * even through it does not have a virtual destructor.  That is because this
 * derived class has no data members that would cause problems in slicing or
 * memory leaks when deleting.
 */
template<class Scalar>
class SubVectorView : public ConstSubVectorView<Scalar> {
public:
  /** \brief . */
  SubVectorView() {}
  /** \brief . */
  SubVectorView(const ArrayRCP<Scalar> &values_in)
    :ConstSubVectorView<Scalar>(values_in)
    {}
  /** \brief . */
  SubVectorView(Ordinal globalOffset_in, Ordinal subDim_in,
    const ArrayRCP<Scalar> &values_in, ptrdiff_t stride_in)
    :ConstSubVectorView<Scalar>(globalOffset_in, subDim_in, values_in, stride_in)
    {}
  /** \brief . */
  SubVectorView(Ordinal subDim_in)
    :ConstSubVectorView<Scalar>(0, subDim_in, Teuchos::arcp<Scalar>(subDim_in), 1)
    {}
  /** \brief . */
  SubVectorView(const SubVectorView<Scalar> & sv)
    :ConstSubVectorView<Scalar>(sv)
    {}
  /** \brief . */
  void initialize(Ordinal globalOffset_in, Ordinal subDim_in,
    const ArrayRCP<Scalar> &values_in, ptrdiff_t stride_in)
    { ConstSubVectorView<Scalar>::initialize(globalOffset_in, subDim_in, values_in, stride_in); }
  /** \brief . */
  const ArrayRCP<Scalar> values() const
    { return Teuchos::arcp_const_cast<Scalar>(ConstSubVectorView<Scalar>::values());  }
  /** \brief Zero-based indexing (Preconditions: <tt>values()!=NULL && (0 <= i
   * < subDim())</tt>). */
  Scalar& operator[](Ordinal i) const
    { return const_cast<Scalar&>(ConstSubVectorView<Scalar>::operator[](i)); }
  /** \brief Zero-based indexing (Preconditions: <tt>values()!=NULL && (0 <= i
   * < subDim())</tt>). */
  Scalar& operator()(Ordinal i) const { return (*this)[i]; }
public:
};


/** \brief . */
template<class Scalar>
void assign_entries( const Ptr<const SubVectorView<Scalar> > &msv,
  const ConstSubVectorView<Scalar> &sv )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(msv->subDim(), sv.subDim());
#endif
  for( int i = 0; i < sv.subDim(); ++i ) {
    (*msv)(i) = sv(i);
  }
}


/** \brief .
 *
 * \relates ConstSubVectorView
 */
template<class Scalar>
std::ostream& operator<<(std::ostream &out, const ConstSubVectorView<Scalar> &sv)
{
  out
    << "{"
    << "globalOffset="<<sv.globalOffset()
    << ",subDim="<<sv.subDim()
    << ",values="<<sv.values()
    << ",stride="<<sv.stride()
    << "}";
  return out;
}


//
// MultiVectorBase subviews
//


/** \brief Class for a non-changeable sub-multi-vector (submatrix).
 *
 * For a sub-multi-vector <tt>mv</tt>, the corresponding entries in the global
 * multi-vector <tt>X(j)</tt> (one based) are as follows:

 \verbatim

   X(mv.globalOffset()+k1,mv.colOffset()+k2) = mv(k1,k2),
       for k1 = 0...mv.subDim()-1, k2 = 0...mv.numSubCols()-1

 \endverbatim

 * Unlike vectors, there can only be a unit stride between vector elements in
 * a particular column and there is a Fortran-like leading dimension
 * <tt>mv.leadingDim()</tt> that separates corresponding elements in each
 * column sub-vector.
 *
 * <b>WARNING!</b> the default copy constructor and assignment operators are
 * allowed which results in only pointer copy, not deep copy!  You have been
 * warned!
 */
template<class Scalar>
class ConstSubMultiVectorView {
public:
  /** \brief . */
  ConstSubMultiVectorView()
    :globalOffset_(0), subDim_(0), colOffset_(0), numSubCols_(0),
     leadingDim_(0)
    {}
  /** \brief . */
  ConstSubMultiVectorView(
    Ordinal globalOffset_in, Ordinal subDim_in,
    Ordinal colOffset_in, Ordinal numSubCols_in,
    const ArrayRCP<const Scalar> &values_in, Ordinal leadingDim_in
    )
    :globalOffset_(0), subDim_(0), colOffset_(0), numSubCols_(0),
     leadingDim_(0)
    {
      initialize(globalOffset_in, subDim_in, colOffset_in, numSubCols_in, values_in,
        leadingDim_in);
    }
  /** \brief . */
  ConstSubMultiVectorView( const ConstSubMultiVectorView<Scalar>& smv )
    :globalOffset_(smv.globalOffset()), subDim_(smv.subDim()),
     colOffset_(smv.colOffset()), numSubCols_(smv.numSubCols()),
     values_(smv.values()), leadingDim_(smv.leadingDim())
    {}
  /** \brief . */
  void initialize(
    Ordinal globalOffset_in, Ordinal subDim_in,
    Ordinal colOffset_in, Ordinal numSubCols_in,
    const ArrayRCP<const Scalar> &values_in, Ordinal leadingDim_in
    )
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT(globalOffset_in >= 0);
      TEUCHOS_ASSERT(colOffset_in >= 0);
      if (!is_null(values_in)) {
        TEUCHOS_ASSERT(subDim_in >= 0);
        TEUCHOS_ASSERT(leadingDim_in >= subDim_in);
        TEUCHOS_ASSERT(numSubCols_in*leadingDim_in - 1 <= values_in.upperOffset());
        TEUCHOS_ASSERT(values_in.lowerOffset() <= 0);
      }
      else {
        TEUCHOS_ASSERT(subDim_in == 0);
      }
#endif
      globalOffset_=globalOffset_in;
      subDim_=subDim_in;
      colOffset_=colOffset_in;
      numSubCols_=numSubCols_in;
      values_=values_in;
      leadingDim_=leadingDim_in;
    }
  /** \brief . */
  void uninitialize()
    {
      globalOffset_ = 0; subDim_=0; colOffset_=0, numSubCols_=0;
      values_=Teuchos::null; leadingDim_=0;
    }
  /** \brief . */
  void setGlobalOffset(Ordinal globalOffset_in)
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT(globalOffset_in >= 0);
#endif
      globalOffset_ = globalOffset_in;
    } 
  /** \brief . */
  Ordinal globalOffset() const { return globalOffset_; }
  /** \brief . */
  Ordinal subDim() const { return subDim_; }
  /** \brief . */
  Ordinal colOffset() const { return colOffset_; }
  /** \brief . */
  Ordinal numSubCols() const { return numSubCols_; }
  /** \brief . */
  const ArrayRCP<const Scalar> values() const { return values_; }
  /** \brief . */
  Ordinal leadingDim() const { return leadingDim_; }
  /** \brief Zero-based indexing (Preconditions: <tt>values()!=NULL &&
   * (0<=i<subDim()) && (0<=j< numSubCols()</tt>).
   */
  const Scalar& operator()(Ordinal i, Ordinal j) const
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(i, 0, subDim_);
      TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(j, 0, numSubCols_ );
#endif
      return values_[ i + leadingDim_*j ];
    }
  /** \brief Return a <tt>ConstSubVectorView</tt> view of the jth sub-column
   * (Preconditions: <tt>values()!=NULL && (0<=j<numSubCols()</tt>).
   */
  ConstSubVectorView<Scalar> col( const Ordinal j ) const
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(j, 0, numSubCols_ );
#endif
      return ConstSubVectorView<Scalar>(
        globalOffset(), subDim(), values().persistingView(j*leadingDim(),subDim()), 1 );
    }
private:
  Ordinal globalOffset_;
  Ordinal subDim_;
  Ordinal colOffset_;
  Ordinal numSubCols_;
  ArrayRCP<const Scalar> values_;
  Ordinal leadingDim_;
public:
};


/** \brief Class for a changeable sub-vector.
 *
 * This class derives from <tt>ConstSubVectorView</tt> and adds methods to
 * change the data.  Note, a <tt>const SubVectorView</tt> object allows
 * clients to change the values in the underlying subvector.  The meaning of
 * <tt>const</tt> in this context is that the view of the data can not change.
 *
 * <b>WARNING!</b> the default copy constructor and assignment operators are
 * allowed which results in only pointer copy, not deep copy!  You have been
 * warned!
 *
 * NOTE: It is perfectly safe to derive this class from
 * ConstSubMultiVectorView even through it does not have a virtual destructor.
 * That is because this derived class has no data members that would cause
 * problems in slicing or memory leaks when deleting.
 */
template<class Scalar>
class SubMultiVectorView : public ConstSubMultiVectorView<Scalar> {
public:
  /** \brief . */
  SubMultiVectorView() {}
  /** \brief . */
  SubMultiVectorView(
    Ordinal numRows_in, Ordinal numCols_in
    )
    :ConstSubMultiVectorView<Scalar>(0, numRows_in, 0, numCols_in,
      Teuchos::arcp<Scalar>(numRows_in*numCols_in), numRows_in)
    {}
  /** \brief . */
  SubMultiVectorView(
    Ordinal globalOffset_in, Ordinal subDim_in,
    Ordinal colOffset_in, Ordinal numSubCols_in,
    const ArrayRCP<Scalar> &values_in, Ordinal leadingDim_in
    )
    :ConstSubMultiVectorView<Scalar>(globalOffset_in, subDim_in,
      colOffset_in, numSubCols_in, values_in, leadingDim_in)
    {}
  /** \brief . */
  SubMultiVectorView( const SubMultiVectorView<Scalar> & smv)
    :ConstSubMultiVectorView<Scalar>(smv)
    {}
  /** \brief . */
 void initialize(
   Ordinal globalOffset_in, Ordinal subDim_in,
   Ordinal colOffset_in, Ordinal numSubCols_in,
   const ArrayRCP<Scalar> &values_in, Ordinal leadingDim_in
   )
   {
     ConstSubMultiVectorView<Scalar>::initialize(globalOffset_in,
       subDim_in, colOffset_in, numSubCols_in, values_in, leadingDim_in);
   }
  /** \brief . */
  const ArrayRCP<Scalar> values() const
    {
      return Teuchos::arcp_const_cast<Scalar>(
        ConstSubMultiVectorView<Scalar>::values());
    }
  /** \brief Zero-based indexing (Preconditions: <tt>values()!=NULL && (0<=i<
   * subDim()) && (0<=j<numSubCols()</tt>).
   */
  Scalar& operator()(Ordinal i, Ordinal j) const
    { return const_cast<Scalar&>(ConstSubMultiVectorView<Scalar>::operator()(i,j)); }
  /** \brief Return a <tt>SubVectorView</tt> view of the jth sub-column
   * (Preconditions: <tt>values()!=NULL && && (0<=j<numSubCols()</tt>).
   */
  SubVectorView<Scalar> col( const Ordinal j ) const
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(j, 0, this->numSubCols());
#endif
      return SubVectorView<Scalar>(this->globalOffset(), this->subDim(),
        values().persistingView(j*this->leadingDim(),this->subDim()), 1);
    }
public:
};


/** \brief . */
template<class Scalar>
void assign_entries( const Ptr<const SubMultiVectorView<Scalar> > &msmv,
  const ConstSubMultiVectorView<Scalar> &smv )
{
#ifdef TEUCHOS_DEBUG
  TEUCHOS_ASSERT_EQUALITY(msmv->subDim(), smv.subDim());
  TEUCHOS_ASSERT_EQUALITY(msmv->numSubCols(), smv.numSubCols());
#endif
  for( Ordinal j = 0; j < smv.numSubCols(); ++j ) {
    for( Ordinal i = 0; i < smv.subDim(); ++i ) {
      (*msmv)(i,j) = smv(i,j);
    }
  }
}


//
// Primitive Type Traits
//


/** \brief A templated traits class for decomposing object into an
 * array of primitive objects.
 *
 * The idea behind this traits class it that it allows an object of
 * semi-complex structure to be externalized into arrays of primitive data
 * types.
 *
 * This default traits class works just fine for types that are
 * already primitive.
 */
template <class Scalar, class ConcreteObj>
class PrimitiveTypeTraits {
public:
  /** \brief . */
  typedef Scalar primitiveType;
  /** \brief . */
  static int numPrimitiveObjs()
    { return Scalar::this_type_is_missing_a_specialization(); }
  /** \brief . */
  static int numIndexObjs()
    { return Scalar::this_type_is_missing_a_specialization(); }
  /** \brief . */
  static int numCharObjs()
    { return Scalar::this_type_is_missing_a_specialization(); }
  /** \brief . */
  static void extractPrimitiveObjs(
    const Scalar &obj,
    const ArrayView<primitiveType> &primitiveObjs,
    const ArrayView<index_type> &indexObjs,
    const ArrayView<char> &charObjs
    )
    {
      Scalar::this_type_is_missing_a_specialization(obj);
    }
  /** \brief . */
  static void loadPrimitiveObjs(
    const ArrayView<const primitiveType> &primitiveObjs,
    const ArrayView<const index_type> &indexObjs,
    const ArrayView<const char> &charObjs,
    const Ptr<Scalar> &obj
    )
    {
      *obj = Scalar::this_type_is_missing_a_specialization();
    }
};



/** \brief Specialization where the scalar type is the same as the concrete
 * object type.
 */
template <class Scalar>
class PrimitiveTypeTraits<Scalar, Scalar> {
public:
  /** \brief . */
  typedef Scalar primitiveType;
  /** \brief . */
  static int numPrimitiveObjs() { return 1; }
  /** \brief . */
  static int numIndexObjs() { return 0; }
  /** \brief . */
  static int numCharObjs() { return 0; }
  /** \brief . */
  static void extractPrimitiveObjs(
    const Scalar &obj,
    const ArrayView<primitiveType> &primitiveObjs,
    const ArrayView<index_type> &indexObjs,
    const ArrayView<char> &charObjs
    )
    {
      assertInput(primitiveObjs, indexObjs, charObjs);
      primitiveObjs[0] = obj;
    }
  /** \brief . */
  static void loadPrimitiveObjs(
    const ArrayView<const primitiveType> &primitiveObjs,
    const ArrayView<const index_type> &indexObjs,
    const ArrayView<const char> &charObjs,
    const Ptr<Scalar> &obj
    )
    {
      assertInput(primitiveObjs, indexObjs, charObjs);
      *obj = primitiveObjs[0];
    }
private:
  static void assertInput(
    const ArrayView<const primitiveType> &primitiveObjs,
    const ArrayView<const index_type> &indexObjs,
    const ArrayView<const char> &charObjs
    )
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT( primitiveObjs.size()!=1 || indexObjs.size()!=0
        || charObjs.size()!=0 );
#endif
    }
};


/** \brief Specialization for index_type concrete object. */
template <class Scalar>
class PrimitiveTypeTraits<Scalar, index_type> {
public:
  /** \brief . */
  typedef PrimitiveTypeTraits<Scalar,Scalar> ScalarPrimitiveTypeTraits;
  /** \brief . */
  typedef typename ScalarPrimitiveTypeTraits::primitiveType primitiveType;
  /** \brief . */
  static int numPrimitiveObjs() { return 0; }
  /** \brief . */
  static int numIndexObjs() { return 1; }
  /** \brief . */
  static int numCharObjs() { return 0; }
  /** \brief . */
  static void extractPrimitiveObjs(
    const index_type &obj,
    const ArrayView<primitiveType> &primitiveObjs,
    const ArrayView<index_type> &indexObjs,
    const ArrayView<char> &charObjs
    )
    {
      assertInput(primitiveObjs, indexObjs, charObjs);
      indexObjs[0] = obj;
    }
  /** \brief . */
  static void loadPrimitiveObjs(
    const ArrayView<const primitiveType> &primitiveObjs,
    const ArrayView<const index_type> &indexObjs,
    const ArrayView<const char> &charObjs,
    const Ptr<index_type> &obj
    )
    {
      assertInput(primitiveObjs, indexObjs, charObjs);
      *obj = indexObjs[0];
    }
private:
  static void assertInput(
    const ArrayView<const primitiveType> &primitiveObjs,
    const ArrayView<const index_type> &indexObjs,
    const ArrayView<const char> &charObjs
    )
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT( primitiveObjs.size()!=0 || indexObjs.size()!=1
        || charObjs.size()!=0 );
#endif
    }
};


#if defined(HAVE_COMPLEX) && defined(HAVE_TEUCHOS_COMPLEX)


/** \brief Partial specialization of <tt>PrimitiveTypeTraits</tt> for
 * <tt>std::complex<Scalar> scalar type and reduction type</tt>.
 */
template <class Scalar>
class PrimitiveTypeTraits<std::complex<Scalar>, std::complex<Scalar> > {
public:
  /** \brief . */
  typedef PrimitiveTypeTraits<Scalar,Scalar> ScalarPrimitiveTypeTraits;
  /** \brief . */
  typedef typename ScalarPrimitiveTypeTraits::primitiveType primitiveType;
  /** \brief . */
  static int numPrimitiveObjs()
    { return 2*ScalarPrimitiveTypeTraits::numPrimitiveObjs(); }
  /** \brief . */
  static int numIndexObjs() { return 0; }
  /** \brief . */
  static int numCharObjs() { return 0; }
  /** \brief . */
  static void extractPrimitiveObjs(
    const std::complex<Scalar> &obj,
    const ArrayView<primitiveType> &primitiveObjs,
    const ArrayView<index_type> &indexObjs,
    const ArrayView<char> &charObjs
    )
    {
      using Teuchos::null;
      const int numScalarPrimitiveObjs =
        ScalarPrimitiveTypeTraits::numPrimitiveObjs();
      assertInput(primitiveObjs, indexObjs, charObjs);
      ScalarPrimitiveTypeTraits::extractPrimitiveObjs(
        obj.real(), primitiveObjs(0,numScalarPrimitiveObjs), null, null );
      ScalarPrimitiveTypeTraits::extractPrimitiveObjs(
        obj.imag(), primitiveObjs(numScalarPrimitiveObjs,numScalarPrimitiveObjs), null, null );
    }
  /** \brief . */
  static void loadPrimitiveObjs(
    const ArrayView<const primitiveType> &primitiveObjs,
    const ArrayView<const index_type> &indexObjs,
    const ArrayView<const char> &charObjs,
    const Ptr<std::complex<Scalar> > &obj
    )
    {
      using Teuchos::null;
      using Teuchos::outArg;
      assertInput(primitiveObjs, indexObjs, charObjs);
      const int numScalarPrimitiveObjs =
        ScalarPrimitiveTypeTraits::numPrimitiveObjs();
      Scalar real, imag;
      ScalarPrimitiveTypeTraits::loadPrimitiveObjs(
        primitiveObjs(0,numScalarPrimitiveObjs), null, null,
        outArg(real) );
      ScalarPrimitiveTypeTraits::loadPrimitiveObjs(
        primitiveObjs(numScalarPrimitiveObjs,numScalarPrimitiveObjs), null, null,
        outArg(imag) );
      *obj = std::complex<Scalar>( real, imag );
    }
private:
  static void assertInput(
    const ArrayView<const primitiveType> &primitiveObjs,
    const ArrayView<const index_type> &indexObjs,
    const ArrayView<const char> &charObjs
    )
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT(
        primitiveObjs.size()!=2*ScalarPrimitiveTypeTraits::numPrimitiveObjs()
        || indexObjs.size()!=0
        || charObjs.size()!=0 );
#endif
    }
};


/** \brief Partial specialization of <tt>PrimitiveTypeTraits</tt> for
 * <tt>std::complex<Scalar> scalar type and Scalar reduction type</tt>.
 */
template <class Scalar>
class PrimitiveTypeTraits<std::complex<Scalar>, Scalar> {
public:
  /** \brief . */
  typedef PrimitiveTypeTraits<Scalar,Scalar> ScalarPrimitiveTypeTraits;
  /** \brief . */
  typedef typename ScalarPrimitiveTypeTraits::primitiveType primitiveType;
  /** \brief . */
  static int numPrimitiveObjs()
    { return ScalarPrimitiveTypeTraits::numPrimitiveObjs(); }
  /** \brief . */
  static int numIndexObjs() { return 0; }
  /** \brief . */
  static int numCharObjs() { return 0; }
  /** \brief . */
  static void extractPrimitiveObjs(
    const Scalar &obj,
    const ArrayView<primitiveType> &primitiveObjs,
    const ArrayView<index_type> &indexObjs,
    const ArrayView<char> &charObjs
    )
    {
      using Teuchos::null;
      assertInput(primitiveObjs, indexObjs, charObjs);
      ScalarPrimitiveTypeTraits::extractPrimitiveObjs(
        obj, primitiveObjs, null, null );
    }
  /** \brief . */
  static void loadPrimitiveObjs(
    const ArrayView<const primitiveType> &primitiveObjs,
    const ArrayView<const index_type> &indexObjs,
    const ArrayView<const char> &charObjs,
    const Ptr<Scalar > &obj
    )
    {
      using Teuchos::null;
      assertInput(primitiveObjs, indexObjs, charObjs);
      ScalarPrimitiveTypeTraits::loadPrimitiveObjs(
        primitiveObjs, null, null, obj );
    }
private:
  static void assertInput(
    const ArrayView<const primitiveType> &primitiveObjs,
    const ArrayView<const index_type> &indexObjs,
    const ArrayView<const char> &charObjs
    )
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_TEST_FOR_EXCEPT(
        primitiveObjs.size()!=ScalarPrimitiveTypeTraits::numPrimitiveObjs()
        || indexObjs.size()!=0
        || charObjs.size()!=0 );
#endif
    }
};


#endif // defined(HAVE_COMPLEX) && defined(HAVE_TEUCHOS_COMPLEX)



//
// Forward declaration for templated types
//


/** \brief . */
template<class Scalar>  class RTOpT;


} // namespace RTOpPack


#endif // RTOPPACK_TYPES_HPP
