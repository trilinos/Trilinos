// @HEADER
// ***********************************************************************
// 
//      Thyra: Interfaces and Support Code for the Interoperability of Abstract Numerical Algorithms 
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

#ifndef RTOPPACK_TYPES_HPP
#define RTOPPACK_TYPES_HPP

#include "RTOp_ConfigDefs.hpp"
#include "Teuchos_Array.hpp"
#include "Teuchos_TestForException.hpp"

namespace RTOpPack {

//
// Basic types
//

/** \brief Depreciated . */
typedef Teuchos_Index index_type;
/** \brief Depreciated . */
typedef char  char_type;
/** \brief . */
typedef index_type Index;

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

//
// VectorBase subviews
//

/** \brief Class for a non-mutable sub-vector.
 *
 * For a sub-vector <tt>vec</tt>, the corresponding entries
 *	in the global vector <tt>x(j)</tt> (one based) are as follows:
 \verbatim

  x( vec.globalOffset() + k ) = v(k), for k = 0,...,vec.subDim()-1
  \endverbatim
 * The stride <tt>vec.stride()</tt> may be positive (>0), negative (<0)
 * or even zero (0).  A negative stride <tt>vec.stride() < 0</tt> allows a
 * reverse traversal of the elements.  A zero stride
 * <tt>vec.stride()</tt> allows a sub-vector with all the elements the same.
 *
 * The raw pointer to the start of the memory can be obtained as
 * <tt>&vec(0)</tt>.
 *
 * Warning! the default copy constructor and assignment operators are
 * allowed which results in only pointer copy, not deep copy!  You
 * have been warned!
 */
template<class Scalar>
class ConstSubVectorView {
public:
  /** \brief . */
  ConstSubVectorView() : globalOffset_(0), subDim_(0), values_(NULL), stride_(0) {}
  /** \brief . */
  ConstSubVectorView(Teuchos_Index globalOffset, Teuchos_Index subDim, const Scalar *values, ptrdiff_t stride)
    :globalOffset_(globalOffset), subDim_(subDim), values_(values), stride_(stride) 
    {}
  /** \brief . */
  ConstSubVectorView( const ConstSubVectorView<Scalar>& sv )
    :globalOffset_(sv.globalOffset()), subDim_(sv.subDim()), values_(sv.values()), stride_(sv.stride()) 
    {}
  /** \brief . */
  void initialize(Teuchos_Index globalOffset, Teuchos_Index subDim, const Scalar *values, ptrdiff_t stride)
    { globalOffset_=globalOffset; subDim_=subDim; values_=values; stride_=stride;  }
  /** \brief . */
  void set_uninitialized()
    { globalOffset_ = 0; subDim_=0; values_=NULL; stride_ = 0; }
  /** \brief . */
  void setGlobalOffset(Teuchos_Index globalOffset) { globalOffset_ = globalOffset; } 
  /** \brief . */
  Teuchos_Index   globalOffset() const { return globalOffset_; }
  /** \brief . */
  Teuchos_Index   subDim()       const { return subDim_;  }
  /** \brief . */
  const Scalar*     values()       const { return values_;  }
  /** \brief . */
  ptrdiff_t         stride()       const { return stride_;  }
  /// Zero-based indexing (Preconditions: <tt>values()!=NULL && (0 <= i < subDim())</tt>)
  const Scalar& operator[](Teuchos_Index i) const
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        !( 0 <= i && i < subDim_ ), std::logic_error
        ,"Error, index i="<<i<<" does not fall in the range [0,"<<(subDim_-1)<<"]!"
        );
#endif
      return values_[ stride_*i ];
    }
  /// Zero-based indexing (Preconditions: <tt>values()!=NULL && (0 <= i < subDim())</tt>)
  const Scalar& operator()(Teuchos_Index i) const { return (*this)[i]; }
private:
  Teuchos_Index     globalOffset_;
  Teuchos_Index     subDim_;
  const Scalar      *values_;
  ptrdiff_t         stride_;
};

/** \brief Class for a mutable sub-vector.
 *
 * This class derives from <tt>ConstSubVectorView</tt> and adds methods to
 * mutate the data.  Note, a <tt>const SubVectorView</tt> object
 * allows clients to change the values in the underlying subvector.
 * The meaning of <tt>const</tt> in this context is that the
 * view of the data can not change.
 *
 * Warning! the default copy constructor and assignment operators are
 * allowed which results in only pointer copy, not deep copy!  You
 * have been warned!
 */
template<class Scalar>
class SubVectorView : public ConstSubVectorView<Scalar> {
public:
  /** \brief . */
  SubVectorView() {}
  /** \brief . */
  SubVectorView(Teuchos_Index globalOffset, Teuchos_Index subDim, Scalar *values, ptrdiff_t stride)
    :ConstSubVectorView<Scalar>(globalOffset, subDim, values, stride)
    {}
  /** \brief . */
  SubVectorView( const SubVectorView<Scalar> & sv)
    :ConstSubVectorView<Scalar>(sv)
    {}
  /** \brief . */
  void initialize(Teuchos_Index globalOffset, Teuchos_Index subDim, Scalar *values, ptrdiff_t stride)
    { ConstSubVectorView<Scalar>::initialize(globalOffset, subDim, values, stride); }
  /** \brief . */
  void set_uninitialized()
    { ConstSubVectorView<Scalar>::set_uninitialized(); }
  /** \brief . */
  Scalar* values() const { return const_cast<Scalar*>(ConstSubVectorView<Scalar>::values());  }
  /// Zero-based indexing (Preconditions: <tt>values()!=NULL && (0 <= i < subDim())</tt>)
  Scalar& operator[](Teuchos_Index i) const { return const_cast<Scalar&>(ConstSubVectorView<Scalar>::operator[](i)); } // Is range changed in subclass!
  /// Zero-based indexing (Preconditions: <tt>values()!=NULL && (0 <= i < subDim())</tt>)
  Scalar& operator()(Teuchos_Index i) const { return (*this)[i]; }
};

/** \brief Class for a (sparse or dense) sub-vector.
 *
 *	Sparse and dense local vectors are supported as follows:
  *
  *	A dense vector <tt>vec</tt> is identified by <tt>vec.subDim() == vec.subNz()</tt>
  * and <tt>vec.indices() == NULL</tt> in which case
  *	<tt>vec.indicesStride()</tt>, <tt>vec.localOffset()</tt> and <tt>vec.isSorted()</tt>
  * are ignored.  For a dense sub-vector <tt>vec</tt>, the corresponding entries
 *	in the global vector <tt>x(j)</tt> (one based) are as follows:
 \verbatim

  x( vec.globalOffset() + k )
    = vec.values()[ vec.valueStride() * k ]

  for k = 0,...,vec.subDim()-1
 \endverbatim
 * The stride member <tt>vec.valueStride()</tt> may be positive (>0), negative (<0)
 * or even zero (0).  A negative stride <tt>vec.valueStride() < 0</tt> allows a
 * reverse traversal of the elements in <tt>vec.values()[]</tt>.  A zero stride
 * <tt>vec.valueStride() == 0</tt> allows a vector with all the elements the same.
 *
 *	A sparse vector is identified by <tt>vec.subDim() > vec.subNz()</tt>
 * or <tt>vec.indices() != NULL</tt>
 * in which case all the fields in the structure are meaningful.
 *	The corresponding elements in the global vector <tt>x(j)</tt>
 * defined as:
 \verbatim

  x( vec.globalOffset() + vec.localOffset() + vec.indices()[vec.indicesStride()*k] )
    = vec.values[vec.valueStride()*k]

  for k = 0,...,vec.subNz()-1
 \endverbatim
 * If <tt>vec.subNz() == 0</tt> then it is allowed for <tt>vec.indices() == NULL</tt>.
 * If <tt>vec.subDim() > vec.subNz() > 0</tt> then <tt>vec.indices() != NULL</tt> must be true.
 *
  * A sparse sub-vector may be sorted (<tt>vec.isSorted()!=0</tt>) or
  * unsorted (<tt>vec.isSorted()==0</tt>) but the indices <tt>vec.indices()[k]</tt>
  * must be unique.  A sorted vector (<tt>vec.isSorted()!=0</tt>) means that
  * the indices are in ascending order:
  \verbatim

  vec.indices()[vec.indicesStride()*k] < vec.indices()[vec.indicesStride()*k]

  for k = 0,...,vec.subNz()-1
 \endverbatim
 * The member <tt>vec.localOffset()</tt> is used to shift the values in <tt>vec.indices()[]</tt>
 * to be in range of the local sub-vector.  In other words:
 \verbatim
  
  0 <= vec.localOffset() + vec.indices()[vec.indicesStride()*k] <= vec.subNz()

  for k = 0...vec.subNz()-1
 \endverbatim
 * The member <tt>vec.valueStride()</tt> may be positive (>0), negative (<0) or zero (0).
 * However, the member <tt>vec.indicesStride()</tt> may be may be positive (>0)
 * or negative (<0) but not zero (0).  Allowing <tt>vec.indicesStride() == 0</tt>
 * would mean that a vector would have <tt>vec.subNz()</tt> nonzero elements with
 * all the same value and all the same indexes and non-unique indices are
 * not allowed.  Allowing non-unique indexes would make some operations
 * (e.g. dot product) very difficult to implement and therefore can not
 * be allowed.  A sparse vector where <tt>vec.valueStride() == 0</tt> is one
 * where all of the non-zeros have the value <tt>vec.values[0]</tt>.  If
 * <tt>vec.subNz() == 0</tt> for a sparse vector then it is allowed for
 * <tt>vec.values == NULL</tt> and <tt>vec.indices() == NULL</tt>.
 *
 *	This specification allows a lot of flexibility in determining
 * how the vectors are laid out in memory.  However, allowing vectors to be
 * sparse and unsorted may make many user defined operations
 * considerably harder and expensive to implement.
 *
 * To avoid making mistakes in setting the members of this struct use
 * one of the helper functions <tt>RTOp_sparse_sub_vector_from_dense()</tt>,
 * <tt>RTOp_sparse_sub_vector()</tt> or <tt>RTOp_sub_vector_null()</tt>.
 */
template<class Scalar>
class SparseSubVectorT {
public:
  /** \brief . */
  SparseSubVectorT()
    :globalOffset_(0),subDim_(0),subNz_(0)
    ,values_(NULL),valuesStride_(0),indices_(NULL)
    ,indicesStride_(0),localOffset_(0),isSorted_(0)
    {}
  /** \brief . */
  SparseSubVectorT(
    Teuchos_Index globalOffset, Teuchos_Index subDim, Teuchos_Index subNz
    ,const Scalar values[], ptrdiff_t valuesStride
    ,const Teuchos_Index indices[], ptrdiff_t indicesStride
    ,ptrdiff_t localOffset, int isSorted
    )
    :globalOffset_(globalOffset),subDim_(subDim),subNz_(subNz)
    ,values_(values),valuesStride_(valuesStride),indices_(indices)
    ,indicesStride_(indicesStride),localOffset_(localOffset),isSorted_(isSorted)
    {}
  /** \brief . */
  SparseSubVectorT(
    Teuchos_Index globalOffset, Teuchos_Index subDim
    ,const Scalar values[], ptrdiff_t valuesStride
    )
    :globalOffset_(globalOffset),subDim_(subDim),subNz_(subDim)
    ,values_(values),valuesStride_(valuesStride),indices_(NULL)
    ,indicesStride_(0),localOffset_(0),isSorted_(0)
    {}
  /** \brief . */
  SparseSubVectorT( const ConstSubVectorView<Scalar>& sv )
    :globalOffset_(sv.globalOffset()), subDim_(sv.subDim()), subNz_(sv.subDim()), values_(sv.values())
    ,valuesStride_(sv.stride()), indices_(NULL),indicesStride_(0),localOffset_(0),isSorted_(0)
    {}
  /** \brief . */
  void initialize(
    Teuchos_Index globalOffset, Teuchos_Index subDim, Teuchos_Index subNz
    ,const Scalar values[], ptrdiff_t valuesStride
    ,const Teuchos_Index indices[], ptrdiff_t indicesStride
    ,ptrdiff_t localOffset, int isSorted
    )
    {
      globalOffset_ = globalOffset; subDim_ = subDim; subNz_ = subNz;
      values_ = values; valuesStride_ = valuesStride; indices_ = indices;
      indicesStride_ = indicesStride; localOffset_ = localOffset; isSorted_ = isSorted;
    }
  /** \brief . */
  void initialize(
    Teuchos_Index globalOffset, Teuchos_Index subDim
    ,const Scalar values[], ptrdiff_t valuesStride
    )
    {
      globalOffset_ = globalOffset; subDim_ = subDim; subNz_ = subDim;
      values_ = values; valuesStride_ = valuesStride; indices_ = NULL;
      indicesStride_ = 0; localOffset_ = 0; isSorted_ = 1;
    }
  /** \brief . */
  void set_uninitialized()
    {
      globalOffset_ = 0; subDim_ = 0; subNz_ = 0;
      values_ = NULL; valuesStride_ = 0; indices_ = NULL;
      indicesStride_ = 0; localOffset_ = 0; isSorted_ = 1;
    }
  /** \brief . */
  void setGlobalOffset(Teuchos_Index globalOffset) { globalOffset_ = globalOffset; } 
  /// Offset for the sub-vector into the global vector
  Teuchos_Index                  globalOffset() const { return globalOffset_; } 
  /// Dimension of the sub-vector
  Teuchos_Index                  subDim() const { return subDim_; }
  /// Number of nonzero elements (<tt>subNz == subDim</tt> for dense vectors)
  Teuchos_Index                  subNz() const { return subNz_; }
  /// Array (size min{|<tt>valueStride*subNz</tt>|,1}) for the values in the vector
  const Scalar*                    values() const { return values_; }
  /// Stride between elements in <tt>values[]</tt>
  ptrdiff_t                        valuesStride() const { return valuesStride_; }
  /** \brief Array (size min{|<tt>indicesStride*subNz</tt>|,1} if not <tt>NULL</tt>) for the
    * indices of the nonzero elements in the vector (sparse vectors only)
    */
  const Teuchos_Index*           indices() const { return indices_; }
  /// Stride between indices in indices[] (sparse vectors only)
  ptrdiff_t                        indicesStride() const { return indicesStride_; }
  /// Offset of indices[] into local sub-vector (sparse vectors only)
  ptrdiff_t                        localOffset() const { return localOffset_; }
  /// If <tt>isSorted == 0</tt> then the vector is not sorted, otherwise it is sorted (sparse vectors only)
  int                              isSorted() const { return isSorted_; }
private:
  Teuchos_Index                  globalOffset_;
  Teuchos_Index                  subDim_;
  Teuchos_Index                  subNz_;
  const Scalar                     *values_;
  ptrdiff_t                        valuesStride_;
  const Teuchos_Index            *indices_;
  ptrdiff_t                        indicesStride_;
  ptrdiff_t                        localOffset_;
  int                              isSorted_;
};

template<class Scalar>
void assign_entries( const SubVectorView<Scalar> *msv, const ConstSubVectorView<Scalar> &sv )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(msv==NULL);
  TEST_FOR_EXCEPT(msv->subDim() != sv.subDim());
#endif
  for( int i = 0; i < sv.subDim(); ++i ) {
    (*msv)(i) = sv(i);
  }
}

//
// MultiVectorBase subviews
//

/** \brief Class for a non-mutable sub-multi-vector (submatrix).
 *
 * For a sub-multi-vector <tt>mv</tt>, the corresponding entries
 *	in the global multi-vector <tt>X(j)</tt> (one based) are as follows:
 \verbatim

  X(mv.globalOffset()+k1,mv.colOffset()+k2) = mv(k1,k2), for k1=0...mv.subDim()-1, k2=0...mv.numSubCols()-1
 \endverbatim
 * Unlike vectors, there can only be a unit stride between vector elements
 * in a particular column and there is a Fortran-like leading dimension
 * <tt>mv.leadingDim()</tt> that separates corresponding elements in 
 * each column sub-vector.
 *
 * The raw pointer to the first element, in the first column can be
 * obtained from the function <tt>mv.values()</tt>.
 *
 * Warning! the default copy constructor and assignment operators are
 * allowed which results in only pointer copy, not deep copy!  You
 * have been warned!
 */
template<class Scalar>
class ConstSubMultiVectorView {
public:
  /** \brief . */
  ConstSubMultiVectorView()
    :globalOffset_(0), subDim_(0), colOffset_(0), numSubCols_(0)
    ,values_(NULL), leadingDim_(0)
    {}
  /** \brief . */
  ConstSubMultiVectorView(
    Teuchos_Index globalOffset, Teuchos_Index subDim
    ,Teuchos_Index colOffset, Teuchos_Index numSubCols
    ,const Scalar *values, Teuchos_Index leadingDim
    )
    :globalOffset_(globalOffset), subDim_(subDim)
    ,colOffset_(colOffset), numSubCols_(numSubCols)
    ,values_(values), leadingDim_(leadingDim)
    {}
  /** \brief . */
  ConstSubMultiVectorView( const ConstSubMultiVectorView<Scalar>& smv )
    :globalOffset_(smv.globalOffset()), subDim_(smv.subDim())
    ,colOffset_(smv.colOffset()), numSubCols_(smv.numSubCols())
    ,values_(smv.values()), leadingDim_(smv.leadingDim())
    {}
  /** \brief . */
  void initialize(
    Teuchos_Index globalOffset, Teuchos_Index subDim
    ,Teuchos_Index colOffset, Teuchos_Index numSubCols
    ,const Scalar *values, Teuchos_Index leadingDim
    )
    { globalOffset_=globalOffset; subDim_=subDim; colOffset_=colOffset; numSubCols_=numSubCols;
    values_=values; leadingDim_=leadingDim; }
  /** \brief . */
  void set_uninitialized()
    { globalOffset_ = 0; subDim_=0; colOffset_=0, numSubCols_=0; values_=NULL; leadingDim_=0; }
  /** \brief . */
  void setGlobalOffset(Teuchos_Index globalOffset) { globalOffset_ = globalOffset; } 
  /** \brief . */
  Teuchos_Index   globalOffset()   const { return globalOffset_; }
  /** \brief . */
  Teuchos_Index   subDim()         const { return subDim_; }
  /** \brief . */
  Teuchos_Index   colOffset()      const { return colOffset_; }
  /** \brief . */
  Teuchos_Index   numSubCols()     const { return numSubCols_; }
  /** \brief . */
  const Scalar*     values()         const { return values_; }
  /** \brief . */
  Teuchos_Index   leadingDim()     const { return leadingDim_;  }
  /// Zero-based indexing (Preconditions: <tt>values()!=NULL && (0<=i<subDim()) && (0<=j< numSubCols()</tt>)
  const Scalar& operator()(Teuchos_Index i, Teuchos_Index j) const
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        !( 0 <= i && i < subDim_ ), std::logic_error
        ,"Error, index i="<<i<<" does not fall in the range [0,"<<(subDim_-1)<<"]!"
        );
      TEST_FOR_EXCEPTION(
        !( 0 <= j && j < numSubCols_ ), std::logic_error
        ,"Error, index j="<<j<<" does not fall in the range [0,"<<(numSubCols_-1)<<"]!"
        );
#endif
      return values_[ i + leadingDim_*j ];
    }
  /// Return a <tt>ConstSubVectorView</tt> view of the jth sub-column (Preconditions: <tt>values()!=NULL && (0<=j<numSubCols()</tt>)
  ConstSubVectorView<Scalar> col( const Teuchos_Index j ) const
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        !( 0 <= j && j < numSubCols_ ), std::logic_error
        ,"Error, index j="<<j<<" does not fall in the range [0,"<<(numSubCols_-1)<<"]!"
        );
#endif
      return ConstSubVectorView<Scalar>(globalOffset(),subDim(),values()+j*leadingDim(),1);
    }
private:
  Teuchos_Index     globalOffset_;
  Teuchos_Index     subDim_;
  Teuchos_Index     colOffset_;
  Teuchos_Index     numSubCols_;
  const Scalar        *values_;
  Teuchos_Index     leadingDim_;
};

/** \brief Class for a mutable sub-vector.
 *
 * This class derives from <tt>ConstSubVectorView</tt> and adds methods to
 * mutate the data.  Note, a <tt>const SubVectorView</tt> object
 * allows clients to change the values in the underlying subvector.
 * The meaning of <tt>const</tt> in this context is that the
 * view of the data can not change.
 *
 * Warning! the default copy constructor and assignment operators are
 * allowed which results in only pointer copy, not deep copy!  You
 * have been warned!
 */
template<class Scalar>
class SubMultiVectorView : public ConstSubMultiVectorView<Scalar> {
public:
  /** \brief . */
  SubMultiVectorView() {}
  /** \brief . */
  SubMultiVectorView(
    Teuchos_Index globalOffset, Teuchos_Index subDim
    ,Teuchos_Index colOffset, Teuchos_Index numSubCols
    ,const Scalar *values, Teuchos_Index leadingDim
    )
    :ConstSubMultiVectorView<Scalar>(globalOffset,subDim,colOffset,numSubCols,values,leadingDim)
    {}
  /** \brief . */
  SubMultiVectorView( const SubMultiVectorView<Scalar> & smv)
    :ConstSubMultiVectorView<Scalar>(smv)
    {}
  /** \brief . */
  void initialize(
    Teuchos_Index globalOffset, Teuchos_Index subDim
    ,Teuchos_Index colOffset, Teuchos_Index numSubCols
    ,const Scalar *values, Teuchos_Index leadingDim
    )
    { ConstSubMultiVectorView<Scalar>::initialize(globalOffset,subDim,colOffset,numSubCols,values,leadingDim); }
  /** \brief . */
  void set_uninitialized()
    { ConstSubMultiVectorView<Scalar>::set_uninitialized(); }
  /** \brief . */
  Scalar* values() const { return const_cast<Scalar*>(ConstSubMultiVectorView<Scalar>::values());  }
  /// Zero-based indexing (Preconditions: <tt>values()!=NULL && (0<=i< subDim()) && (0<=j<numSubCols()</tt>)
  Scalar& operator()(Teuchos_Index i, Teuchos_Index j) const
    { return const_cast<Scalar&>(ConstSubMultiVectorView<Scalar>::operator()(i,j)); } // Is range checked in subclass
  /// Return a <tt>SubVectorView</tt> view of the jth sub-column (Preconditions: <tt>values()!=NULL && && (0<=j<numSubCols()</tt>)
  SubVectorView<Scalar> col( const Teuchos_Index j ) const
    {
#ifdef TEUCHOS_DEBUG
      TEST_FOR_EXCEPTION(
        !( 0 <= j && j < this->numSubCols() ), std::logic_error
        ,"Error, index j="<<j<<" does not fall in the range [0,"<<(this->numSubCols()-1)<<"]!"
        );
#endif
      return SubVectorView<Scalar>(this->globalOffset(),this->subDim(),values()+j*this->leadingDim(),1);
    }
};

template<class Scalar>
void assign_entries( const SubMultiVectorView<Scalar> *msmv, const ConstSubMultiVectorView<Scalar> &smv )
{
#ifdef TEUCHOS_DEBUG
  TEST_FOR_EXCEPT(msmv==NULL);
  TEST_FOR_EXCEPT(msmv->subDim() != smv.subDim());
  TEST_FOR_EXCEPT(msmv->numSubCols() != smv.numSubCols());
#endif
  for( Teuchos_Index j = 0; j < smv.numSubCols(); ++j ) {
    for( Teuchos_Index i = 0; i < smv.subDim(); ++i ) {
      (*msmv)(i,j) = smv(i,j);
    }
  }
}

//
// Templated types
//

/** \brief . */
template<class Scalar>  class RTOpT;

} // namespace RTOpPack

#endif // RTOPPACK_TYPES_HPP

