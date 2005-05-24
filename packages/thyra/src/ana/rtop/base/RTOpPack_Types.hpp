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

#include "RTOp_config.h"
#include "Teuchos_TestForException.hpp"

namespace RTOpPack {

//
// Basic types
//

/// Depreciated
typedef RTOp_value_type value_type;
///Depreciated
typedef RTOp_index_type index_type;
/// Depreciated
typedef RTOp_char_type  char_type;
/** \brief . */
typedef RTOp_index_type Index;

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

  x( vec.globalOffset() + k ) = v.(k), for k = 1,...,vec.subDim()
  \endverbatim
 * The stride <tt>vec.stride()</tt> may be positive (>0), negative (<0)
 * or even zero (0).  A negative stride <tt>vec.stride() < 0</tt> allows a
 * reverse traversal of the elements.  A zero stride
 * <tt>vec.stride()</tt> allows a sub-vector with all the elements the same.
 *
 * The raw pointer to the start of the memory can be obtained as
 * <tt>&vec(1)</tt>.
 *
 * Warning! the default copy constructor and assignment operators are
 * allowed which results in only pointer copy, not deep copy!  You
 * have been warned!
 */
template<class Scalar>
class SubVectorT {
public:
  /** \brief . */
  SubVectorT() : globalOffset_(0), subDim_(0), values_(NULL), stride_(0) {}
  /** \brief . */
  SubVectorT(RTOp_index_type globalOffset, RTOp_index_type subDim, const Scalar *values, ptrdiff_t stride)
    :globalOffset_(globalOffset), subDim_(subDim), values_(values), stride_(stride) 
    {}
  /** \brief . */
  SubVectorT( const SubVectorT<Scalar>& sv )
    :globalOffset_(sv.globalOffset()), subDim_(sv.subDim()), values_(sv.values()), stride_(sv.stride()) 
    {}
  /** \brief . */
  void initialize(RTOp_index_type globalOffset, RTOp_index_type subDim, const Scalar *values, ptrdiff_t stride)
    { globalOffset_=globalOffset; subDim_=subDim; values_=values; stride_=stride;  }
  /** \brief . */
  void set_uninitialized()
    { globalOffset_ = 0; subDim_=0; values_=NULL; stride_ = 0; }
  /** \brief . */
  void setGlobalOffset(RTOp_index_type globalOffset) { globalOffset_ = globalOffset; } 
  /** \brief . */
  RTOp_index_type   globalOffset() const { return globalOffset_; }
  /** \brief . */
  RTOp_index_type   subDim()       const { return subDim_;  }
  /** \brief . */
  const Scalar*     values()       const { return values_;  }
  /** \brief . */
  ptrdiff_t         stride()       const { return stride_;  }
  /// Zero-based indexing (Preconditions: <tt>values()!=NULL && (0 <= i <= subDim()-1)</tt>)
  const Scalar& operator[](RTOp_index_type i) const { return values_[ stride_*i ]; }
  /// One-based indexing (Preconditions: <tt>values()!=NULL && (1 <= i <= subDim())</tt>)
  const Scalar& operator()(RTOp_index_type i) const { return values_[ stride_*(i-1) ]; }
private:
  RTOp_index_type     globalOffset_;
  RTOp_index_type     subDim_;
  const Scalar        *values_;
  ptrdiff_t           stride_;
};

/** \brief Class for a mutable sub-vector.
 *
 * This class derives from <tt>SubVectorT</tt> and adds methods to
 * mutate the data.  Note, a <tt>const MutableSubVectorT</tt> object
 * allows clients to change the values in the underlying subvector.
 * The meaning of <tt>const</tt> in this context is that the
 * view of the data can not change.
 *
 * Warning! the default copy constructor and assignment operators are
 * allowed which results in only pointer copy, not deep copy!  You
 * have been warned!
 */
template<class Scalar>
class MutableSubVectorT : public SubVectorT<Scalar> {
public:
  /** \brief . */
  MutableSubVectorT() {}
  /** \brief . */
  MutableSubVectorT(RTOp_index_type globalOffset, RTOp_index_type subDim, Scalar *values, ptrdiff_t stride)
    :SubVectorT<Scalar>(globalOffset, subDim, values, stride)
    {}
  /** \brief . */
  MutableSubVectorT( const MutableSubVectorT<Scalar> & sv)
    :SubVectorT<Scalar>(sv)
    {}
  /** \brief . */
  void initialize(RTOp_index_type globalOffset, RTOp_index_type subDim, Scalar *values, ptrdiff_t stride)
    { SubVectorT<Scalar>::initialize(globalOffset, subDim, values, stride); }
  /** \brief . */
  void set_uninitialized()
    { SubVectorT<Scalar>::set_uninitialized(); }
  /** \brief . */
  Scalar* values() const { return const_cast<Scalar*>(SubVectorT<Scalar>::values());  }
  /// Zero-based indexing (Preconditions: <tt>values()!=NULL && (0 <= i <= subDim()-1)</tt>)
  Scalar& operator[](RTOp_index_type i) const { return const_cast<Scalar&>(SubVectorT<Scalar>::operator[](i)); }
  /// One-based indexing (Preconditions: <tt>values()!=NULL && (1 <= i <= subDim())</tt>)
  Scalar& operator()(RTOp_index_type i) const { return const_cast<Scalar&>(SubVectorT<Scalar>::operator()(i)); }
};

/** \brief Class for a (sparse or dense) sub-vector.
 *
 *	Sparse and dense local vectors are supported as follows:
  *
  *	A dense vector <tt>vec</tt> is identified by <tt>vec.subDim() == vec.sub_nz</tt>
  * and <tt>vec.indices() == NULL</tt> in which case
  *	<tt>vec.indicesStride()</tt>, <tt>vec.localOffset()</tt> and <tt>vec.isSorted()</tt>
  * are ignored.  For a dense sub-vector <tt>vec</tt>, the corresponding entries
 *	in the global vector <tt>x(j)</tt> (one based) are as follows:
 \verbatim

  x( vec.globalOffset() + k )
    = vec.values()[ vec.valueStride() * (k-1) ]

  for k = 1,...,vec.subDim()
 \endverbatim
 * The stride member <tt>vec.valueStride()()</tt> may be positive (>0), negative (<0)
 * or even zero (0).  A negative stride <tt>vec.valueStride() < 0</tt> allows a
 * reverse traversal of the elements in <tt>vec.values()[]</tt>.  A zero stride
 * <tt>vec.valueStride()() == 0</tt> allows a vector with all the elements the same.
 *
 *	A sparse vector is identified by <tt>vec.subDim() > vec.sub_nz()</tt>
 * or <tt>vec.indices() != NULL</tt>
 * in which case all the fields in the structure are meaningful.
 *	The corresponding elements in the global vector <tt>x(j)</tt>
 * defined as:
 \verbatim

  x( vec.globalOffset() + vec.localOffset() + vec.indices()[vec.indicesStride()*(k-1)] )
    = vec.values[vec.valueStride()*(k-1)]

  for k = 1,...,vec.sub_nz
 \endverbatim
 * If <tt>vec.sub_nz == 0</tt> then it is allowed for <tt>vec.indices() == NULL</tt>.
 * If <tt>vec.subDim() > vec.sub_nz > 0</tt> then <tt>vec.indices() != NULL</tt> must be true.
 *
  * A sparse sub-vector may be sorted (<tt>vec.isSorted()!=0</tt>) or
  * unsorted (<tt>vec.isSorted()==0</tt>) but the indices <tt>vec.indices()[k]</tt>
  * must be unique.  A sorted vector (<tt>vec.isSorted()!=0</tt>) means that
  * the indices are in ascending order:
  \verbatim

  vec.indices()[vec.indicesStride()*(k-1)] < vec.indices()[vec.indicesStride()*(k)]

  for k = 1,...,vec.sub_nz-1
 \endverbatim
 * The member <tt>vec.localOffset()</tt> is used to shift the values in <tt>vec.indices()[]</tt>
 * to be in range of the local sub-vector.  In other words:
 \verbatim
  
  1 <= vec.localOffset() + vec.indices()[vec.indicesStride()*(k-1)] <= vec.sub_nz

  for k = 1...vec.sub_nz
 \endverbatim
 * The member <tt>vec.valueStride()</tt> may be positive (>0), negative (<0) or zero (0).
 * However, the member <tt>vec.indicesStride()</tt> may be may be positive (>0)
 * or negative (<0) but not zero (0).  Allowing <tt>vec.indicesStride() == 0</tt>
 * would mean that a vector would have <tt>vec.sub_nz</tt> nonzero elements with
 * all the same value and all the same indexes and non-unique indices are
 * not allowed.  Allowing non-unique indexes would make some operations
 * (e.g. dot product) very difficult to implement and therefore can not
 * be allowed.  A sparse vector where <tt>vec.valueStride() == 0</tt> is one
 * where all of the non-zeros have the value <tt>vec.values[0]</tt>.  If
 * <tt>vec.sub_nz == 0</tt> for a sparse vector then it is allowed for
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
    RTOp_index_type globalOffset, RTOp_index_type subDim, RTOp_index_type subNz
    ,const Scalar values[], ptrdiff_t valuesStride
    ,const RTOp_index_type indices[], ptrdiff_t indicesStride
    ,ptrdiff_t localOffset, int isSorted
    )
    :globalOffset_(globalOffset),subDim_(subDim),subNz_(subNz)
    ,values_(values),valuesStride_(valuesStride),indices_(indices)
    ,indicesStride_(indicesStride),localOffset_(localOffset),isSorted_(isSorted)
    {}
  /** \brief . */
  SparseSubVectorT(
    RTOp_index_type globalOffset, RTOp_index_type subDim
    ,const Scalar values[], ptrdiff_t valuesStride
    )
    :globalOffset_(globalOffset),subDim_(subDim),subNz_(subDim)
    ,values_(values),valuesStride_(valuesStride),indices_(NULL)
    ,indicesStride_(0),localOffset_(0),isSorted_(0)
    {}
  /** \brief . */
  SparseSubVectorT( const SubVectorT<Scalar>& sv )
    :globalOffset_(sv.globalOffset()), subDim_(sv.subDim()), subNz_(sv.subDim()), values_(sv.values())
    ,valuesStride_(sv.stride()), indices_(NULL),indicesStride_(0),localOffset_(0),isSorted_(0)
    {}
  /** \brief . */
  void initialize(
    RTOp_index_type globalOffset, RTOp_index_type subDim, RTOp_index_type subNz
    ,const Scalar values[], ptrdiff_t valuesStride
    ,const RTOp_index_type indices[], ptrdiff_t indicesStride
    ,ptrdiff_t localOffset, int isSorted
    )
    {
      globalOffset_ = globalOffset; subDim_ = subDim; subNz_ = subNz;
      values_ = values; valuesStride_ = valuesStride; indices_ = indices;
      indicesStride_ = indicesStride; localOffset_ = localOffset; isSorted_ = isSorted;
    }
  /** \brief . */
  void initialize(
    RTOp_index_type globalOffset, RTOp_index_type subDim
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
  void setGlobalOffset(RTOp_index_type globalOffset) { globalOffset_ = globalOffset; } 
  /// Offset for the sub-vector into the global vector
  RTOp_index_type                  globalOffset() const { return globalOffset_; } 
  /// Dimension of the sub-vector
  RTOp_index_type                  subDim() const { return subDim_; }
  /// Number of nonzero elements (<tt>subNz == subDim</tt> for dense vectors)
  RTOp_index_type                  subNz() const { return subNz_; }
  /// Array (size min{|<tt>valueStride*subNz</tt>|,1}) for the values in the vector
  const Scalar*                    values() const { return values_; }
  /// Stride between elements in <tt>values[]</tt>
  ptrdiff_t                        valuesStride() const { return valuesStride_; }
  /** \brief Array (size min{|<tt>indicesStride*subNz</tt>|,1} if not <tt>NULL</tt>) for the
    * indices of the nonzero elements in the vector (sparse vectors only)
    */
  const RTOp_index_type*           indices() const { return indices_; }
  /// Stride between indices in indices[] (sparse vectors only)
  ptrdiff_t                        indicesStride() const { return indicesStride_; }
  /// Offset of indices[] into local sub-vector (sparse vectors only)
  ptrdiff_t                        localOffset() const { return localOffset_; }
  /// If <tt>isSorted == 0</tt> then the vector is not sorted, otherwise it is sorted (sparse vectors only)
  int                              isSorted() const { return isSorted_; }
private:
  RTOp_index_type                  globalOffset_;
  RTOp_index_type                  subDim_;
  RTOp_index_type                  subNz_;
  const Scalar                     *values_;
  ptrdiff_t                        valuesStride_;
  const RTOp_index_type            *indices_;
  ptrdiff_t                        indicesStride_;
  ptrdiff_t                        localOffset_;
  int                              isSorted_;
};


template<class Scalar>
void assign_entries( const MutableSubVectorT<Scalar> *msv, const SubVectorT<Scalar> &sv )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(msv==NULL);
  TEST_FOR_EXCEPT(msv->subDim() != sv.subDim());
#endif
  for( int i = 1; i <= sv.subDim(); ++i ) {
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

  X(mv.globalOffset()+k1,mv.colOffset()+k2) = mv(k1,k2), for k1=1...mv.subDim(), k2=1...mv.numSubCols()
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
class SubMultiVectorT {
public:
  /** \brief . */
  SubMultiVectorT()
    :globalOffset_(0), subDim_(0), colOffset_(0), numSubCols_(0)
    ,values_(NULL), leadingDim_(0)
    {}
  /** \brief . */
  SubMultiVectorT(
    RTOp_index_type globalOffset, RTOp_index_type subDim
    ,RTOp_index_type colOffset, RTOp_index_type numSubCols
    ,const Scalar *values, RTOp_index_type leadingDim
    )
    :globalOffset_(globalOffset), subDim_(subDim)
    ,colOffset_(colOffset), numSubCols_(numSubCols)
    ,values_(values), leadingDim_(leadingDim)
    {}
  /** \brief . */
  SubMultiVectorT( const SubMultiVectorT<Scalar>& smv )
    :globalOffset_(smv.globalOffset()), subDim_(smv.subDim())
    ,colOffset_(smv.colOffset()), numSubCols_(smv.numSubCols())
    ,values_(smv.values()), leadingDim_(smv.leadingDim())
    {}
  /** \brief . */
  void initialize(
    RTOp_index_type globalOffset, RTOp_index_type subDim
    ,RTOp_index_type colOffset, RTOp_index_type numSubCols
    ,const Scalar *values, RTOp_index_type leadingDim
    )
    { globalOffset_=globalOffset; subDim_=subDim; colOffset_=colOffset; numSubCols_=numSubCols;
    values_=values; leadingDim_=leadingDim; }
  /** \brief . */
  void set_uninitialized()
    { globalOffset_ = 0; subDim_=0; colOffset_=0, numSubCols_=0; values_=NULL; leadingDim_=0; }
  /** \brief . */
  void setGlobalOffset(RTOp_index_type globalOffset) { globalOffset_ = globalOffset; } 
  /** \brief . */
  RTOp_index_type   globalOffset()   const { return globalOffset_; }
  /** \brief . */
  RTOp_index_type   subDim()         const { return subDim_; }
  /** \brief . */
  RTOp_index_type   colOffset()      const { return colOffset_; }
  /** \brief . */
  RTOp_index_type   numSubCols()     const { return numSubCols_; }
  /** \brief . */
  const Scalar*     values()         const { return values_; }
  /** \brief . */
  RTOp_index_type   leadingDim()     const { return leadingDim_;  }
  /// One-based indexing (Preconditions: <tt>values()!=NULL && (1<=i<=subDim()) && (1<=j<= numSubCols()</tt>)
  const Scalar& operator()(RTOp_index_type i, RTOp_index_type j) const
    { return values_[ (i-1) + leadingDim_*(j-1) ]; }
  /// Return a <tt>SubVectorT</tt> view of the jth sub-column (Preconditions: <tt>values()!=NULL && (1<=i<=subDim()) && (1<=j<= numSubCols()</tt>)
  SubVectorT<Scalar> col( const RTOp_index_type j ) const
    { return SubVectorT<Scalar>(globalOffset(),subDim(),values()+(j-1)*leadingDim(),1); }
private:
  RTOp_index_type     globalOffset_;
  RTOp_index_type     subDim_;
  RTOp_index_type     colOffset_;
  RTOp_index_type     numSubCols_;
  const Scalar        *values_;
  RTOp_index_type     leadingDim_;
};

/** \brief Class for a mutable sub-vector.
 *
 * This class derives from <tt>SubVectorT</tt> and adds methods to
 * mutate the data.  Note, a <tt>const MutableSubVectorT</tt> object
 * allows clients to change the values in the underlying subvector.
 * The meaning of <tt>const</tt> in this context is that the
 * view of the data can not change.
 *
 * Warning! the default copy constructor and assignment operators are
 * allowed which results in only pointer copy, not deep copy!  You
 * have been warned!
 */
template<class Scalar>
class MutableSubMultiVectorT : public SubMultiVectorT<Scalar> {
public:
  /** \brief . */
  MutableSubMultiVectorT() {}
  /** \brief . */
  MutableSubMultiVectorT(
    RTOp_index_type globalOffset, RTOp_index_type subDim
    ,RTOp_index_type colOffset, RTOp_index_type numSubCols
    ,const Scalar *values, RTOp_index_type leadingDim
    )
    :SubMultiVectorT<Scalar>(globalOffset,subDim,colOffset,numSubCols,values,leadingDim)
    {}
  /** \brief . */
  MutableSubMultiVectorT( const MutableSubMultiVectorT<Scalar> & smv)
    :SubMultiVectorT<Scalar>(smv)
    {}
  /** \brief . */
  void initialize(
    RTOp_index_type globalOffset, RTOp_index_type subDim
    ,RTOp_index_type colOffset, RTOp_index_type numSubCols
    ,const Scalar *values, RTOp_index_type leadingDim
    )
    { SubMultiVectorT<Scalar>::initialize(globalOffset,subDim,colOffset,numSubCols,values,leadingDim); }
  /** \brief . */
  void set_uninitialized()
    { SubMultiVectorT<Scalar>::set_uninitialized(); }
  /** \brief . */
  Scalar* values() const { return const_cast<Scalar*>(SubMultiVectorT<Scalar>::values());  }
  /// One-based indexing (Preconditions: <tt>values()!=NULL && (1 <= i <= subDim()) && (1<= j <= numSubCols()</tt>)
  Scalar& operator()(RTOp_index_type i, RTOp_index_type j) const
    { return const_cast<Scalar&>(SubMultiVectorT<Scalar>::operator()(i,j)); }
  /// Return a <tt>MutableSubVectorT</tt> view of the jth sub-column (Preconditions: <tt>values()!=NULL && (1<=i<=subDim()) && (1<=j<= numSubCols()</tt>)
  MutableSubVectorT<Scalar> col( const RTOp_index_type j ) const
    { return MutableSubVectorT<Scalar>(SubMultiVectorT<Scalar>::globalOffset(),SubMultiVectorT<Scalar>::subDim(),values()+(j-1)*SubMultiVectorT<Scalar>::leadingDim(),1); }
};

template<class Scalar>
void assign_entries( const MutableSubMultiVectorT<Scalar> *msmv, const SubMultiVectorT<Scalar> &smv )
{
#ifdef _DEBUG
  TEST_FOR_EXCEPT(msmv==NULL);
  TEST_FOR_EXCEPT(msmv->subDim() != smv.subDim());
  TEST_FOR_EXCEPT(msmv->numSubCols() != smv.numSubCols());
#endif
  for( RTOp_index_type j = 1; j <= smv.numSubCols(); ++j ) {
    for( RTOp_index_type i = 1; i <= smv.subDim(); ++i ) {
      (*msmv)(i,j) = smv(i,j);
    }
  }
}

//
// Templated types
//

/** \brief . */
template<class Scalar>  class RTOpT;

//
// Typedefs
//

/** \brief . */
typedef SubVectorT<RTOp_value_type>              SubVector;
/** \brief . */
typedef MutableSubVectorT<RTOp_value_type>       MutableSubVector;
/** \brief . */
typedef SparseSubVectorT<RTOp_value_type>        SparseSubVector;
/** \brief . */
typedef SubMultiVectorT<RTOp_value_type>         SubMultiVector;
/** \brief . */
typedef MutableSubMultiVectorT<RTOp_value_type>  MutableSubMultiVector;
/** \brief . */
typedef RTOpT<RTOp_value_type>                   RTOp;

} // namespace RTOpPack

#endif // RTOPPACK_TYPES_HPP

