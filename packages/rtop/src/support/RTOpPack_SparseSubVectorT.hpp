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


#ifndef RTOPPACK_SPARSE_SUB_VECTOR_T_HPP
#define RTOPPACK_SPARSE_SUB_VECTOR_T_HPP


#include "RTOpPack_Types.hpp"


namespace RTOpPack {


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
  /** \brief Offset for the sub-vector into the global vector. */
  Teuchos_Index globalOffset() const { return globalOffset_; } 
  /** \brief Dimension of the sub-vector. */
  Teuchos_Index subDim() const { return subDim_; }
  /** \brief Number of nonzero elements (<tt>subNz == subDim</tt> for dense
   * vectors). */
  Teuchos_Index subNz() const { return subNz_; }
  /** \brief Array (size min{|<tt>valueStride*subNz</tt>|,1}) for the values
   * in the vector. */
  const Scalar* values() const { return values_; }
  /** \brief Stride between elements in <tt>values[]</tt>. */
  ptrdiff_t valuesStride() const { return valuesStride_; }
  /** \brief Array (size min{|<tt>indicesStride*subNz</tt>|,1} if not <tt>NULL</tt>) for the
   * indices of the nonzero elements in the vector (sparse vectors only)
   */
  const Teuchos_Index* indices() const { return indices_; }
  /** \brief Stride between indices in indices[] (sparse vectors only). */
  ptrdiff_t indicesStride() const { return indicesStride_; }
  /** \brief Offset of indices[] into local sub-vector (sparse vectors
   * only). */
  ptrdiff_t localOffset() const { return localOffset_; }
  /** \brief If <tt>isSorted == 0</tt> then the vector is not sorted,
   * otherwise it is sorted (sparse vectors only). */
  int isSorted() const { return isSorted_; }
private:
  Teuchos_Index globalOffset_;
  Teuchos_Index subDim_;
  Teuchos_Index subNz_;
  const Scalar *values_;
  // 2007/11/14: rabartl: ToDo: Replace with ArrayRCP<Scalar>
  ptrdiff_t valuesStride_;
  const Teuchos_Index *indices_;
  ptrdiff_t indicesStride_;
  ptrdiff_t localOffset_;
  int isSorted_;
};


} // namespace RTOpPack


#endif // RTOPPACK_SPARSE_SUB_VECTOR_T_HPP
