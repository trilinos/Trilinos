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


#ifndef RTOPPACK_SPARSE_SUB_VECTOR_T_HPP
#define RTOPPACK_SPARSE_SUB_VECTOR_T_HPP


#include "RTOpPack_Types.hpp"
#include "Teuchos_Assert.hpp"


namespace RTOpPack {


/** \brief Class for a (sparse or dense) sub-vector.
 *
 * Warning! This class is really nothing more than a dumb data container!
 *
 * Sparse and dense local vectors are supported as follows:
 *
 * A dense vector <tt>vec</tt> is identified by <tt>vec.subDim() ==
 * vec.subNz()</tt> and <tt>vec.indices() == Teuchos::null</tt> in which case
 * <tt>vec.indicesStride()</tt>, <tt>vec.localOffset()</tt> and
 * <tt>vec.isSorted()</tt> are ignored.  For a dense sub-vector <tt>vec</tt>,
 * the corresponding entries in the global vector <tt>x(j)</tt> (zero-based)
 * are as follows:

 \verbatim

 x(vec.globalOffset() + k) == vec.values()[vec.valueStride() * k],

     for k = 0,...,vec.subDim()-1

 \endverbatim

 * The stride member <tt>vec.valueStride()</tt> may be positive (>0) or
 * negative (<0) but not zero (0).  A negative stride <tt>vec.valueStride() <
 * 0</tt> allows a reverse traversal of the elements in <tt>vec.values()</tt>.
 *
 * A sparse vector is identified by <tt>vec.subDim() > vec.subNz()</tt> or
 * <tt>vec.indices() != Teuchos::null</tt> in which case all the fields in the
 * structure are meaningful.  The corresponding elements in the global vector
 * <tt>x(j)</tt> defined as:

 \verbatim

 x( vec.globalOffset() + vec.localOffset() + vec.indices()[vec.indicesStride()*k] )
 == vec.values[vec.valueStride()*k],

     for k = 0,...,vec.subNz()-1

 \endverbatim

 * If <tt>vec.subNz() == 0</tt> then <tt>vec.indices() == Teuchos::null</tt>.
 * If <tt>vec.subDim() > vec.subNz() > 0</tt>, then <tt>vec.indices() !=
 * Teuchos::null</tt> must be true.
 *
 * A sparse sub-vector may be sorted (<tt>vec.isSorted()!=0</tt>) or unsorted
 * (<tt>vec.isSorted()==0</tt>) but the indices <tt>vec.indices()[k]</tt> must
 * be unique.  A sorted vector (<tt>vec.isSorted()!=0</tt>) means that the
 * indices are in ascending order:

 \verbatim

 vec.indices()[vec.indicesStride()*k] < vec.indices()[vec.indicesStride()*k],

     for k = 0,...,vec.subNz()-1

 \endverbatim

 * The member <tt>vec.localOffset()</tt> is used to shift the values in
 * <tt>vec.indices()[]</tt> to be in range of the local sub-vector.  In other
 * words:

 \verbatim
  
 0 <= vec.localOffset() + vec.indices()[vec.indicesStride()*k] < vec.subDim(),

     for k = 0...vec.subNz()-1

 \endverbatim

 * The member <tt>vec.valueStride()</tt> may be positive (>0) or negative (<0)
 * but not zero (0).  Also, the member <tt>vec.indicesStride()</tt> may be may
 * be positive (>0) or negative (<0) but not zero (0).  If <tt>vec.subNz() ==
 * 0</tt> for a sparse vector then <tt>vec.values == Teuchos::null</tt> and
 * <tt>vec.indices() == Teuchos::null</tt>.
 *
 * This specification allows a lot of flexibility in determining how the
 * vectors are laid out in memory.  However, allowing vectors to be sparse and
 * unsorted may make many user defined operations considerably harder and
 * expensive to implement.
 */
template<class Scalar>
class SparseSubVectorT {
public:
  /** \brief . */
  SparseSubVectorT()
    :globalOffset_(0), subDim_(0), subNz_(0),
     valuesStride_(0), indicesStride_(0), localOffset_(0),
     isSorted_(false)
    {}
  /** \brief . */
  SparseSubVectorT(
    Teuchos_Index globalOffset_in, Teuchos_Index subDim_in, Teuchos_Index subNz_in,
    const ArrayRCP<const Scalar> &values_in, ptrdiff_t valuesStride_in,
    const ArrayRCP<const Teuchos_Index> &indices_in, ptrdiff_t indicesStride_in,
    ptrdiff_t localOffset_in, bool isSorted_in
    )
    :globalOffset_(globalOffset_in), subDim_(subDim_in), subNz_(subNz_in),
     values_(values_in), valuesStride_(valuesStride_in), indices_(indices_in),
     indicesStride_(indicesStride_in), localOffset_(localOffset_in), isSorted_(isSorted_in)
    {
#ifdef TEUCHOS_DEBUG
      // Call initialize(...) just to check the preconditions
      initialize(globalOffset_in, subDim_in, subNz_in, values_in, valuesStride_in,
        indices_in, indicesStride_in, localOffset_in, isSorted_in);
#endif
    }
  /** \brief . */
  SparseSubVectorT(
    Teuchos_Index globalOffset_in, Teuchos_Index subDim_in,
    const ArrayRCP<const Scalar> &values_in, ptrdiff_t valuesStride_in
    )
    :globalOffset_(globalOffset_in), subDim_(subDim_in), subNz_(subDim_in),
     values_(values_in), valuesStride_(valuesStride_in),
     indicesStride_(0), localOffset_(0), isSorted_(true)
    {
#ifdef TEUCHOS_DEBUG
      // Call initialize(...) just to check the preconditions
      initialize(globalOffset_in, subDim_in, values_in, valuesStride_in);
#endif
    }
  /** \brief . */
  SparseSubVectorT( const ConstSubVectorView<Scalar>& sv )
    :globalOffset_(sv.globalOffset()), subDim_(sv.subDim()), subNz_(sv.subDim()),
     values_(sv.values()),  valuesStride_(sv.stride()), indicesStride_(0),
     localOffset_(0), isSorted_(true)
    {}
  /** \brief . */
  void initialize(
    Teuchos_Index globalOffset_in, Teuchos_Index subDim_in, Teuchos_Index subNz_in,
    const ArrayRCP<const Scalar> &values_in, ptrdiff_t valuesStride_in,
    const ArrayRCP<const Teuchos_Index> &indices_in, ptrdiff_t indicesStride_in,
    ptrdiff_t localOffset_in, bool isSorted_in
    )
    {
#ifdef TEUCHOS_DEBUG
      TEUCHOS_ASSERT(globalOffset_in >= 0);
      TEUCHOS_ASSERT(subDim_in > 0);
      TEUCHOS_ASSERT_IN_RANGE_UPPER_EXCLUSIVE(subNz_in, 0, subDim_in+1);
      TEUCHOS_ASSERT_EQUALITY(values_in.lowerOffset(), 0);
      TEUCHOS_ASSERT(valuesStride_in != 0);
      TEUCHOS_ASSERT_EQUALITY(values_in.size(), subNz_in*std::abs(valuesStride_in));
      if (!is_null(indices_in)) {
        TEUCHOS_ASSERT(indicesStride_in != 0);
        TEUCHOS_ASSERT_EQUALITY(indices_in.size(), subNz_in*std::abs(indicesStride_in));
        // Note: localOffset can be +, -, or 0 so there is nothing to assert!
        if (isSorted_in) {
          for (int k = 0; k < subNz_in-1; ++k) {
            const Teuchos_Index idx_k = indices_in[k*indicesStride_in];
            const Teuchos_Index idx_kp1 = indices_in[(k+1)*indicesStride_in];
            TEST_FOR_EXCEPTION( !(idx_k < idx_kp1), std::out_of_range,
              "Error indices["<<k<<"]="<<idx_k<<" >= indices["<<k+1<<"]="<<idx_kp1
              <<"!" );
          }
        }
      }
#endif
      globalOffset_ = globalOffset_in;
      subDim_ = subDim_in;
      subNz_ = subNz_in;
      values_ = values_in;
      valuesStride_ = valuesStride_in;
      indices_ = indices_in;
      indicesStride_ = indicesStride_in;
      localOffset_ = localOffset_in;
      isSorted_ = isSorted_in;
    }
  /** \brief . */
  void initialize(
    Teuchos_Index globalOffset_in, Teuchos_Index subDim_in,
    const ArrayRCP<const Scalar> &values_in, ptrdiff_t valuesStride_in
    )
    {
      initialize(globalOffset_in, subDim_in, subDim_in, values_in, valuesStride_in,
        Teuchos::null, 0, 0, true);
    }
  /** \brief . */
  void uninitialize()
    {
      globalOffset_ = 0; subDim_ = 0; subNz_ = 0;
      values_ = Teuchos::null; valuesStride_ = 0; indices_ = Teuchos::null;
      indicesStride_ = 0; localOffset_ = 0; isSorted_ = false;
    }
  /** \brief . */
  void setGlobalOffset(Teuchos_Index globalOffset_in) { globalOffset_ = globalOffset_in; } 
  /** \brief Offset for the sub-vector into the global vector. */
  Teuchos_Index globalOffset() const { return globalOffset_; } 
  /** \brief Dimension of the sub-vector. */
  Teuchos_Index subDim() const { return subDim_; }
  /** \brief Number of nonzero elements (<tt>subNz == subDim</tt> for dense
   * vectors). */
  Teuchos_Index subNz() const { return subNz_; }
  /** \brief Array (size min{|<tt>valueStride*subNz</tt>|,1}) for the values
   * in the vector. */
  const ArrayRCP<const Scalar> values() const { return values_; }
  /** \brief Stride between elements in <tt>values[]</tt>. */
  ptrdiff_t valuesStride() const { return valuesStride_; }
  /** \brief Array (size <tt>indicesStride*subNz</tt>) if not
   * <tt>Teuchos::null</tt>) for the indices of the nonzero elements in the
   * vector (sparse vectors only).
   */
  const ArrayRCP<const Teuchos_Index> indices() const { return indices_; }
  /** \brief Stride between indices in indices[] (sparse vectors only). */
  ptrdiff_t indicesStride() const { return indicesStride_; }
  /** \brief Offset of indices[] into local sub-vector (sparse vectors
   * only). */
  ptrdiff_t localOffset() const { return localOffset_; }
  /** \brief If <tt>isSorted == false</tt> then the vector is not sorted,
   * otherwise it is sorted (sparse vectors only). */
  bool isSorted() const { return isSorted_; }
private:
  Teuchos_Index globalOffset_;
  Teuchos_Index subDim_;
  Teuchos_Index subNz_;
  ArrayRCP<const Scalar> values_;
  ptrdiff_t valuesStride_;
  ArrayRCP<const Teuchos_Index> indices_;
  ptrdiff_t indicesStride_;
  ptrdiff_t localOffset_;
  bool isSorted_;
public:
  /** \brief Deprecated. */
  RTOP_DEPRECATED SparseSubVectorT(
    Teuchos_Index globalOffset_in, Teuchos_Index subDim_in, Teuchos_Index subNz_in,
    const Scalar values_in[], ptrdiff_t valuesStride_in,
    const Teuchos_Index indices_in[], ptrdiff_t indicesStride_in,
    ptrdiff_t localOffset_in, int isSorted_in
    )
    {
      initialize(globalOffset, subDim, subNz, values, valuesStride,
        indices, indicesStride, localOffset, isSorted);
    }
  /** \brief Deprecated. */
  RTOP_DEPRECATED SparseSubVectorT(
    Teuchos_Index globalOffset_in, Teuchos_Index subDim_in,
    const Scalar values_in[], ptrdiff_t valuesStride_in
    )
    {
      initialize(globalOffset_in, subDim_in, values_in, valuesStride_in);
    }
  /** \brief Deprecated. */
  RTOP_DEPRECATED void initialize(
    Teuchos_Index globalOffset_in, Teuchos_Index subDim_in, Teuchos_Index subNz_in,
    const Scalar values_in[], ptrdiff_t valuesStride_in,
    const Teuchos_Index indices_in[], ptrdiff_t indicesStride_in,
    ptrdiff_t localOffset_in, int isSorted_in
    )
    {
      initialize(globalOffset_in, subDim_in, subNz_in,
        Teuchos::arcp(values_in, 0, subNz_in*std::abs(valuesStride_in), false), valuesStride_in,
        Teuchos::arcp(indices_in, 0, subNz_in*std::abs(indicesStride_in), false), indicesStride_in,
        localOffset_in, isSorted_in);
    }
  /** \brief Deprecated. */
  RTOP_DEPRECATED void initialize(
    Teuchos_Index globalOffset_in, Teuchos_Index subDim_in,
    const Scalar values_in[], ptrdiff_t valuesStride_in
    )
    {
      initialize(globalOffset_in, subDim_in,
        Teuchos::arcp(values_in, 0, subDim_in*std::abs(valuesStride_in), false),
        valuesStride_in);
    }
  /** \brief Deprecated. */
  RTOP_DEPRECATED void set_uninitialized()
    { uninitialize(); }
};


} // namespace RTOpPack


#endif // RTOPPACK_SPARSE_SUB_VECTOR_T_HPP
