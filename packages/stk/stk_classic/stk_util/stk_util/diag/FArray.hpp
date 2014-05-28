/**
 *    Copyright 2005-2009 Sandia Corporation.
 *    Under the terms of Contract DE-AC04-94AL85000, there is a
 *    non-exclusive license for use of this work by or on behalf
 *    of the U.S. Government.  Export of this program may require
 *    a license from the United States Government.
 */

/**
 * @file
 * @author H. Carter Edwards
 * @date   August 2005
 */

#ifndef STK_UTIL_DIAG_FArray_h
#define STK_UTIL_DIAG_FArray_h

#include <cstddef>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <typeinfo>
#include <iterator>

#include <stk_util/util/StaticAssert.hpp>

#ifndef NDEBUG
#  define SIERRA_ARRAY_BOUNDS_CHECK
#endif

namespace sierra {

///
/// @addtogroup FArrayDetail
/// @{
///


void
array_dimension_error(
  const std::type_info &	typeinfo,
  unsigned			dimension,
  unsigned			value,
  unsigned			upper);


/**
 * @class FArray
 *
 * @brief Multidimensional array of contiguous memory.  The memory is not owned by the
 * array, but container access semantics are enforced, i.e. const Array<> elements cannot
 * be assigned to.
 *
 * @par Indexing via 'operator()'
 *
 * Elements of an array can be accessed via the () operator, where an index is provided
 * for each dimension of the array.  The number of dimensions is enforced at compile
 * time; however, enforcement of the range of the indices is only performed at runtime in
 * debug mode, i.e. 'NDEBUG' is not defined.
 *
 *    FArray<double,4> a(mem, n0, n1, n2, n3);
 *    a(i0, i1, i2, i3) = 24 ;
 *
 * @par Construction with a block of memory
 *
 * An array is created by providing it with memory address and dimensions by which that
 * memory is to be viewed.
 *
 *    FArray(pointer, ndim0, ndim1, ndim2, ...);
 *
 * @par Shallow copy construction
 *
 * The copy constructor is a shallow copy.  The constructed array is just another a view
 * into the same memory as the input array.  This is the preferred method of passing
 * arrays as access required less memory indirection.
 */
template<class ElementType, int Dimension>
class FArray;

/**
 * @class sierra::FArrayContainer
 *
 * @brief Extend FArray with deep copy assignment and resize operations.
 */
template<class ElementType, int Dimension, class A = std::allocator<ElementType> >
class FArrayContainer;

//----------------------------------------------------------------------
//----------------------------------------------------------------------
// Helpers

template<unsigned N> struct ArrayHelper;

//----------------------------------------------------------------------
//----------------------------------------------------------------------

template<>
struct ArrayHelper<0>
{
  template<typename T>
  inline static void fill(const T, T *)
  {}

  template<typename T>
  inline static void copy(const T *, T *)
  {}

  template<typename T>
  inline static void prefix(const T *, T *)
  {}

  inline static void get_index(const unsigned * const, unsigned *index, int offset) {
    index[0] = offset;
  }

  template<typename T>
  inline static T *index(T * ptr, const unsigned * const, const unsigned * const, const unsigned * const index) {
    return ptr + index[0];
  }
};

template<>
struct ArrayHelper<1>
{
  template<typename T>
  inline static bool equal(const T * x, const T * y) {
    return *x == *y;
  }

  template<typename T>
  inline static void fill(const T x, T * y) {
    *y = x;
  }

  template<typename T>
  inline static void copy(const T * x, T * y) {
    *y = *x;
  }

  template<typename T>
  inline static void prefix(const T * x, T * p) {
    *(p+1) = *p * *x;
  }

  template<typename T>
  static void copy(const unsigned * const dim,
		   T       * a_ptr, const unsigned * const a_inc,
		   T const * b_ptr, const unsigned * const b_inc)
  {
    const unsigned ia = *a_inc;
    const unsigned ib = *b_inc;
    T * const a_end = a_ptr + *dim * ia;
    while (a_end != a_ptr) {
      *a_ptr = *b_ptr;
      a_ptr += ia;
      b_ptr += ib;
    }
  }

  template<typename T>
  static void fill(const unsigned * const dim,
		   T       * a_ptr, const unsigned * const a_inc,
		   const T & value)
  {
    const unsigned ia = *a_inc;
    T * const a_end = a_ptr + *dim * ia;
    while (a_end != a_ptr) {
      *a_ptr = value;
      a_ptr += ia;
    }
  }

  template<typename T>
  static bool equal(const unsigned * const dim,
		    const T * a_ptr, const unsigned * const a_inc,
		    const T * b_ptr, const unsigned * const b_inc)
  {
    const unsigned ia = *a_inc;
    const unsigned ib = *b_inc;
    const T * const a_end = a_ptr + *dim * ia;
    bool result = true;
    while (a_end != a_ptr && (result = *a_ptr == *b_ptr)) {
      a_ptr += ia;
      b_ptr += ib;
    }
    return result;
  }

  inline static void get_index(const unsigned * const inc, unsigned *index, int offset) {
    index[1] = offset/inc[1];
    index[0] = offset%inc[1];
  }

  template<typename T>
  inline static T *index(T * ptr, const unsigned * const inc, const unsigned * const, const unsigned * const index) {
    return ptr + index[0] + index[1]*inc[1];
  }
};


template<unsigned N>
struct ArrayHelper
{
private:
  typedef ArrayHelper<N-1> M;

public:
  template<typename T>
  inline static bool equal(const T * x, const T * y) {
    return *x == *y && M::equal(x + 1, y + 1);
  }

  template<typename T>
  inline static void fill(const T x, T * y) {
    *y = x;
    M::fill(x, y+1);
  }

  template<typename T>
  inline static void copy(const T * x, T * y) {
    *y = *x;
    M::copy(x+1, y+1);
  }

  template<typename T>
  inline static void prefix(const T * x, T * p) {
    *(p+1) = *p * *x;
    M::prefix(x + 1, p + 1);
  }

  template<typename T>
  static void copy(const unsigned * const dim,
		   T       * a_ptr, const unsigned * const a_inc,
		   T const * b_ptr, const unsigned * const b_inc)
  {
    const unsigned ia = *a_inc;
    const unsigned ib = *b_inc;
    T * const a_end = a_ptr + *dim * ia;
    while (a_end != a_ptr) {
      M::copy(dim + 1, a_ptr, a_inc + 1, b_ptr, b_inc + 1);
      a_ptr += ia;
      b_ptr += ib;
    }
  }

  template<typename T>
  static void fill(const unsigned * const dim,
		   T       * a_ptr, const unsigned * const a_inc,
		   const T & value)
  {
    const unsigned ia = *a_inc;
    T * const a_end = a_ptr + *dim * ia;
    while (a_end != a_ptr) {
      M::fill(dim + 1, a_ptr, a_inc + 1, value);
      a_ptr += ia;
    }
  }

  template<typename T>
  static bool equal(const unsigned * const dim,
		    const T * a_ptr, const unsigned * const a_inc,
		    const T * b_ptr, const unsigned * const b_inc)
  {
    const unsigned ia = *a_inc;
    const unsigned ib = *b_inc;
    const T * const a_end = a_ptr + *dim * ia;
    bool result = true;
    while (a_end != a_ptr &&
	   (result = M::equal(dim+1, a_ptr, a_inc+1, b_ptr, b_inc+1))){
      a_ptr += ia;
      b_ptr += ib;
    }
    return result;
  }

  inline static void get_index(const unsigned * const inc, unsigned *index, int offset) {
    index[N] = offset/inc[N];
    M::get_index(inc, index, offset%inc[N]);
  }

  template<typename T>
  inline static T *index(T * ptr, const unsigned * const inc, const unsigned * const dim, const unsigned * const index) {
    return M::index(ptr, inc, dim, index) + index[N]*inc[N];
  }
};

//----------------------------------------------------------------------
// A typeless array is invalid...

template<int Dimension>
class FArray<void, Dimension> {};

template <class ElementType>
class FArray<ElementType, 0>
{
public:
  typedef ElementType value_type;

  typedef FArray<value_type, 0> SelfType;

  typedef FArray<const value_type, 0> Const;

  enum { NumDim = 0};
};

class FArrayBootstrap
{
public:
  ~FArrayBootstrap();
};

template<class ElementType, int Dimension>
class FArray : public FArrayBootstrap
{
public:
  enum { NumDim = Dimension };

  typedef ElementType           value_type;
  typedef size_t                size_type;
  typedef ptrdiff_t             difference_type;

  typedef value_type *          pointer;
  typedef const value_type *    const_pointer;
  typedef value_type &          reference;
  typedef const value_type &    const_reference;

  typedef pointer               iterator;
  typedef const_pointer         const_iterator;
  typedef typename std::reverse_iterator<iterator> reverse_iterator;
  typedef typename std::reverse_iterator<const_iterator> const_reverse_iterator;

  typedef FArray<ElementType, Dimension> SelfType;
  typedef FArrayContainer<ElementType, Dimension> Container;
  typedef FArray<ElementType, Dimension - 1> Trunc;

  class Index
  {
  public:
    const unsigned &operator[](unsigned i) const {
      return m_index[i];
    }

    unsigned &operator[](unsigned i) {
      return m_index[i];
    }

  private:
    unsigned    m_index[NumDim];
  };

  const unsigned * dimension() const {
    return m_dim;
  }

  const unsigned * stride() const {
    return m_stride;
  }

#ifdef SIERRA_ARRAY_BOUNDS_CHECK
  inline void array_dimension_verify(unsigned l_dimension, unsigned value, unsigned upper) const {
    if (value >= upper)
    {
      array_dimension_error(typeid(*this), l_dimension, value, upper);
    }
  }
#else
  inline void array_dimension_verify(unsigned, unsigned, unsigned) const {}
#endif

  template<unsigned I>
  unsigned dimension() const {
    enum { ok = stk::StaticAssert<I < NumDim>::OK };
    return m_dim[I];
  }

  template<unsigned I>
  unsigned stride() const {
    enum { ok = stk::StaticAssert<I < NumDim>::OK };
    return m_stride[I];
  }

  /** Dimensions defined at construction */
  unsigned dimension(const unsigned i) const {
    array_dimension_verify(0, i, NumDim);
    return m_dim[i];
  }

  unsigned stride(const unsigned i) const {
    array_dimension_verify(0, i, NumDim);
    return m_stride[i];
  }

  unsigned size() const {
    return m_stride[NumDim];
  }

  //----------------------------------------

  bool operator==(const SelfType & a) const {
    return ArrayHelper<NumDim>::equal(m_dim, a.m_dim) &&
      ArrayHelper<NumDim>::equal(m_dim, m_ptr, m_stride, a.m_ptr, a.m_stride);
  }

  template<typename T>
  bool operator==(const FArray<T, Dimension> & a) const {
    return ArrayHelper<NumDim>::equal(m_dim, a.dimension()) &&
      ArrayHelper<NumDim>::equal(m_dim, m_ptr, m_stride, a.ptr(), a.stride());
  }

  bool operator!=(const SelfType & a) const {
    return ! operator==(a); }

  template<typename T>
  bool operator!=(const FArray<T, Dimension> & a) const {
    return ! operator==(a);
  }

  value_type &operator()(const Index &index) {
    for (unsigned i = 0; i < NumDim; ++i)
      array_dimension_verify(i, index[i], m_dim[i]);

    value_type *l_ptr = ArrayHelper<NumDim - 1>::index(m_ptr, m_stride, m_dim, &index[0]);

    return *l_ptr;
  }

  value_type & operator()(const unsigned i0, const unsigned i1,
			    const unsigned i2, const unsigned i3,
			    const unsigned i4, const unsigned i5,
			    const unsigned i6, const unsigned i7)
  {
    enum { ok = stk::StaticAssert<8 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);
    array_dimension_verify(2, i2, m_dim[2]);
    array_dimension_verify(3, i3, m_dim[3]);
    array_dimension_verify(4, i4, m_dim[4]);
    array_dimension_verify(5, i5, m_dim[5]);
    array_dimension_verify(6, i6, m_dim[6]);
    array_dimension_verify(7, i7, m_dim[7]);

    return *(m_ptr +    i0 + m_stride[1] * i1 +
	     m_stride[2] * i2 + m_stride[3] * i3 +
	     m_stride[4] * i4 + m_stride[5] * i5 +
	     m_stride[6] * i6 + m_stride[7] * i7);
  }

  value_type & operator()(const unsigned i0, const unsigned i1,
			    const unsigned i2, const unsigned i3,
			    const unsigned i4, const unsigned i5,
			    const unsigned i6)
  {
    enum { ok = stk::StaticAssert<7 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);
    array_dimension_verify(2, i2, m_dim[2]);
    array_dimension_verify(3, i3, m_dim[3]);
    array_dimension_verify(4, i4, m_dim[4]);
    array_dimension_verify(5, i5, m_dim[5]);
    array_dimension_verify(6, i6, m_dim[6]);

    return *(m_ptr +    i0 + m_stride[1] * i1 +
	     m_stride[2] * i2 + m_stride[3] * i3 +
	     m_stride[4] * i4 + m_stride[5] * i5 +
	     m_stride[6] * i6);
  }

  value_type & operator()(const unsigned i0, const unsigned i1,
			    const unsigned i2, const unsigned i3,
			    const unsigned i4, const unsigned i5)
  {
    enum { ok = stk::StaticAssert<6 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);
    array_dimension_verify(2, i2, m_dim[2]);
    array_dimension_verify(3, i3, m_dim[3]);
    array_dimension_verify(4, i4, m_dim[4]);
    array_dimension_verify(5, i5, m_dim[5]);

    return *(m_ptr +    i0 + m_stride[1] * i1 +
	     m_stride[2] * i2 + m_stride[3] * i3 +
	     m_stride[4] * i4 + m_stride[5] * i5);
  }

  value_type & operator()(const unsigned i0, const unsigned i1,
			    const unsigned i2, const unsigned i3,
			    const unsigned i4)
  {
    enum { ok = stk::StaticAssert<5 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);
    array_dimension_verify(2, i2, m_dim[2]);
    array_dimension_verify(3, i3, m_dim[3]);
    array_dimension_verify(4, i4, m_dim[4]);

    return *(m_ptr +    i0 + m_stride[1] * i1 +
	     m_stride[2] * i2 + m_stride[3] * i3 +
	     m_stride[4] * i4);
  }

  value_type & operator()(const unsigned i0, const unsigned i1,
			    const unsigned i2, const unsigned i3)
  {
    enum { ok = stk::StaticAssert<4 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);
    array_dimension_verify(2, i2, m_dim[2]);
    array_dimension_verify(3, i3, m_dim[3]);

    return *(m_ptr +    i0 + m_stride[1] * i1 +
	     m_stride[2] * i2 + m_stride[3] * i3);
  }

  value_type & operator()(const unsigned i0, const unsigned i1,
			    const unsigned i2)
  {
    enum { ok = stk::StaticAssert<3 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);
    array_dimension_verify(2, i2, m_dim[2]);

    return *(m_ptr +    i0 + m_stride[1] * i1 +
	     m_stride[2] * i2);
  }

  value_type & operator()(const unsigned i0, const unsigned i1) {
    enum { ok = stk::StaticAssert<2 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);

    return *(m_ptr +    i0 + m_stride[1] * i1);
  }

  value_type & operator()(const unsigned i0) {
    enum { ok = stk::StaticAssert<1 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);

    return *(m_ptr +    i0);
  }

  value_type * ptr() {
    return m_ptr;
  }

  value_type & operator[](unsigned i) {
    array_dimension_verify(0, i, m_stride[NumDim]);
    return m_ptr[i];
  }

  const value_type &operator()(const Index &index) const {
    for (unsigned i = 0; i < NumDim; ++i)
      array_dimension_verify(i, index[i], m_dim[i]);

    const value_type *l_ptr = ArrayHelper<NumDim - 1>::index(m_ptr, m_stride, m_dim, &index[0]);

    return *l_ptr;
  }

  const value_type & operator()(const unsigned i0, const unsigned i1,
				  const unsigned i2, const unsigned i3,
				  const unsigned i4, const unsigned i5,
				  const unsigned i6, const unsigned i7) const
  {
    enum { ok = stk::StaticAssert<8 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);
    array_dimension_verify(2, i2, m_dim[2]);
    array_dimension_verify(3, i3, m_dim[3]);
    array_dimension_verify(4, i4, m_dim[4]);
    array_dimension_verify(5, i5, m_dim[5]);
    array_dimension_verify(6, i6, m_dim[6]);
    array_dimension_verify(7, i7, m_dim[7]);

    return *(m_ptr +    i0 + m_stride[1] * i1 +
	     m_stride[2] * i2 + m_stride[3] * i3 +
	     m_stride[4] * i4 + m_stride[5] * i5 +
	     m_stride[6] * i6 + m_stride[7] * i7);
  }

  const value_type & operator()(const unsigned i0, const unsigned i1,
				  const unsigned i2, const unsigned i3,
				  const unsigned i4, const unsigned i5,
				  const unsigned i6) const
  {
    enum { ok = stk::StaticAssert<7 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);
    array_dimension_verify(2, i2, m_dim[2]);
    array_dimension_verify(3, i3, m_dim[3]);
    array_dimension_verify(4, i4, m_dim[4]);
    array_dimension_verify(5, i5, m_dim[5]);
    array_dimension_verify(6, i6, m_dim[6]);

    return *(m_ptr +    i0 + m_stride[1] * i1 +
	     m_stride[2] * i2 + m_stride[3] * i3 +
	     m_stride[4] * i4 + m_stride[5] * i5 +
	     m_stride[6] * i6);
  }

  const value_type & operator()(const unsigned i0, const unsigned i1,
				  const unsigned i2, const unsigned i3,
				  const unsigned i4, const unsigned i5) const
  {
    enum { ok = stk::StaticAssert<6 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);
    array_dimension_verify(2, i2, m_dim[2]);
    array_dimension_verify(3, i3, m_dim[3]);
    array_dimension_verify(4, i4, m_dim[4]);
    array_dimension_verify(5, i5, m_dim[5]);

    return *(m_ptr +    i0 + m_stride[1] * i1 +
	     m_stride[2] * i2 + m_stride[3] * i3 +
	     m_stride[4] * i4 + m_stride[5] * i5);
  }

  const value_type & operator()(const unsigned i0, const unsigned i1,
				  const unsigned i2, const unsigned i3,
				  const unsigned i4) const
  {
    enum { ok = stk::StaticAssert<5 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);
    array_dimension_verify(2, i2, m_dim[2]);
    array_dimension_verify(3, i3, m_dim[3]);
    array_dimension_verify(4, i4, m_dim[4]);

    return *(m_ptr +    i0 + m_stride[1] * i1 +
	     m_stride[2] * i2 + m_stride[3] * i3 +
	     m_stride[4] * i4);
  }

  const value_type & operator()(const unsigned i0, const unsigned i1,
				  const unsigned i2, const unsigned i3) const
  {
    enum { ok = stk::StaticAssert<4 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);
    array_dimension_verify(2, i2, m_dim[2]);
    array_dimension_verify(3, i3, m_dim[3]);

    return *(m_ptr +    i0 + m_stride[1] * i1 +
	     m_stride[2] * i2 + m_stride[3] * i3);
  }

  const value_type & operator()(const unsigned i0, const unsigned i1,
				  const unsigned i2) const
  {
    enum { ok = stk::StaticAssert<3 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);
    array_dimension_verify(2, i2, m_dim[2]);

    return *(m_ptr +    i0 + m_stride[1] * i1 +
	     m_stride[2] * i2);
  }

  const value_type & operator()(const unsigned i0, const unsigned i1) const {
    enum { ok = stk::StaticAssert<2 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);
    array_dimension_verify(1, i1, m_dim[1]);

    return *(m_ptr +    i0 + m_stride[1] * i1);
  }

  const value_type & operator()(const unsigned i0) const {
    enum { ok = stk::StaticAssert<1 == NumDim>::OK };

    array_dimension_verify(0, i0, m_dim[0]);

    return *(m_ptr +    i0);
  }

  value_type * ptr() const {
    return m_ptr;
  }

  const value_type & operator[](unsigned i) const {
    array_dimension_verify(0, i, m_stride[NumDim]);
    return m_ptr[i];
  }

  //----------------------------------------

  bool verify_dimension(const unsigned n0, const unsigned n1,
			const unsigned n2, const unsigned n3,
			const unsigned n4, const unsigned n5,
			const unsigned n6, const unsigned n7) const
  {
    enum { ok = stk::StaticAssert<8 == NumDim>::OK };
    return m_dim[0] == n0 && m_dim[1] == n1 &&
      m_dim[2] == n2 && m_dim[3] == n3 &&
      m_dim[4] == n4 && m_dim[5] == n5 &&
      m_dim[6] == n6 && m_dim[7] == n7;
  }

  bool verify_dimension(const unsigned n0, const unsigned n1,
			const unsigned n2, const unsigned n3,
			const unsigned n4, const unsigned n5,
			const unsigned n6) const
  {
    enum { ok = stk::StaticAssert<7 == NumDim>::OK };
    return m_dim[0] == n0 && m_dim[1] == n1 &&
      m_dim[2] == n2 && m_dim[3] == n3 &&
      m_dim[4] == n4 && m_dim[5] == n5 &&
      m_dim[6] == n6;
  }

  bool verify_dimension(const unsigned n0, const unsigned n1,
			const unsigned n2, const unsigned n3,
			const unsigned n4, const unsigned n5) const
  {
    enum { ok = stk::StaticAssert<6 == NumDim>::OK };
    return m_dim[0] == n0 && m_dim[1] == n1 &&
      m_dim[2] == n2 && m_dim[3] == n3 &&
      m_dim[4] == n4 && m_dim[5] == n5;
  }

  bool verify_dimension(const unsigned n0, const unsigned n1,
			const unsigned n2, const unsigned n3,
			const unsigned n4) const
  {
    enum { ok = stk::StaticAssert<5 == NumDim>::OK };
    return m_dim[0] == n0 && m_dim[1] == n1 &&
      m_dim[2] == n2 && m_dim[3] == n3 &&
      m_dim[4] == n4;
  }

  bool verify_dimension(const unsigned n0, const unsigned n1,
			const unsigned n2, const unsigned n3) const
  {
    enum { ok = stk::StaticAssert<4 == NumDim>::OK };
    return m_dim[0] == n0 && m_dim[1] == n1 &&
      m_dim[2] == n2 && m_dim[3] == n3;
  }

  bool verify_dimension(const unsigned n0, const unsigned n1,
			const unsigned n2) const
  {
    enum { ok = stk::StaticAssert<3 == NumDim>::OK };
    return m_dim[0] == n0 && m_dim[1] == n1 &&
      m_dim[2] == n2;
  }

  bool verify_dimension(const unsigned n0, const unsigned n1) const {
    enum { ok = stk::StaticAssert<2 == NumDim>::OK };
    return m_dim[0] == n0 && m_dim[1] == n1;
  }

  bool verify_dimension(const unsigned n0) const {
    enum { ok = stk::StaticAssert<1 == NumDim>::OK };
    return m_dim[0] == n0;
  }

  unsigned set_dim(const unsigned d[]) {
    m_stride[0] = 1;
    ArrayHelper<NumDim>::copy(d, m_dim);
    ArrayHelper<NumDim>::prefix(d, m_stride);
    return m_stride[NumDim];
  }

  unsigned set_dim(const unsigned n0, const unsigned n1,
		   const unsigned n2, const unsigned n3,
		   const unsigned n4, const unsigned n5,
		   const unsigned n6, const unsigned n7)
  {
    enum { ok = stk::StaticAssert<8 == NumDim>::OK };

    m_stride[0] = 1;
    m_stride[1] = m_stride[0] * (m_dim[0] = n0);
    m_stride[2] = m_stride[1] * (m_dim[1] = n1);
    m_stride[3] = m_stride[2] * (m_dim[2] = n2);
    m_stride[4] = m_stride[3] * (m_dim[3] = n3);
    m_stride[5] = m_stride[4] * (m_dim[4] = n4);
    m_stride[6] = m_stride[5] * (m_dim[5] = n5);
    m_stride[7] = m_stride[6] * (m_dim[6] = n6);
    m_stride[8] = m_stride[7] * (m_dim[7] = n7);
    return m_stride[NumDim];
  }

  unsigned set_dim(const unsigned n0, const unsigned n1,
		   const unsigned n2, const unsigned n3,
		   const unsigned n4, const unsigned n5,
		   const unsigned n6)
  {
    enum { ok = stk::StaticAssert<7 == NumDim>::OK };

    m_stride[0] = 1;
    m_stride[1] = m_stride[0] * (m_dim[0] = n0);
    m_stride[2] = m_stride[1] * (m_dim[1] = n1);
    m_stride[3] = m_stride[2] * (m_dim[2] = n2);
    m_stride[4] = m_stride[3] * (m_dim[3] = n3);
    m_stride[5] = m_stride[4] * (m_dim[4] = n4);
    m_stride[6] = m_stride[5] * (m_dim[5] = n5);
    m_stride[7] = m_stride[6] * (m_dim[6] = n6);
    return m_stride[NumDim];
  }

  unsigned set_dim(const unsigned n0, const unsigned n1,
		   const unsigned n2, const unsigned n3,
		   const unsigned n4, const unsigned n5)
  {
    enum { ok = stk::StaticAssert<6 == NumDim>::OK };

    m_stride[0] = 1;
    m_stride[1] = m_stride[0] * (m_dim[0] = n0);
    m_stride[2] = m_stride[1] * (m_dim[1] = n1);
    m_stride[3] = m_stride[2] * (m_dim[2] = n2);
    m_stride[4] = m_stride[3] * (m_dim[3] = n3);
    m_stride[5] = m_stride[4] * (m_dim[4] = n4);
    m_stride[6] = m_stride[5] * (m_dim[5] = n5);
    return m_stride[NumDim];
  }

  unsigned set_dim(const unsigned n0, const unsigned n1,
		   const unsigned n2, const unsigned n3,
		   const unsigned n4)
  {
    enum { ok = stk::StaticAssert<5 == NumDim>::OK };

    m_stride[0] = 1;
    m_stride[1] = m_stride[0] * (m_dim[0] = n0);
    m_stride[2] = m_stride[1] * (m_dim[1] = n1);
    m_stride[3] = m_stride[2] * (m_dim[2] = n2);
    m_stride[4] = m_stride[3] * (m_dim[3] = n3);
    m_stride[5] = m_stride[4] * (m_dim[4] = n4);
    return m_stride[NumDim];
  }

  unsigned set_dim(const unsigned n0, const unsigned n1,
		   const unsigned n2, const unsigned n3)
  {
    enum { ok = stk::StaticAssert<4 == NumDim>::OK };

    m_stride[0] = 1;
    m_stride[1] = m_stride[0] * (m_dim[0] = n0);
    m_stride[2] = m_stride[1] * (m_dim[1] = n1);
    m_stride[3] = m_stride[2] * (m_dim[2] = n2);
    m_stride[4] = m_stride[3] * (m_dim[3] = n3);
    return m_stride[NumDim];
  }

  unsigned set_dim(const unsigned n0, const unsigned n1,
		   const unsigned n2)
  {
    enum { ok = stk::StaticAssert<3 == NumDim>::OK };

    m_stride[0] = 1;
    m_stride[1] = m_stride[0] * (m_dim[0] = n0);
    m_stride[2] = m_stride[1] * (m_dim[1] = n1);
    m_stride[3] = m_stride[2] * (m_dim[2] = n2);
    return m_stride[NumDim];
  }

  unsigned set_dim(const unsigned n0, const unsigned n1)
  {
    enum { ok = stk::StaticAssert<2 == NumDim>::OK };

    m_stride[0] = 1;
    m_stride[1] = m_stride[0] * (m_dim[0] = n0);
    m_stride[2] = m_stride[1] * (m_dim[1] = n1);
    return m_stride[NumDim];
  }

  unsigned set_dim(const unsigned n0)
  {
    enum { ok = stk::StaticAssert<1 == NumDim>::OK };

    m_stride[0] = 1;
    m_stride[1] = m_stride[0] * (m_dim[0] = n0);
    return m_stride[NumDim];
  }

  unsigned set_dim(const SelfType & a)
  {
    ArrayHelper<NumDim>::copy(a.m_dim, m_dim);
    ArrayHelper<NumDim+1>::copy(a.m_stride, m_stride);
    return m_stride[NumDim];
  }

public:
  ~FArray()
  {
    m_ptr = NULL;
    ArrayHelper<NumDim>::fill((unsigned) 0, m_dim);
    ArrayHelper<NumDim+1>::fill((unsigned) 0, m_stride);
  }

  FArray()
    : m_ptr(NULL)
  {
    ArrayHelper<NumDim>::fill((unsigned) 0, m_dim);
    ArrayHelper<NumDim+1>::fill((unsigned) 0, m_stride);
  }

  FArray(const SelfType & a)
    : m_ptr(a.m_ptr)
  {
    ArrayHelper<NumDim>::copy(a.m_dim, m_dim);
    ArrayHelper<NumDim+1>::copy(a.m_stride, m_stride);
  }

  SelfType & operator=(SelfType const & a)
  {
    if (this != &a) {

      m_ptr = a.m_ptr;

      ArrayHelper<NumDim>::copy(a.m_dim, m_dim);
      ArrayHelper<NumDim+1>::copy(a.m_stride, m_stride);
    }
    return *this;
  }

  template<typename T>
  FArray(const FArray<T, Dimension> & a)
    : m_ptr(a.ptr())
  {
    ArrayHelper<NumDim>::copy(a.dimension(), m_dim);
    ArrayHelper<NumDim+1>::copy(a.stride(), m_stride);
  }

  FArray(value_type * const l_ptr,
	 const unsigned n0, const unsigned n1,
	 const unsigned n2, const unsigned n3,
	 const unsigned n4, const unsigned n5,
	 const unsigned n6, const unsigned n7)
    : m_ptr(l_ptr)
  {
    set_dim(n0, n1, n2, n3, n4, n5, n6, n7);
  }

  FArray(value_type * const l_ptr,
	 const unsigned n0, const unsigned n1,
	 const unsigned n2, const unsigned n3,
	 const unsigned n4, const unsigned n5,
	 const unsigned n6)
    : m_ptr(l_ptr)
  {
    set_dim(n0, n1, n2, n3, n4, n5, n6);
  }

  FArray(value_type * const l_ptr,
	 const unsigned n0, const unsigned n1,
	 const unsigned n2, const unsigned n3,
	 const unsigned n4, const unsigned n5)
    : m_ptr(l_ptr)
  {
    set_dim(n0, n1, n2, n3, n4, n5);
  }

  FArray(value_type * const l_ptr,
	 const unsigned n0, const unsigned n1,
	 const unsigned n2, const unsigned n3,
	 const unsigned n4)
    : m_ptr(l_ptr)
  {
    set_dim(n0, n1, n2, n3, n4);
  }

  FArray(value_type * const l_ptr,
	 const unsigned n0, const unsigned n1,
	 const unsigned n2, const unsigned n3)
    : m_ptr(l_ptr)
  {
    set_dim(n0, n1, n2, n3);
  }

  FArray(value_type * const l_ptr,
	 const unsigned n0, const unsigned n1,
	 const unsigned n2)
    : m_ptr(l_ptr)
  {
    set_dim(n0, n1, n2);
  }

  FArray(value_type * const l_ptr,
	 const unsigned n0, const unsigned n1)
    : m_ptr(l_ptr)
  {
    set_dim(n0, n1);
  }

  FArray(value_type * const l_ptr,
	 const unsigned n0)
    : m_ptr(l_ptr)
  {
    set_dim(n0);
  }

  FArray(value_type * const l_ptr,
	 const unsigned n[NumDim])
    : m_ptr(l_ptr)
  {
    set_dim(n);
  }

  void set(const SelfType & a) {
    m_ptr = a.m_ptr;

    ArrayHelper<NumDim>::copy(a.m_dim, m_dim);
    ArrayHelper<NumDim+1>::copy(a.m_stride, m_stride);
  }

  void set(value_type * const l_ptr,
	   const unsigned n0, const unsigned n1,
	   const unsigned n2, const unsigned n3,
	   const unsigned n4, const unsigned n5,
	   const unsigned n6, const unsigned n7)
  {
    m_ptr = l_ptr;
    set_dim(n0, n1, n2, n3, n4, n5, n6, n7);
  }

  void set(value_type * const l_ptr,
	   const unsigned n0, const unsigned n1,
	   const unsigned n2, const unsigned n3,
	   const unsigned n4, const unsigned n5,
	   const unsigned n6)
  {
    m_ptr = l_ptr;
    set_dim(n0, n1, n2, n3, n4, n5, n6);
  }

  void set(value_type * const l_ptr,
	   const unsigned n0, const unsigned n1,
	   const unsigned n2, const unsigned n3,
	   const unsigned n4, const unsigned n5)
  {
    m_ptr = l_ptr;
    set_dim(n0, n1, n2, n3, n4, n5);
  }

  void set(value_type * const l_ptr,
	   const unsigned n0, const unsigned n1,
	   const unsigned n2, const unsigned n3,
	   const unsigned n4) {
    m_ptr = l_ptr;
    set_dim(n0, n1, n2, n3, n4);
  }

  void set(value_type * const l_ptr,
	   const unsigned n0, const unsigned n1,
	   const unsigned n2, const unsigned n3)
  {
    m_ptr = l_ptr;
    set_dim(n0, n1, n2, n3);
  }

  void set(value_type * const l_ptr,
	   const unsigned n0, const unsigned n1,
	   const unsigned n2)
  {
    m_ptr = l_ptr;
    set_dim(n0, n1, n2);
  }

  void set(value_type * const l_ptr,
	   const unsigned n0, const unsigned n1) {
    m_ptr = l_ptr;
    set_dim(n0, n1);
  }

  void set(value_type * const l_ptr,
	   const unsigned n0) {
    m_ptr = l_ptr;
    set_dim(n0);
  }

  void set(value_type * const l_ptr,
	   const unsigned n[NumDim]) {
    m_ptr = l_ptr;
    set_dim(n);
  }

  Trunc dive(int i) {
    array_dimension_verify(0, i, m_dim[NumDim - 1]);

    value_type *l_ptr = m_ptr + i*m_stride[NumDim - 1];

    return Trunc(l_ptr, m_dim);
  }

  const Trunc dive(int i) const {
    array_dimension_verify(0, i, m_dim[NumDim - 1]);

    value_type *l_ptr = m_ptr + i*m_stride[NumDim - 1];

    return Trunc(l_ptr, m_dim);
  }

  iterator begin() {
    return m_ptr;
  }

  iterator end() {
    return m_ptr + m_stride[NumDim];
  }

  const_iterator begin() const {
    return m_ptr;
  }

  const_iterator end() const {
    return m_ptr + m_stride[NumDim];
  }

  void dimensions(const_iterator it, Index &index) const {
    ArrayHelper<NumDim - 1>::get_index(m_stride, &index[0], it - m_ptr);
  }

  //----------------------------------------
  // Copy contents of a compatible array

  template<typename T>
  void copy(const FArray<T, Dimension> & a) {
    ArrayHelper<NumDim>::copy(m_dim, m_ptr, m_stride, a.ptr(), a.stride());
  }

  template<typename T>
  void fill(const T & value) {
    ArrayHelper<NumDim>::fill(m_dim, m_ptr, m_stride, value);
  }

private:
  // Mutation (non-const methods) is not allowed so as to
  // provide derived classes with complete control over mutation.

protected:
  value_type *	        m_ptr;
  unsigned		m_dim[NumDim];
  unsigned		m_stride[NumDim + 1];
};


template<class ElementType, int Dimension, class A>
class FArrayContainer
  : public FArray<ElementType, Dimension>
{
public:
//  typedef typename FArray<ElementType, Dimension>::value_type value_type;

  typedef typename A::value_type value_type;
  typedef typename A::size_type size_type;
  typedef typename A::difference_type difference_type;

  typedef typename A::pointer pointer;
  typedef typename A::const_pointer const_pointer;
  typedef typename A::reference reference;
  typedef typename A::const_reference const_reference;

  typedef typename A::pointer iterator;
  typedef typename A::const_pointer const_iterator;

  typedef typename std::reverse_iterator<iterator> reverse_iterator;
  typedef typename std::reverse_iterator<const_iterator> const_reverse_iterator;

  typedef FArrayContainer<ElementType, Dimension> SelfType;

  typedef FArray<ElementType, Dimension> BaseType;

//  typedef typename BaseType::value_type value_type;

  enum { NumDim = BaseType::NumDim };

private:
  using BaseType::m_ptr;
  using BaseType::m_dim;
  using BaseType::m_stride;

  void resize_memory(const unsigned new_size) {
    if (m_capacity < new_size) {
      if ( m_capacity )
        m_allocator.deallocate(m_ptr, m_capacity);
      m_capacity = new_size;
      m_ptr = m_allocator.allocate(m_capacity);
    }
  }

public:
  ~FArrayContainer() {
    if (m_capacity) {
      m_allocator.deallocate(m_ptr, m_capacity);
      m_capacity = 0;
    }
  }

  //----------------------------------------
  // Constructors for initial view of contiguous memory.

  FArrayContainer()
    : BaseType(),
      m_capacity(0)
  {}

  FArrayContainer(const SelfType & a)
    : BaseType(),
      m_capacity(0)
  {
    resize_memory(BaseType::set_dim(a.m_dim));
    this->copy(a);
  }

  template<typename T>
  FArrayContainer(const FArray<T, Dimension> & a)
    : BaseType(),
      m_capacity(0)
  {
    resize_memory(BaseType::set_dim(a.dimension()));
    this->copy(a);
  }

  SelfType & operator=(const SelfType & a)
  {
    resize_memory(BaseType::set_dim(a.dimension()));
    this->copy(a);
    return *this;
  }

  template<typename T>
  SelfType & operator=(const FArray<T, Dimension> & a)
  {
    resize_memory(BaseType::set_dim(a.dimension()));
    this->copy(a);
    return *this;
  }

  //----------------------------------------

  FArrayContainer(const unsigned n0, const unsigned n1,
		  const unsigned n2, const unsigned n3,
		  const unsigned n4, const unsigned n5,
		  const unsigned n6, const unsigned n7)
    : BaseType(NULL, n0, n1, n2, n3, n4, n5, n6, n7),
      m_capacity(0)
  {
    resize_memory(m_stride[NumDim]);
    std::fill(m_ptr, m_ptr + m_stride[NumDim], 0);
  }

  FArrayContainer(const unsigned n0, const unsigned n1,
		  const unsigned n2, const unsigned n3,
		  const unsigned n4, const unsigned n5,
		  const unsigned n6)
    : BaseType(NULL, n0, n1, n2, n3, n4, n5, n6),
      m_capacity(0)
  {
    resize_memory(m_stride[NumDim]);
    std::fill(m_ptr, m_ptr + m_stride[NumDim], 0);
  }


  FArrayContainer(const unsigned n0, const unsigned n1,
		  const unsigned n2, const unsigned n3,
		  const unsigned n4, const unsigned n5)
    : BaseType(NULL, n0, n1, n2, n3, n4, n5),
      m_capacity(0)
  {
    resize_memory(m_stride[NumDim]);
    std::fill(m_ptr, m_ptr + m_stride[NumDim], 0);
  }


  FArrayContainer(const unsigned n0, const unsigned n1,
		  const unsigned n2, const unsigned n3,
		  const unsigned n4)
    : BaseType(NULL, n0, n1, n2, n3, n4),
      m_capacity(0)
  {
    resize_memory(m_stride[NumDim]);
    std::fill(m_ptr, m_ptr + m_stride[NumDim], 0);
  }


  FArrayContainer(const unsigned n0, const unsigned n1,
		  const unsigned n2, const unsigned n3)
    : BaseType(NULL, n0, n1, n2, n3),
      m_capacity(0)
  {
    resize_memory(m_stride[NumDim]);
    std::fill(m_ptr, m_ptr + m_stride[NumDim], 0);
  }


  FArrayContainer(const unsigned n0, const unsigned n1,
		  const unsigned n2)
    : BaseType(NULL, n0, n1, n2),
      m_capacity(0)
  {
    resize_memory(m_stride[NumDim]);
    std::fill(m_ptr, m_ptr + m_stride[NumDim], 0);
  }


  FArrayContainer(const unsigned n0, const unsigned n1)
    : BaseType(NULL, n0, n1),
      m_capacity(0)
  {
    resize_memory(m_stride[NumDim]);
    std::fill(m_ptr, m_ptr + m_stride[NumDim], 0);
  }


  FArrayContainer(const unsigned n0)
    : BaseType(NULL, n0),
      m_capacity(0)
  {
    resize_memory(m_stride[NumDim]);
    std::fill(m_ptr, m_ptr + m_stride[NumDim], 0);
  }


  FArrayContainer(const unsigned n[])
    : BaseType(NULL, n),
      m_capacity(0)
  {
    resize_memory(m_stride[NumDim]);
    std::fill(m_ptr, m_ptr + m_stride[NumDim], 0);
  }


  //----------------------------------------

  template<typename T>
  SelfType & resize(const FArray<T, Dimension> & a)
  {
    resize_memory(BaseType::set_dim(a.dimension()));
    return *this;
  }

  SelfType & resize(const unsigned n0, const unsigned n1,
		    const unsigned n2, const unsigned n3,
		    const unsigned n4, const unsigned n5,
		    const unsigned n6, const unsigned n7)
  {
    resize_memory(BaseType::set_dim(n0, n1, n2, n3, n4, n5, n6, n7));
    return *this;
  }

  SelfType & resize(const unsigned n0, const unsigned n1,
		    const unsigned n2, const unsigned n3,
		    const unsigned n4, const unsigned n5,
		    const unsigned n6)
  {
    resize_memory(BaseType::set_dim(n0, n1, n2, n3, n4, n5, n6));
    return *this;
  }

  SelfType & resize(const unsigned n0, const unsigned n1,
		    const unsigned n2, const unsigned n3,
		    const unsigned n4, const unsigned n5)
  {
    resize_memory(BaseType::set_dim(n0, n1, n2, n3, n4, n5));
    return *this;
  }

  SelfType & resize(const unsigned n0, const unsigned n1,
		    const unsigned n2, const unsigned n3,
		    const unsigned n4)
  {
    resize_memory(BaseType::set_dim(n0, n1, n2, n3, n4));
    return *this;
  }

  SelfType & resize(const unsigned n0, const unsigned n1,
		    const unsigned n2, const unsigned n3)
  {
    resize_memory(BaseType::set_dim(n0, n1, n2, n3));
    return *this;
  }


  SelfType & resize(const unsigned n0, const unsigned n1,
		    const unsigned n2)
  {
    resize_memory(BaseType::set_dim(n0, n1, n2));
    return *this;
  }

  SelfType & resize(const unsigned n0, const unsigned n1) {
    resize_memory(BaseType::set_dim(n0, n1));
    return *this;
  }

  SelfType & resize(const unsigned n0) {
    resize_memory(BaseType::set_dim(n0));
    return *this;
  }

  SelfType & resize(const unsigned n[]) {
    resize_memory(BaseType::set_dim(n));
    return *this;
  }

private:
  A             m_allocator;
  unsigned      m_capacity;
};

///
/// @}
///

} // namespace sierra

#endif // STK_UTIL_DIAG_FArray_h
