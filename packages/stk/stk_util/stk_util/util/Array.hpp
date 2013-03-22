/**
 *    Copyright 2005 - 2008 Sandia Corporation.
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

#ifndef STK_UTIL_UTIL_Array_h
#define STK_UTIL_UTIL_Array_h

#include <cstddef>
#include <utility>
#include <algorithm>
#include <stdexcept>
#include <stk_util/util/TypeList.hpp>
#include <stk_util/diag/FArray.hpp>

namespace sierra {

///
/// @addtogroup ArrayDetail
/// @{
///

/**
 * @class Array
 *
 * @brief Multidimensional array view of contiguous memory.
 *
 * @par Dimension tagging
 *
 * @brief Multidimensional array of contiguous memory.  The memory is not owned by the
 * array, but container access semantics are enforced, i.e. const Array<> elements cannot
 * be assigned to.
 *
 *    class X {}; class Y {}; class Z {};
 *    Array<double, X, Y, Z> grid;
 *
 * describes a three dimensional array with the three indices associated with 'X', 'Y',
 * and 'Z' respectively.
 *
 * @par Indexing via 'operator()'
 *
 * Elements of an array can be accessed via the () operator, where an index is provided
 * for each dimension of the array.  The number of dimensions is enforced at compile time;
 * however, enforcement of the range of the indices is only performed at runtime in debug
 * mode, i.e. 'NDEBUG' is not defined.
 *
 *    FArray<double, W, X, Y, Z> a(mem, n0, n1, n2, n3);
 *    a(i0, i1, i2, i3) = 24;
 *
 * @par Construction with a block of memory
 *
 *  An array is created by providing it with memory and dimensions by which that memory is
 *  to be viewed.
 *
 *    FArray(pointer, ndim0, ndim1, ndim2, ...);
 *
 * @par Shallow copy construction
 *
 * The copy constructor is a shallow copy.  The constructed array is just another a view
 * into the same memory as the input array.  This is the preferred method of passing
 * arrays as access required less memory indirection.
 */

template< class ElementType,
	  class Tag0,
	  class Tag1 = TypeListEnd,
	  class Tag2 = TypeListEnd,
	  class Tag3 = TypeListEnd,
	  class Tag4 = TypeListEnd,
	  class Tag5 = TypeListEnd,
	  class Tag6 = TypeListEnd,
	  class Tag7 = TypeListEnd>
class Array;

//----------------------------------------------------------------------
/**
 * @class sierra::ArrayContainer
 *
 * @brief Extend Array with deep copy assignment and resize operations.
 */

template< class ElementType,
	  class Tag0,
	  class Tag1 = TypeListEnd,
	  class Tag2 = TypeListEnd,
	  class Tag3 = TypeListEnd,
	  class Tag4 = TypeListEnd,
	  class Tag5 = TypeListEnd,
	  class Tag6 = TypeListEnd,
	  class Tag7 = TypeListEnd,
          class A = std::allocator<ElementType> >
class ArrayContainer;


//----------------------------------------------------------------------
// A typeless array is invalid...

template< class Tag0,
	  class Tag1,
	  class Tag2,
	  class Tag3,
	  class Tag4,
	  class Tag5,
	  class Tag6,
	  class Tag7>
class Array<void, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7> {};

template <class ElementType>
class Array<ElementType, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd>
{
public:
  typedef ElementType element_type;

  typedef Array< element_type, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd> SelfType;

//  typedef Array< const element_type, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd > Const;

  typedef typename MakeTypeList<TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd, TypeListEnd>::type TagList;

  /** Number of dimensions */

  enum { NumDim = 0 };
};



//----------------------------------------------------------------------
//----------------------------------------------------------------------
template< class ElementType,
	  class Tag0,
	  class Tag1,
	  class Tag2,
	  class Tag3,
	  class Tag4,
	  class Tag5,
	  class Tag6,
	  class Tag7>
class Array : public FArray<ElementType, TypeListLength<typename MakeTypeList<Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7>::type>::value>
{
public:
  typedef ElementType element_type;

  typedef Array< element_type, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7> SelfType;

  typedef Array< const element_type, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7 > Const;

  typedef typename MakeTypeList<Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7>::type TagList;

  typedef FArray<ElementType, TypeListLength<typename MakeTypeList<Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7>::type>::value> BaseType;

  /** Number of dimensions */

  enum { NumDim = TypeListLength<TagList>::value };

  typedef Array< ElementType,
		 typename TypeListAt<typename TypeListEraseAt<TagList, NumDim - 1>::list_type, 0>::type,
		 typename TypeListAt<typename TypeListEraseAt<TagList, NumDim - 1>::list_type, 1>::type,
		 typename TypeListAt<typename TypeListEraseAt<TagList, NumDim - 1>::list_type, 2>::type,
		 typename TypeListAt<typename TypeListEraseAt<TagList, NumDim - 1>::list_type, 3>::type,
		 typename TypeListAt<typename TypeListEraseAt<TagList, NumDim - 1>::list_type, 4>::type,
		 typename TypeListAt<typename TypeListEraseAt<TagList, NumDim - 1>::list_type, 5>::type,
		 typename TypeListAt<typename TypeListEraseAt<TagList, NumDim - 1>::list_type, 6>::type,
		 TypeListEnd> Trunc;

  //----------------------------------------
  /** Dimensions defined at contruction */

  const unsigned * dimension() const {
    return BaseType::dimension();
  }

  const unsigned * stride() const {
    return BaseType::stride();
  }

  template<unsigned I>
  unsigned dimension() const {
    enum { ok = stk::StaticAssert< I < NumDim >::OK };

    return BaseType::m_dim[I];
  }

  unsigned dimension( const unsigned i ) const {
    this->array_dimension_verify(0, i, NumDim );
    return BaseType::m_dim[i];
  }

  unsigned stride( const unsigned i ) const {
    this->array_dimension_verify(0, i, NumDim );
    return BaseType::m_stride[i];
  }

  template<unsigned I>
  unsigned stride() const {
    enum { ok = stk::StaticAssert< I < NumDim >::OK };
    return BaseType::m_stride[I];
  }

  template<class Tag>
  unsigned dimension() const{
    enum { I  = TypeListIndex<TagList, Tag>::value };
    enum { ok = stk::StaticAssert< 0 <= I >::OK };
    return BaseType::m_dim[I];
  }

  template<class Tag, unsigned Ordinal>
  unsigned dimension() const {
    enum { I  = TypeListIndex<TagList, Tag, Ordinal>::value };
    enum { ok = stk::StaticAssert< 0 <= I >::OK };
    return BaseType::m_dim[I];
  }

  template<class Tag>
  unsigned stride() const {
    enum { I  = TypeListIndex<TagList, Tag>::value };
    enum { ok = stk::StaticAssert< 0 <= I >::OK };
    return BaseType::m_stride[I];
  }

  template<class Tag, unsigned Ordinal>
  unsigned stride() const {
    enum { I  = TypeListIndex<TagList, Tag, Ordinal>::value };
    enum { ok = stk::StaticAssert< 0 <= I >::OK };
    return BaseType::m_stride[I];
  }

  bool operator == ( const SelfType & a ) const {
    return ArrayHelper<NumDim>::equal( BaseType::m_dim, a.m_dim ) &&
      ArrayHelper<NumDim>::equal(BaseType::m_dim, BaseType::m_ptr, BaseType::m_stride, a.m_ptr, a.m_stride);
  }

  template<typename T>
  bool operator == (const Array<T, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7> & a) const {
    return ArrayHelper<NumDim>::equal( BaseType::m_dim, a.dimension() ) &&
      ArrayHelper<NumDim>::equal(BaseType::m_dim, BaseType::m_ptr, BaseType::m_stride, a.ptr(), a.stride());
  }

  bool operator != ( const SelfType & a ) const {
    return ! operator == ( a );
  }

  template<typename T>
  bool operator != ( const Array<T, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7> & a ) const {
    return ! operator == ( a );
  }

public:
  ~Array()
  {}

  Array() : BaseType()
  {}

  Array( const SelfType & a )
    : BaseType( a )
  {}

  template<typename T>
  Array( const Array<T, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7> & a )
    : BaseType( a )
  {}

  Array( element_type * const in_ptr,
	 const unsigned n0, const unsigned n1,
	 const unsigned n2, const unsigned n3,
	 const unsigned n4, const unsigned n5,
	 const unsigned n6, const unsigned n7 )
    : BaseType( in_ptr, n0, n1, n2, n3, n4, n5, n6, n7)
  {}

  Array( element_type * const in_ptr,
	 const unsigned n0, const unsigned n1,
	 const unsigned n2, const unsigned n3,
	 const unsigned n4, const unsigned n5,
	 const unsigned n6 )
    : BaseType( in_ptr, n0, n1, n2, n3, n4, n5, n6)
  {}

  Array( element_type * const in_ptr,
	 const unsigned n0, const unsigned n1,
	 const unsigned n2, const unsigned n3,
	 const unsigned n4, const unsigned n5 )
    : BaseType( in_ptr, n0, n1, n2, n3, n4, n5)
  {}

  Array( element_type * const in_ptr,
	 const unsigned n0, const unsigned n1,
	 const unsigned n2, const unsigned n3,
	 const unsigned n4 )
    : BaseType( in_ptr, n0, n1, n2, n3, n4)
  {}

  Array( element_type * const in_ptr,
	 const unsigned n0, const unsigned n1,
	 const unsigned n2, const unsigned n3 )
    : BaseType( in_ptr, n0, n1, n2, n3)
  {}

  Array( element_type * const in_ptr,
	 const unsigned n0, const unsigned n1,
	 const unsigned n2 )
    : BaseType( in_ptr, n0, n1, n2)
  {}

  Array( element_type * const in_ptr,
	 const unsigned n0, const unsigned n1 )
    : BaseType( in_ptr, n0, n1)
  {}

  Array( element_type * const in_ptr,
	 const unsigned n0 )
    : BaseType( in_ptr, n0)
  {}

  Array( element_type * const in_ptr,
	 const unsigned n[NumDim] )
    : BaseType( in_ptr, n){}

  void set( const SelfType & a ) {
    BaseType::m_ptr = a.m_ptr;

    ArrayHelper<NumDim  >::copy( a.m_dim, BaseType::m_dim );
    ArrayHelper<NumDim+1>::copy( a.m_stride, BaseType::m_stride );
  }

  void set( element_type * const in_ptr,
	    const unsigned n0, const unsigned n1,
	    const unsigned n2, const unsigned n3,
	    const unsigned n4, const unsigned n5,
	    const unsigned n6, const unsigned n7 ) {
    BaseType::set(in_ptr, n0, n1, n2, n3, n4, n5, n6, n7);
  }

  void set( element_type * const in_ptr,
	    const unsigned n0, const unsigned n1,
	    const unsigned n2, const unsigned n3,
	    const unsigned n4, const unsigned n5,
	    const unsigned n6 ) {
    BaseType::set(in_ptr, n0, n1, n2, n3, n4, n5, n6);
  }

  void set( element_type * const in_ptr,
	    const unsigned n0, const unsigned n1,
	    const unsigned n2, const unsigned n3,
	    const unsigned n4, const unsigned n5 ) {
    BaseType::set(in_ptr, n0, n1, n2, n3, n4, n5);
  }

  void set( element_type * const in_ptr,
	    const unsigned n0, const unsigned n1,
	    const unsigned n2, const unsigned n3,
	    const unsigned n4 ) {
    BaseType::set(in_ptr, n0, n1, n2, n3, n4);
  }

  void set( element_type * const in_ptr,
	    const unsigned n0, const unsigned n1,
	    const unsigned n2, const unsigned n3 ) {
    BaseType::set(in_ptr, n0, n1, n2, n3);
  }

  void set( element_type * const in_ptr,
	    const unsigned n0, const unsigned n1,
	    const unsigned n2 ) {
    BaseType::set(in_ptr, n0, n1, n2);
  }

  void set( element_type * const in_ptr,
	    const unsigned n0, const unsigned n1 ) {
    BaseType::set(in_ptr, n0, n1);
  }

  void set( element_type * const in_ptr,
	    const unsigned n0 ) {
    BaseType::set(in_ptr, n0);
  }

  void set( element_type * const in_ptr,
	    const unsigned n[NumDim] ) {
    BaseType::set(in_ptr, n);
  }

  Trunc dive(int i) {
    this->array_dimension_verify(0, i, BaseType::m_dim[NumDim - 1] );

    element_type *calc_ptr = BaseType::m_ptr + i*BaseType::m_stride[NumDim - 1];

    return Trunc(calc_ptr, BaseType::m_dim);
  }

  const Trunc dive(int i) const {
    this->array_dimension_verify(0, i, BaseType::m_dim[NumDim - 1] );

    element_type *calc_ptr = BaseType::m_ptr + i*BaseType::m_stride[NumDim - 1];

    return Trunc(calc_ptr, BaseType::m_dim);
  }

  template<typename T>
  void copy( const Array<T, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7> & a ) {
    ArrayHelper<NumDim>::copy(BaseType::m_dim, BaseType::m_ptr, BaseType::m_stride, a.ptr(), a.stride() );
  }

  template<typename T>
  void fill( const T & value ) {
    ArrayHelper<NumDim>::fill(BaseType::m_dim, BaseType::m_ptr, BaseType::m_stride, value);
  }

private:
  // Mutation (non-const methods) is not allowed so as to
  // provide derived classes with complete control over mutation.

  SelfType & operator = ( SelfType const & a );
};

template< class ElementType,
	  class Tag0,
	  class Tag1,
	  class Tag2,
	  class Tag3,
	  class Tag4,
	  class Tag5,
	  class Tag6,
	  class Tag7,
          class A>
class ArrayContainer
  : public Array<ElementType, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7>
{
public:
  typedef ArrayContainer<ElementType, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7> SelfType;

  typedef Array<ElementType, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7> BaseType;

  typedef typename BaseType::element_type element_type;

  typedef typename BaseType::TagList TagList;

  enum { NumDim = BaseType::NumDim };

private:
  using BaseType::m_ptr;
  using BaseType::m_dim;
  using BaseType::m_stride;

  void resize_memory( const unsigned new_size ) {
    if ( m_capacity < new_size ) {
      if ( m_capacity )
        m_allocator.deallocate(m_ptr, m_capacity);
      m_capacity = new_size;
      m_ptr = m_allocator.allocate(m_capacity);
    }
  }

public:
  ~ArrayContainer() {
    if ( m_capacity ) {
      m_allocator.deallocate(m_ptr, m_capacity);
      m_capacity = 0;
    }
  }

  //----------------------------------------
  // Constructors for initial view of contiguous memory.

  ArrayContainer( )
    : BaseType( ),
      m_capacity(0)
  {}

  ArrayContainer( const SelfType & a )
    : BaseType(),
      m_capacity(0)
  {
    resize_memory( BaseType::set_dim( a.m_dim ) );
    this->copy(a);
  }

  template<typename T>
  ArrayContainer( const Array<T, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7> & a )
    : BaseType(), m_capacity(0)
  {
    resize_memory( BaseType::set_dim( a.dimension() ) );
    this->copy(a);
  }

  SelfType & operator = ( const SelfType & a ) {
    resize_memory( BaseType::set_dim( a.dimension() ) );
    this->copy(a);
    return *this;
  }

  template<typename T>
  SelfType & operator =( const Array<T, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7> & a ) {
    resize_memory( BaseType::set_dim( a.dimension() ) );
    this->copy(a);
    return *this;
  }

  //----------------------------------------

  ArrayContainer( const unsigned n0, const unsigned n1,
		  const unsigned n2, const unsigned n3,
		  const unsigned n4, const unsigned n5,
		  const unsigned n6, const unsigned n7 )
    : BaseType(NULL, n0, n1, n2, n3, n4, n5, n6, n7 ),
      m_capacity(0)
  {
    resize_memory( m_stride[NumDim] );
  }

  ArrayContainer( const unsigned n0, const unsigned n1,
		  const unsigned n2, const unsigned n3,
		  const unsigned n4, const unsigned n5,
		  const unsigned n6 )
    : BaseType( NULL, n0, n1, n2, n3, n4, n5, n6 ),
      m_capacity(0)
  {
    resize_memory( m_stride[NumDim] );
  }

  ArrayContainer( const unsigned n0, const unsigned n1,
		  const unsigned n2, const unsigned n3,
		  const unsigned n4, const unsigned n5 )
    : BaseType( NULL, n0, n1, n2, n3, n4, n5 ),
      m_capacity(0)
  {
    resize_memory( m_stride[NumDim] );
  }

  ArrayContainer( const unsigned n0, const unsigned n1,
		  const unsigned n2, const unsigned n3,
		  const unsigned n4 )
    : BaseType( NULL, n0, n1, n2, n3, n4 ),
      m_capacity(0)
  {
    resize_memory( m_stride[NumDim] );
  }

  ArrayContainer( const unsigned n0, const unsigned n1,
		  const unsigned n2, const unsigned n3 )
    : BaseType( NULL, n0, n1, n2, n3 ),
      m_capacity(0)
  {
    resize_memory( m_stride[NumDim] );
  }

  ArrayContainer( const unsigned n0, const unsigned n1,
		  const unsigned n2 )
    : BaseType( NULL, n0, n1, n2 ),
      m_capacity(0)
  {
    resize_memory( m_stride[NumDim] );
  }

  ArrayContainer( const unsigned n0, const unsigned n1 )
    : BaseType( NULL, n0, n1 ),
      m_capacity(0)
  {
    resize_memory( m_stride[NumDim] );
  }

  ArrayContainer( const unsigned n0 )
    : BaseType( NULL, n0 ),
      m_capacity(0)
  {
    resize_memory( m_stride[NumDim] );
  }

  ArrayContainer( const unsigned n[] )
    : BaseType( NULL, n ),
      m_capacity(0)
  {
    resize_memory( m_stride[NumDim] );
  }

  //----------------------------------------

  template<typename T>
  SelfType & resize( const Array<T, Tag0, Tag1, Tag2, Tag3, Tag4, Tag5, Tag6, Tag7> & a ) {
    resize_memory( BaseType::set_dim( a.dimension() ) );
    return *this;
  }

  SelfType & resize( const unsigned n0, const unsigned n1,
		     const unsigned n2, const unsigned n3,
		     const unsigned n4, const unsigned n5,
		     const unsigned n6, const unsigned n7 )
  {
    resize_memory(BaseType::set_dim(n0, n1, n2, n3, n4, n5, n6, n7));
    return *this;
  }

  SelfType & resize( const unsigned n0, const unsigned n1,
		     const unsigned n2, const unsigned n3,
		     const unsigned n4, const unsigned n5,
		     const unsigned n6 )
  {
    resize_memory(BaseType::set_dim(n0, n1, n2, n3, n4, n5, n6));
    return *this;
  }

  SelfType & resize( const unsigned n0, const unsigned n1,
		     const unsigned n2, const unsigned n3,
		     const unsigned n4, const unsigned n5 )
  {
    resize_memory(BaseType::set_dim(n0, n1, n2, n3, n4, n5));
    return *this;
  }

  SelfType & resize( const unsigned n0, const unsigned n1,
		     const unsigned n2, const unsigned n3,
		     const unsigned n4 )
  {
    resize_memory(BaseType::set_dim(n0, n1, n2, n3, n4));
    return *this;
  }

  SelfType & resize( const unsigned n0, const unsigned n1,
		     const unsigned n2, const unsigned n3 )
  {
    resize_memory(BaseType::set_dim(n0, n1, n2, n3));
    return *this;
  }


  SelfType & resize( const unsigned n0, const unsigned n1,
		     const unsigned n2 )
  {
    resize_memory(BaseType::set_dim(n0, n1, n2));
    return *this;
  }

  SelfType & resize( const unsigned n0, const unsigned n1 )
  {
    resize_memory(BaseType::set_dim(n0, n1));
    return *this;
  }

  SelfType & resize( const unsigned n0 )
  {
    resize_memory(BaseType::set_dim(n0));
    return *this;
  }

  SelfType & resize( const unsigned n[] )
  {
    resize_memory( BaseType::set_dim(n) );
    return *this;
  }

  //----------------------------------------

private:
  A             m_allocator;
  unsigned      m_capacity;
};

///
/// @}
///

} // namespace sierra

#endif // STK_UTIL_UTIL_Array_h
