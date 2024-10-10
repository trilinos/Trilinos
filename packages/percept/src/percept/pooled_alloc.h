#ifndef INCL_POOLED_ALLOC
#define INCL_POOLED_ALLOC

#include "pool.h"
#include <iostream>
#include <typeinfo>
#include <cassert>
#include <list>
#include <set>
#include <map>

//! A standards-compliant pooled allocator.
template<typename T>
class PooledAllocator
{
public:
	typedef size_t		size_type;			//!< A type that can represent the size of the largest object in the allocation model.
	typedef ptrdiff_t	difference_type;	//!< A type that can represent the difference between any two pointers in the allocation model.

	typedef T			value_type;			//!< Identical to T.
	typedef T*			pointer;			//!< Pointer to T;
	typedef T const*	const_pointer;		//!< Pointer to const T.
	typedef T&			reference;			//!< Reference to T.
	typedef T const&	const_reference;	//!< Reference to const T.

	//! A struct to construct an allocator for a different type.
	template<typename U>
	struct rebind { typedef PooledAllocator<U> other; };

	//! Creates a pooled allocator to the given pool.
	/*! This is non-explicit for ease of use.
	*/
	PooledAllocator( Pool* pool = 0 ) : m_pool( pool )
	{
    if (sizeof( T ) > m_pool->GetGranularity() )
      {
        std::cout << " sizeof( T ) = " <<  sizeof( T ) << "  <= m_pool->GetGranularity() " <<  m_pool->GetGranularity()
                  << " typeid= " << typeid(T).name()
                  << std::endl;
        throw std::runtime_error("Pool granularity too small");
      }
		assert( !m_pool || sizeof( T ) <= m_pool->GetGranularity() );
	}

	//! Creates a pooled allocator to the argument's pool.
	/*! If the argument has no pool, then this allocator will allocate off the heap.
	*/
	template<typename U>
	PooledAllocator( PooledAllocator<U> const& arg ) : m_pool( arg.m_pool )
	{
    if (sizeof( T ) > m_pool->GetGranularity() )
      {
#ifndef NDEBUG
        std::cout << " sizeof( T ) = " <<  sizeof( T ) << "  <= m_pool->GetGranularity() " <<  m_pool->GetGranularity()
                  << " typeid= " << typeid(T).name()
                  << std::endl;
#endif
        m_pool = new Pool(sizeof(T), arg.m_pool->GetSize());
#ifndef NDEBUG
        std::cout << " after sizeof( T ) = " <<  sizeof( T ) << "  <= m_pool->GetGranularity() " <<  m_pool->GetGranularity()
                  << " typeid= " << typeid(T).name()
                  << std::endl;
#endif
      }
		assert( !m_pool || sizeof( T ) <= m_pool->GetGranularity() );
	}

	//! The largest value that can meaningfully passed to allocate.
	size_type max_size() const { return 0xffffffff; }

	//! Memory is allocated for \c count objects of type \c T but objects are not constructed.
	/*! This function may raise an appropriate exception. The result is a random access iterator.
	*/
	pointer allocate( size_type count, typename std::allocator_traits<T>::const_void_pointer /*hint*/ = 0 ) const
	{
		if( count == 1 && m_pool )
		{
			return reinterpret_cast<T*>( m_pool->Allocate() );
		}
		else
		{
			return reinterpret_cast<T*>( new char[ count*sizeof( T ) ] );
		}
	}

	//! Deallocates memory allocated by allocate.
	/*! All \c count objects in the area pointed by \c block must be destroyed prior to this call.
		\c count must match the value passed to allocate to obtain this memory. Does not throw exceptions.
		\c block must not be null.
	*/
	void deallocate( pointer block, size_type count ) const throw()
	{
		assert( block && "null pointer argument" );
		if( count == 1 && m_pool )
		{
			m_pool->Deallocate( block );
		}
		else
		{
			delete[] reinterpret_cast<char*>( block );
		}
	}

	//! Constructs an element of \c T at the given pointer.
	/*! Effect: new( element ) T( arg )
	*/
	void construct( pointer element, T const& arg )
	{
		new( element ) T( arg );
	}

	//! Destroys an element of \c T at the given pointer.
	/*! Effect: element->~T()
	*/
	void destroy( pointer element )
	{
		element->~T();
		(void)element; // FIXME: genuinely bizarre 'unused parameter' bug with vs2003
	}

	//! Returns the address of the given reference.
	pointer address( reference element ) const
	{
		return &element;
	}

	//! Returns the address of the given reference.
	const_pointer address( const_reference element ) const
	{
		return &element;
	}

	//! The pool for this allocator.
	Pool* m_pool;
};

//! A specialisation of the pooled allocator for the void type.
template<>
class PooledAllocator<void>
{
public:
	typedef size_t		size_type;			//!< A type that can represent the size of the largest object in the allocation model.
	typedef ptrdiff_t	difference_type;	//!< A type that can represent the difference between any two pointers in the allocation model.

	typedef void		value_type;			//!< Identical to void.
	typedef void*		pointer;			//!< Pointer to void;
	typedef void const*	const_pointer;		//!< Pointer to const void.

	//! A struct to construct an allocator for a different type.
	template<typename U>
	struct rebind { typedef PooledAllocator<U> other; };

	//! Creates a pooled allocator with no pool.
	/*! This allocator will allocate off the heap.
	*/
	PooledAllocator() : m_pool( 0 ) {}

	//! Creates a pooled allocator to the given pool.
	PooledAllocator( Pool* pool ) : m_pool( pool ) {}

	//! The pool for this allocator.
	Pool* m_pool;
};

//! Returns true if objects allocated from one pool can be deallocated from the other.
template<typename T, typename U>
bool operator==( PooledAllocator<T> const& left, PooledAllocator<U> const& right )
{
	return left.m_pool == right.m_pool;
}

//! Returns true if objects allocated from one pool cannot be deallocated from the other.
template<typename T, typename U>
bool operator!=( PooledAllocator<T> const& left, PooledAllocator<U> const& right )
{
	return left.m_pool != right.m_pool;
}

//! Template typedef std::map<..., PooledAllocator> to PooledList<Key, Value>::Type.
template<typename Key, typename Value, class Traits = std::less<Key> >
struct PooledMap
{
	typedef std::map<Key, Value, Traits, PooledAllocator<std::pair<Key, Value> > > Type;
};

//! Template typedef std::list<..., PooledAllocator> to PooledList<Value>::Type.
template<typename Value>
struct PooledList
{
	typedef std::list<Value, PooledAllocator<Value> > Type;
};

//! Template typedef std::set<..., PooledAllocator> to PooledList<Value>::Type.
template<typename Value, typename Traits = std::less<Value> >
struct PooledSet
{
	typedef std::set<Value, Traits, PooledAllocator<Value> > Type;
};

#endif // ndef INCL_POOLED_ALLOC
