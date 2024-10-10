#ifndef INCL_CORE_POOL
#define INCL_CORE_POOL

#include <cstddef>
#include <assert.h>
#include <vector>

//! Simple pool class.
class Pool
{
public:
	//! Allocates a pool with size elements each of granularity bytes.
	Pool( size_t granularity, size_t size );

	//! Checks for emptiness before destructing.
	~Pool();

        void Reset() { m_used = 0; m_overflow = 0; }

	//! Gets the pool granularity.
	size_t GetGranularity() const { return m_granularity; }

	//! Gets the pool size in elements.
	size_t GetSize() const { return m_size; }

	//! Gets the number of elements allocated from pooled storage.
	size_t GetUsed() const { return m_used; }

	//! Gets the number of elements allocated from non-pooled storage.
	size_t GetOverflow() const { return m_overflow; }


	//! Allocates memory from the pool without construction.
	void* Allocate();

	//! Deallocates memory from the pool without destruction.
	void Deallocate( void* block );


	//! Constructs an object from the pool.
	template<typename T>
	T* Construct()
	{
		assert( sizeof( T ) <= m_granularity );
		T* block = reinterpret_cast<T*>( Allocate() );
		return new( block ) T;
	}

	//! Destructs an object back into the pool.
	template<typename T>
	void Destroy( T* instance )
	{
		assert( sizeof( T ) <= m_granularity );
		instance->~T();
		Deallocate( instance );
	}

private:
	//! Returns true if the given instance is from pooled storage.
	bool IsFromPool( void const* instance ) const
	{
		char const* block = reinterpret_cast<char const*>( instance );
		return m_storage.data() <= block && block < ( m_storage.data() + m_size*m_granularity );
	}

	size_t m_granularity;	//!< The size of each element in the pool in bytes.
	size_t m_size;			//!< The number of elements in pooled storage.
	size_t m_used;			//!< The number of pooled allocations.
	size_t m_overflow;		//!< The number of non-pooled allocations.

        std::vector<char> m_storage;	//!< The pool storage.
        std::vector<void*> m_slots;	//!< The free list.
};

#endif // ndef INCL_CORE_POOL
