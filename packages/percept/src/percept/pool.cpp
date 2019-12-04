// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#include "pool.h"

Pool::Pool( size_t granularity, size_t size )
  : m_granularity( granularity ), m_size( size ), m_used( 0 ), m_overflow( 0 )
{
	if( m_size > 0 )
	{
		m_storage.resize(m_size*granularity);
		m_slots.resize(m_size);

		for( size_t i = 0; i < m_size; ++i )
			m_slots[i] = reinterpret_cast<void*>( m_storage.data() + i*granularity );
	}
}

Pool::~Pool()
{
  // can't destroy a pool with outstanding allocations
  assert( m_used == 0 && m_overflow == 0);
}

void* Pool::Allocate()
{
	if( m_used < m_size )
	{
		return m_slots[m_used++];
	}
	else
	{
		++m_overflow;
		return reinterpret_cast<void*>( new char[m_granularity] );
	}
}

void Pool::Deallocate( void* block )
{
        // null pointer argument
	assert( block );
	if( IsFromPool( block ) )
	{
    	        // internal error
		assert( m_used > 0 );
		m_slots[--m_used] = block;
	}
	else
	{
    	        // internal error
		assert( m_overflow > 0 );
		delete[] reinterpret_cast<char*>( block );
		--m_overflow;
	}
}

