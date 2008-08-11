// $Id$ 
// $Source$ 
// @HEADER
// ***********************************************************************
// 
//                           Sacado Package
//                 Copyright (2006) Sandia Corporation
// 
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact David M. Gay (dmgay@sandia.gov) or Eric T. Phipps
// (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#ifndef SACADO_FAD_MEMPOOL_HPP
#define SACADO_FAD_MEMPOOL_HPP

namespace Sacado {

  namespace Fad {

    //! Memory pool
    class MemPool {

    public:

      /*! 
       * \brief Constructor.  \c elem_size is the size of elements, \c n_elem 
       * is the number of elements per chunk.  \c pre_alloc is the number
       * of chunks to pre-allocate.
       */
      MemPool(unsigned int elem_size, unsigned int n_elem, 
	      unsigned int pre_alloc = 0);

      //! Destructor
      ~MemPool();

      //! Allocate a new element
      void* alloc();

      //! Free an element
      void free(void *b);

      //! Return number of allocated chunks
      unsigned int numChunks() const { return num_chunks; }

    private:
      
      //! Private to prohibit copying
      MemPool(const MemPool&);

      //! Private to prohibit copying
      MemPool& operator=(const MemPool&);

      //! Allocate a new chunk
      inline void grow();

    protected:

      //! Represents a memory element
      struct Link { 
	Link* next; 
      };

      //! Represents a memory chunk
      struct Chunk {
	Chunk *next;
	char *mem;
	//~Chunk() { delete [] mem; }
      };

      //! Size of elements in a chunk
      const unsigned int esize;

      //! Number of elements per chunk
      const unsigned int n;

      //! Size of memory chunks
      const unsigned int csize;

      //! Pointer to memory chunks
      Chunk *chunks;

      //! Pointer to first free link
      Link *head;

      //! Number of allocated chunks
      unsigned int num_chunks;
    };

  } // namespace Fad

} // namespace Sacado

// Include implementation
#include "Sacado_Fad_MemPoolImp.hpp"

#endif // SACADO_FAD_MEMPOOL_HPP
