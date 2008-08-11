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

#ifndef SACADO_FAD_MEMPOOLMANAGER_HPP
#define SACADO_FAD_MEMPOOLMANAGER_HPP

#include <map>

#include "Sacado_Fad_MemPool.hpp"

namespace Sacado {

  namespace Fad {

    //! Class to manage memory pools for different Fad dimensions
    template <typename T>
    class MemPoolManager {

    public:
      
      //! Constructor
      MemPoolManager(unsigned int nfad);

      //! Destructor
      ~MemPoolManager();

      //! Get memory pool for supplied dimension \c dim
      MemPool* getMemoryPool(unsigned int dim);

    private:

      //! Private to prohibit copying
      MemPoolManager(const MemPoolManager&);

      //! Private to prohibit copying
      MemPoolManager& operator=(const MemPoolManager&);
      
    protected:

      //! Number of Fad objects per chunk
      unsigned int num_fad;

      //! Typename of memory pool map
      typedef std::map<unsigned int, MemPool*> MapType;

      //! Map of memory pools
      MapType poolMap;

    };

  } // namespace Fad

} // namespace Sacado

// Include implementation
#include "Sacado_Fad_MemPoolManagerImp.hpp"

#endif // SACADO_FAD_MEMPOOLMANAGER_HPP
