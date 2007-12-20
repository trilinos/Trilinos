/*--------------------------------------------------------------------*/
/*    Copyright 2005 Sandia Corporation.                              */
/*    Under the terms of Contract DE-AC04-94AL85000, there is a       */
/*    non-exclusive license for use of this work by or on behalf      */
/*    of the U.S. Government.  Export of this program may require     */
/*    a license from the United States Government.                    */
/*--------------------------------------------------------------------*/

#include "fei_macros.hpp"

#include "fei_Pool_alloc.hpp"
#include "fei_Pool.hpp"

#include <cstdlib>
#include <map>

struct fei_map_holder {
  std::map<size_t,fei_Pool*> fei_mem_pools;

  ~fei_map_holder()
  {
    std::map<size_t,fei_Pool*>::iterator
      iter = fei_mem_pools.begin(),
      iter_end = fei_mem_pools.end();
    for(; iter != iter_end; ++iter) {
      delete iter->second;
    }
  }
};

static fei_map_holder fmh;

fei_Pool* get_fei_mem_pool(size_t n)
{
  std::map<size_t,fei_Pool*>::iterator
    iter = fmh.fei_mem_pools.find(n);

  if (iter == fmh.fei_mem_pools.end()) {
    fei_Pool* new_pool = new fei_Pool(n);
    fmh.fei_mem_pools.insert(std::make_pair(n,new_pool));
    return(new_pool);
  }

  return(iter->second);
}

