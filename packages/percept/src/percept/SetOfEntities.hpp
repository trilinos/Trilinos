// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef SetOfEntities_hpp
#define SetOfEntities_hpp

#define PERCEPT_USE_STD_SET 0
#define PERCEPT_USE_STD_POOLED_SET 1
#define PERCEPT_USE_STD_USET 0

#if PERCEPT_USE_STD_USET
#include <unordered_set>
#include <stk_mesh/base/HashEntityAndEntityKey.hpp>
#endif

#if PERCEPT_USE_STD_POOLED_SET
#include <boost/scoped_ptr.hpp>
#include <boost/pool/pool_alloc.hpp>
#endif

  namespace percept {

#if PERCEPT_USE_STD_SET
    //typedef std::set<stk::mesh::Entity, stk::mesh::EntityLess> SetOfEntitiesBase;
    typedef std::set<stk::mesh::Entity> SetOfEntitiesBase;
    struct SetOfEntities : public SetOfEntitiesBase
    {
      SetOfEntities() : SetOfEntitiesBase() {}
      SetOfEntities(stk::mesh::BulkData& bulk) : SetOfEntitiesBase() {}
    };
#endif

#if PERCEPT_USE_STD_USET
    typedef std::unordered_set<stk::mesh::Entity> SetOfEntitiesBase;
    struct SetOfEntities : public SetOfEntitiesBase
    {
      SetOfEntities() : SetOfEntitiesBase() {}
      SetOfEntities(stk::mesh::BulkData& bulk) : SetOfEntitiesBase() {}
    };

#endif

#if PERCEPT_USE_STD_POOLED_SET

    // this is the expected number of elements that are node neighbors of any element
    enum { PERCEPT_POOLED_SET_POOL_SIZE = 100 };


    template<typename Value, typename Less = std::less<Value> >
    struct SetOfEntitiesBase
    {
      typedef std::set<Value, Less,
                       boost::fast_pool_allocator<Value, boost::default_user_allocator_new_delete,
                                                  boost::details::pool::null_mutex,
                                                  PERCEPT_POOLED_SET_POOL_SIZE
                                                  >
                       > Type;

      typedef Less less;

    };

    struct SetOfEntities : public SetOfEntitiesBase<stk::mesh::Entity >::Type
    {
      SetOfEntities() {}
      SetOfEntities(stk::mesh::BulkData& bulk) {}
    };


    //typedef std::set<stk::mesh::Entity, stk::mesh::EntityLess> ElementUnrefineCollection;
    //typedef elements_to_be_destroyed_type ElementUnrefineCollection;

#endif

  }

#endif
