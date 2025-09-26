// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef SetOfEntities_hpp
#define SetOfEntities_hpp

  namespace percept {

    template<typename Value, typename Less = std::less<Value>, typename Allocator = std::allocator<Value> >
    struct SetOfEntitiesBase
    {
      typedef std::set<Value, Less, Allocator> Type;

      typedef Less less;
    };

    struct SetOfEntities : public SetOfEntitiesBase<stk::mesh::Entity >::Type
    {
      SetOfEntities() {}
      SetOfEntities(stk::mesh::BulkData& /*bulk*/) {}
    };

  }

#endif
