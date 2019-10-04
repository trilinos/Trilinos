// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_ElementOp_hpp
#define percept_ElementOp_hpp

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>


  namespace percept
  {

    class ElementOp
    {
    public:
      /// innermost operation of an element-based loop; return value of true forces the enclosing loop to terminate and this class'
      ///   derived classes can return info back to the loop invoker
      virtual bool operator()(const stk::mesh::Entity element, stk::mesh::FieldBase *field,  const stk::mesh::BulkData& bulkData)=0;
      virtual void init_elementOp()=0;
      virtual void fini_elementOp()=0;
      virtual ~ElementOp() {}
    };

  }//namespace percept

#endif
