// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef percept_PartOp_hpp
#define percept_PartOp_hpp

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>


  namespace percept
  {

    class PartOp
    {
    public:
      /// innermost operation of an part-based loop; return value of true forces the enclosing loop to terminate and this class'
      ///   derived classes can return info back to the loop invoker
      virtual bool operator()(const stk::mesh::Part& part, 
                              stk::mesh::Selector& select_owned,   // select which buckets to use
                              stk::mesh::FieldBase *field,  
                              const stk::mesh::BulkData& bulkData)=0;
    };

  }//namespace percept

#endif
