#ifndef stk_percept_PartOp_hpp
#define stk_percept_PartOp_hpp

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>


namespace stk
{
  namespace percept
  {

    class PartOp
    {
    public:
      /// innermost operation of an part-based loop; return value of true forces the enclosing loop to terminate and this class'
      ///   derived classes can return info back to the loop invoker
      virtual bool operator()(const stk::mesh::Part& part, 
                              mesh::Selector& select_owned,   // select which buckets to use
                              stk::mesh::FieldBase *field,  
                              const mesh::BulkData& bulkData)=0;
    };

  }//namespace percept
}//namespace stk

#endif
