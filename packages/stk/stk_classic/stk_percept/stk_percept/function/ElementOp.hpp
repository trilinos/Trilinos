#ifndef stk_percept_ElementOp_hpp
#define stk_percept_ElementOp_hpp

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Entity.hpp>


namespace stk
{
  namespace percept
  {

    class ElementOp
    {
    public:
      /// innermost operation of an element-based loop; return value of true forces the enclosing loop to terminate and this class'
      ///   derived classes can return info back to the loop invoker
      virtual bool operator()(const stk::mesh::Entity& element, stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)=0;
      virtual void init_elementOp()=0;
      virtual void fini_elementOp()=0;
      virtual ~ElementOp() {}
    };

  }//namespace percept
}//namespace stk

#endif
