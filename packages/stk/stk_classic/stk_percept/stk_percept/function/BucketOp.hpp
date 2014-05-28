#ifndef stk_percept_BucketOp_hpp
#define stk_percept_BucketOp_hpp

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Bucket.hpp>


namespace stk
{
  namespace percept
  {

    class BucketOp
    {
    public:
      /// innermost operation of an bucket-based loop; return value of true forces the enclosing loop to terminate and this class'
      ///   derived classes can return info back to the loop invoker
      virtual bool operator()(const stk::mesh::Bucket& bucket, 
                              stk::mesh::FieldBase *field,  const mesh::BulkData& bulkData)=0;
      virtual ~BucketOp() {}
    };

  }//namespace percept
}//namespace stk

#endif
