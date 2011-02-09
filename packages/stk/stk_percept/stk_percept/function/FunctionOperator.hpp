#ifndef stk_encr_FunctionOperator_hpp
#define stk_encr_FunctionOperator_hpp

//#define HAVE_INTREPID_DEBUG 1

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>

#include <stk_percept/function/Function.hpp>

using namespace Intrepid;

namespace stk
{
  namespace percept
  {
    //class FunctionWithIntrepidRequest;

    /**
     */
    class FunctionOperator
    {
    protected:
      mesh::BulkData& m_bulkData;
      mesh::Part *m_part;
    public:
      FunctionOperator(mesh::BulkData& bulkData, mesh::Part *part = 0) : m_bulkData(bulkData), m_part(part)
      {
        if (!part)
          {
            m_part = &stk::mesh::MetaData::get(bulkData).universal_part();
          }
      }
      virtual ~FunctionOperator() {}

      virtual void operator()(Function& integrand, Function& result) = 0;
#if 0
      void norm_l2( Function& integrand, Function& result);
      void integrate( Function& integrand, Function& result);


      void integrate(FunctionWithIntrepidRequest& integrand, Function& result);
#endif

    };

  }//namespace percept
}//namespace stk
#endif
