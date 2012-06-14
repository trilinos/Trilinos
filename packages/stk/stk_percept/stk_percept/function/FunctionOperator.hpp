#ifndef stk_encr_FunctionOperator_hpp
#define stk_encr_FunctionOperator_hpp

//#define HAVE_INTREPID_DEBUG 1

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/fem/FEMMetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>

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
      stk::mesh::Selector *m_selector;
      mesh::Part *m_part;
      bool m_own_selector;
    public:
      FunctionOperator(mesh::BulkData& bulkData, mesh::Part *part = 0) : m_bulkData(bulkData), m_selector(0), m_part(part), m_own_selector(false)
      {
        if (!part)
          {
            m_part = &stk::mesh::fem::FEMMetaData::get(bulkData).universal_part();
          }
        m_selector = new stk::mesh::Selector(*m_part);
        m_own_selector = true;
      }

      FunctionOperator(mesh::BulkData& bulkData, mesh::Selector *selector) : m_bulkData(bulkData), m_selector(selector), m_part(0), m_own_selector(false)
      {
        if (!selector)
          {
            m_part = &stk::mesh::fem::FEMMetaData::get(bulkData).universal_part();
            m_selector = new stk::mesh::Selector(*m_part);
            m_own_selector = true;
          }
      }

      virtual ~FunctionOperator() {
        if (m_own_selector)
          delete m_selector;
      }

      //stk::mesh::Selector *get_selector() { return m_selector; }

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
