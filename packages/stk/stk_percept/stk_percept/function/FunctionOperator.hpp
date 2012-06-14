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
      //mesh::Part *m_part;
      bool m_own_selector;
    public:

      FunctionOperator(mesh::BulkData& bulkData, mesh::Part *part = 0) : m_bulkData(bulkData), m_selector(0), m_own_selector(false)
      {
        init(part);
      }

      FunctionOperator(mesh::BulkData& bulkData, mesh::Selector *selector) : m_bulkData(bulkData), m_selector(selector), m_own_selector(false)
      {
        init(selector);
      }

      void init(mesh::Part *part)
      {
        if (!part)
          {
            m_selector = new stk::mesh::Selector(stk::mesh::fem::FEMMetaData::get(m_bulkData).universal_part());
          }
        else
          {
            m_selector = new stk::mesh::Selector(*part);
          }
        m_own_selector = true;
      }

      void init(mesh::Selector *selector)
      {
        if (!selector)
          {
            m_selector = new stk::mesh::Selector(stk::mesh::fem::FEMMetaData::get(m_bulkData).universal_part());
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
