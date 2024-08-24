// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef stk_encr_FunctionOperator_hpp
#define stk_encr_FunctionOperator_hpp

#include <stk_mesh/base/Field.hpp>
#include <stk_mesh/base/MetaData.hpp>
#include <stk_mesh/base/BulkData.hpp>
#include <stk_mesh/base/Part.hpp>
#include <stk_mesh/base/Selector.hpp>

#include <percept/function/Function.hpp>

  namespace percept
  {
    /**
     */
    class FunctionOperator
    {
    protected:
      stk::mesh::BulkData& m_bulkData;
      stk::mesh::Selector *m_selector;
      //mesh::Part *m_part;
      bool m_own_selector;
    public:

      FunctionOperator(stk::mesh::BulkData& bulkData, stk::mesh::Part *part = 0) : m_bulkData(bulkData), m_selector(0), m_own_selector(false)
      {
        init(part);
      }

      FunctionOperator(stk::mesh::BulkData& bulkData, stk::mesh::Selector *selector) : m_bulkData(bulkData), m_selector(selector), m_own_selector(false)
      {
        init(selector);
      }

      void init(stk::mesh::Part *part)
      {
        if (m_own_selector)
          {
            VERIFY_OP_ON(m_selector, !=, 0, "FunctionOperator::init");
            delete m_selector;
          }
        if (!part)
          {
            m_selector = new stk::mesh::Selector(m_bulkData.mesh_meta_data().universal_part());
          }
        else
          {
            m_selector = new stk::mesh::Selector(*part);
          }
        m_own_selector = true;
      }

      void init(stk::mesh::Selector *selector)
      {
        if (!selector)
          {
            if (m_own_selector)
              {
                VERIFY_OP_ON(m_selector, !=, 0, "FunctionOperator::init");
                delete m_selector;
              }
            m_selector = new stk::mesh::Selector(m_bulkData.mesh_meta_data().universal_part());
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

#endif
