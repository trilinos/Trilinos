// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef stk_encr_ConstantFunction_hpp
#define stk_encr_ConstantFunction_hpp

#include <percept/function/Function.hpp>
#include <percept/function/internal/HasValue.hpp>

  namespace percept
  {
    class ConstantFunction : public Function, public HasValue<double>
    {
      double m_value;
    public:
      ConstantFunction(double value = 0.0,
                       const char *name = NULL,
                       Dimensions domain_dimensions = Dimensions(),
                       Dimensions codomain_dimensions = Dimensions(),
                       unsigned integration_order = 0) : Function(name, domain_dimensions, codomain_dimensions, integration_order), m_value(value) {}

      using Function::operator();
      virtual void operator()(MDArray& domain, MDArray& codomain, double time_value_optional=0.0)
      {
        // set all values (regardless of rank of codomain) to m_value
        Kokkos::deep_copy(codomain, m_value);
      }

      double& getValue() { return m_value; }
      void setValue(double& v) { m_value = v; }

    };

    class ConstantFunctionVec : public Function, public HasValue<std::vector<double> >
    {
      std::vector<double> m_value;
    public:
      ConstantFunctionVec(std::vector<double>& value,
                          const char *name,
                          Dimensions domain_dimensions = Dimensions(),
                          Dimensions codomain_dimensions = Dimensions(),
                          unsigned integration_order = 0) : Function(name, domain_dimensions, codomain_dimensions, integration_order), m_value(value) {}

      using Function::operator();
      virtual void operator()(MDArray& domain, MDArray& codomain, double time_value_optional=0.0)
      {
        if (codomain.rank() <= 0)
          {
            throw std::runtime_error("ConstantFunctionVec::operator() codomain rank is <= 0");
          }
        int stride = codomain.extent_int(codomain.rank()-1);
        if (stride <= 0)
          {
            throw std::runtime_error("ConstantFunctionVec::operator() codomain stride is <= 0");
          }
        if (stride != (int)m_value.size())
          {
            throw std::runtime_error("ConstantFunctionVec::operator() codomain stride is not same as value");
          }
        int sz = codomain.size();
        int k = 0;
        for (int i = 0; i < sz/stride; i++)
          {
            for (int j = 0; j < stride; j++)
              {
                codomain[k+j] = m_value[j];
              }
            k += stride;
          }
      }

      std::vector<double>& getValue() { return m_value; }
      void setValue(std::vector<double>& v) { m_value = v; }

    };

  }//namespace percept

#endif
