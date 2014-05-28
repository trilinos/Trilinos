#ifndef stk_encr_ConstantFunction_hpp
#define stk_encr_ConstantFunction_hpp

#include <stk_percept/function/Function.hpp>
#include <stk_percept/function/internal/HasValue.hpp>

namespace stk
{
  namespace percept
  {
    class ConstantFunction : public Function, public HasValue<double>
    {
      double m_value;
    public:
      ConstantFunction(double value,
                       const char *name,
                       Dimensions domain_dimensions = Dimensions(),
                       Dimensions codomain_dimensions = Dimensions(),
                       unsigned integration_order = 0) : Function(name, domain_dimensions, codomain_dimensions, integration_order), m_value(value) {}

      using Function::operator();
      virtual void operator()(MDArray& domain, MDArray& codomain, double time_value_optional=0.0)
      {
        // set all values (regardless of rank of codomain) to m_value
        codomain.initialize(m_value);
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
        int stride = codomain.dimension(codomain.rank()-1);
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
}//namespace stk
#endif
