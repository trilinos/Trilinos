// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef stk_encr_CompositeFunction_hpp
#define stk_encr_CompositeFunction_hpp

#include <string>
#include <sstream>
#include <stdexcept>
#include <iostream>

#include <percept/function/Function.hpp>

  namespace percept
  {
    //typedef GenericVector<GenericFunction> GenericFunctionVector;

    //typedef std::vector<GenericFunction *> GenericFunctionVector;

    class CompositeFunction : public Function 
    {
    public:
      /// compose two functions to be able to apply in turn as in func_2(func_1(x)), or more specifically as:
      ///   func_1(domain, codomain_temp);
      //    func_2(codomain_temp, codomain);

      /// Note that since this is a Function also, one can make multiple compositions e.g. h(g(f(x))) by
      /// CompositeFunction g_of_f (f, g)
      /// CompositeFunction h_of_g_of_f (g_of_f, h);
      /// The first function in the list is always applied first.
      CompositeFunction(const char *name, Function& func_1, Function& func_2,
                        Dimensions domain_dimensions = Dimensions(),
                        Dimensions codomain_dimensions = Dimensions(),
                        unsigned integration_order = 0) : Function(name, domain_dimensions, codomain_dimensions, integration_order),
                                                          m_func1(func_1), m_func2(func_2) 
      {
        EXCEPTWATCH;
//         std::cout << "func_1 = " << func_1 << std::endl;
//         std::cout << "func_2 = " << func_2 << std::endl;
//         std::cout << "domain_dimensions= " << domain_dimensions << std::endl;
//         std::cout << "codomain_dimensions= " << codomain_dimensions << std::endl;
        setDomainDimensions(func_1.getDomainDimensions());
        setCodomainDimensions(func_2.getCodomainDimensions());
      }

      virtual void operator()(MDArray& domain, MDArray& codomain, double time_value_optional=0.0)
      {
        EXCEPTWATCH;
        m_func1(domain, codomain, time_value_optional);
        MDArray input("input",codomain.layout());
        Kokkos::deep_copy(input,codomain);
        m_func2(input, codomain, time_value_optional);
      }
      
      virtual void operator()(MDArray& domain, MDArray& codomain, const stk::mesh::Entity element, const MDArray& parametric_coords, double time_value_optional=0.0)
      {
        EXCEPTWATCH;
        m_func1(domain, codomain, element, parametric_coords, time_value_optional);
        MDArray input("input",codomain.layout());
        Kokkos::deep_copy(input,codomain);
        m_func2( input, codomain, element, parametric_coords, time_value_optional);
      }


      virtual void operator()(MDArray& domain, MDArray& codomain, const stk::mesh::Bucket& bucket, const MDArray& parametric_coords, double time_value_optional=0.0)
      {
        EXCEPTWATCH;
        m_func1(domain, codomain, bucket, parametric_coords, time_value_optional);
        MDArray input("input",codomain.layout());
        Kokkos::deep_copy(input,codomain);
        m_func2( input, codomain, bucket, parametric_coords, time_value_optional);
      }

    protected:

      Function& m_func1;
      Function& m_func2;

    };

  }//namespace percept

#endif
