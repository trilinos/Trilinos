// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#ifndef stk_encr_StringFunction_hpp
#define stk_encr_StringFunction_hpp

#include <string>

#include <stk_expreval/Evaluator.hpp>
#include <map>
#include <vector>

#include <percept/function/Function.hpp>
#include <percept/Name.hpp>

  namespace percept
  {

    class StringFunction : public Function, public stk::expreval::VariableMap::Resolver
    {
    public:


      //========================================================================================================================
      // high-level interface

      StringFunction(const char *function_string,
                     Name name = Name("noname"),
                     int domain_dimension = 3,
                     int codomain_dimension = 1,
                     unsigned integration_order = 0);


      std::string
      getFunctionString() { return m_func_string; }

      void
      set_gradient_strings(std::string gstring[3], int len);

      void
      set_gradient_strings(MDArrayString& gstring);

      virtual Teuchos::RCP<Function >
      derivative(MDArrayString& deriv_spec) override;

      Teuchos::RCP<Function > gradient(int spatialDim=3) override;

      void set_current_time(double time) {current_time=time; using_current_time=true;}
      void unset_time() {using_current_time=true;}
      //========================================================================================================================
      // low-level interface
      StringFunction(const char *function_string,
                     Name name,
                     Dimensions domain_dimensions,
                     Dimensions codomain_dimensions,
                     unsigned integration_order = 0);

      StringFunction(const StringFunction& s);

      void resolve(stk::expreval::VariableMap::iterator & var_it) override;

      Teuchos::RCP<Function > derivative_test(MDArrayString& deriv_spec);
      Teuchos::RCP<Function > derivative_test_fd(MDArrayString& deriv_spec, double eps=1.e-6);

      virtual void operator()(MDArray& in, MDArray& out, double time_value_optional=0.0) override;
      virtual void operator()(MDArray& in, MDArray& out, const stk::mesh::Entity element, const MDArray& parametric_coords, double time_value_optional=0.0) override;
      virtual void operator()(MDArray& in, MDArray& out, const stk::mesh::Bucket& bucket, const MDArray& parametric_coords, double time_value_optional=0.0) override;

    private:
      void evalFunctions(MDArray& in, double time_value_optional=0.0);

      std::string m_func_string;
      stk::expreval::Eval m_functionExpr;
      //std::vector<std::string> m_gradient_string;
      std::string m_gradient_string;
      //Expr::Eval gradientExpr;
      //Expr::Eval dotExpr;

      //Local variables
      double m_x,m_y,m_z,m_t;
      std::vector<double> m_v;

      //Map of function names to their evalutation values
      std::map<Function *, std::vector<double> > m_func_to_value;

      stk::mesh::Entity m_element;
      const stk::mesh::Bucket * m_bucket;
      MDArray m_parametric_coordinates;
      bool m_have_element;
      bool m_have_bucket;

      int m_spatialDim;

      double current_time;
      bool using_current_time;

      void init();



    };
#if 1

    // new StringFunction = lhs OP rhs
    StringFunction operator-(StringFunction& lhs, StringFunction& rhs);
    StringFunction operator+(StringFunction& lhs, StringFunction& rhs);
    StringFunction operator/(StringFunction& lhs, StringFunction& rhs);
    StringFunction operator*(StringFunction& lhs, StringFunction& rhs);

    // unary ops
    StringFunction operator-(StringFunction& lhs);


#endif

  }//namespace percept

#endif
