#ifndef stk_encr_StringFunction_hpp
#define stk_encr_StringFunction_hpp

#include <string>

#include <stk_expreval/Evaluator.hpp>
#include <map>
#include <vector>

#include <stk_percept/function/Function.hpp>
#include <stk_percept/Name.hpp>

namespace stk
{
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
      derivative(MDArrayString& deriv_spec);

      Teuchos::RCP<Function > gradient(int spatialDim=3);

      //========================================================================================================================
      // low-level interface
      StringFunction(const char *function_string, 
                     Name name,
                     Dimensions domain_dimensions,
                     Dimensions codomain_dimensions,
                     unsigned integration_order = 0);

      StringFunction(const StringFunction& s);

      void resolve(stk::expreval::VariableMap::iterator & var_it);

      Teuchos::RCP<Function > derivative_test(MDArrayString& deriv_spec);
      Teuchos::RCP<Function > derivative_test_fd(MDArrayString& deriv_spec, double eps=1.e-6);

      virtual void operator()(MDArray& in, MDArray& out, double time_value_optional=0.0);
      virtual void operator()(MDArray& in, MDArray& out, const stk::mesh::Entity& element, const MDArray& parametric_coords, double time_value_optional=0.0);
      virtual void operator()(MDArray& in, MDArray& out, const stk::mesh::Bucket& bucket, const MDArray& parametric_coords, double time_value_optional=0.0);

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

      const stk::mesh::Entity * m_element;
      const stk::mesh::Bucket * m_bucket;
      MDArray m_parametric_coordinates;
      bool m_have_element;
      bool m_have_bucket;

      int m_spatialDim;

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
}//namespace stk
#endif
