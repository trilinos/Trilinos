#include <string>
#include <sstream>
#include <stdexcept>
#include <sstream>
#include <stk_percept/function/StringFunction.hpp>
#include <stk_percept/ExceptionWatch.hpp>
#include <stk_mesh/base/Bucket.hpp>
#include <boost/lexical_cast.hpp>

namespace stk
{
  namespace percept
  {
    StringFunction::StringFunction(const char *function_string, 
                                   Name name,
                                   int domain_dimension,
                                   int codomain_dimension,
                                   unsigned integration_order) :
      Function(name.getName().c_str(), Dimensions(domain_dimension), Dimensions(codomain_dimension), integration_order),
      m_func_string(function_string), m_functionExpr(*this), m_gradient_string("") 
    {
      init();
    }

    StringFunction::StringFunction(const char *function_string, 
                                   Name name,
                                   Dimensions domain_dimensions,
                                   Dimensions codomain_dimensions,
                                   unsigned integration_order) :
      Function(name.getName().c_str(), domain_dimensions, codomain_dimensions, integration_order),
      m_func_string(function_string), m_functionExpr(*this), m_gradient_string("") 
    {
      init();
    }

    StringFunction::StringFunction(const StringFunction& s) : 
      Function(s.m_name.c_str(), s.m_domain_dimensions, s.m_codomain_dimensions, s.m_integration_order),
      m_func_string(s.m_func_string), m_functionExpr(*this), m_gradient_string(s.m_gradient_string)
    {
      init();
    }

    void StringFunction::init()
    {
      // FIXME for tensor-valued
      int nOutDim = m_codomain_dimensions.back();
      m_v.resize(nOutDim);
      //m_g.resize(nOutDim);
      m_functionExpr.setExpression(m_func_string.c_str()).parse();
      m_element=0;
      m_bucket=0;
      m_parametric_coordinates = MDArray(1,3); 
      m_have_element = false;
      m_have_bucket = false;
    }

    void StringFunction::setGradientStrings(std::string gstring[3], int len)
    {
      m_gradient_string = "";
      for (int i = 0; i < len; i++)
        m_gradient_string += "v["+boost::lexical_cast<std::string>(i)+"]= "+gstring[i]+";";
    }

    // new StringFunction = lhs OP rhs
    StringFunction operator-(StringFunction& lhs, StringFunction& rhs)
    {
      std::string newFuncString = "(" + lhs.getFunctionString() + ") - (" + rhs.getFunctionString() + ")";
      return StringFunction(newFuncString.c_str());
    }

    StringFunction operator+(StringFunction& lhs, StringFunction& rhs)
    {
      std::string newFuncString = "(" + lhs.getFunctionString() + ") + (" + rhs.getFunctionString() + ")";
      return StringFunction(newFuncString.c_str());
    }

    StringFunction operator*(StringFunction& lhs, StringFunction& rhs)
    {
      std::string newFuncString = "(" + lhs.getFunctionString() + ") * (" + rhs.getFunctionString() + ")";
      return StringFunction(newFuncString.c_str());
    }

    StringFunction operator/(StringFunction& lhs, StringFunction& rhs)
    {
      std::string newFuncString = "(" + lhs.getFunctionString() + ") / (" + rhs.getFunctionString() + ")";
      return StringFunction(newFuncString.c_str());
    }

    // unary minus
    StringFunction operator-(StringFunction& lhs)
    {
      std::string newFuncString = "-(" + lhs.getFunctionString() + ")";
      return StringFunction(newFuncString.c_str());
    }

    void StringFunction::evalFunctions(MDArray& inp, double time_value_optional)
    {
      EXCEPTWATCH;
      //static MDArray out(3); // FIXME
      std::map<Function *, std::vector<double> >::iterator it = m_func_to_value.begin();
      const std::map<Function *, std::vector<double> >::iterator it_end = m_func_to_value.end();

      //Evaluate the function and store the answer in it's bound location
      for(; it != it_end; ++it)
        {
          Function& func = *it->first;

          // check to avoid infinite recursion
          if (&func == this)
            {
              std::ostringstream msg;
              msg << "StringFunction::evalFunctions: infinite recursion by self-referential string function, name= " << m_name;
              VERIFY_1(msg.str());
            }

          VERIFY_OP(func.getNewDomain().rank(), ==, 1, "StringFunction::evalFunctions: functions must be defined with domain/codomain of rank 1 for now");
          VERIFY_OP(func.getNewCodomain().rank(), ==, 1, "StringFunction::evalFunctions: functions must be defined with domain/codomain of rank 1 for now");
          int numCells = 1;
          //                  double max_val = std::numeric_limits<double>::max();

          MDArray f_inp = MDArray(numCells, func.getNewDomain().dimension(0));
          MDArray f_out = MDArray(numCells, func.getNewCodomain().dimension(0));

          int nInDim = last_dimension(f_inp);
          // FIXME - do we really need the extra array?
          for (int iDim = 0; iDim < nInDim; iDim++)
            {
              f_inp(0, iDim) = inp( iDim);
            }

          //VERIFY_OP(!m_have_element, || , !m_have_bucket, "StringFunction::evalFunctions: can't have both element and bucket");
          if (m_have_element)
            {
              func(f_inp, f_out, *m_element, m_parametric_coordinates, time_value_optional); 
            }
          //           else if (m_have_bucket)
          //             {
          //               func(f_inp, f_out, *m_bucket, m_parametric_coordinates, time_value_optional); 
          //             }
          else
            {
              func(f_inp, f_out); 
            }

          int nOutDim = last_dimension(f_out);

          for (int iDim = 0; iDim < nOutDim; iDim++)
            {
              (it->second)[iDim] = f_out(0, iDim);
            }

        }
    }

    static void first_dimensions(MDArray& arr, int arr_offset, int *n_points)
    {
      for (int ii = 0; ii < 3; ii++)
        {
          n_points[ii] = 1;
        }
      for (int ii = 0; ii < arr_offset; ii++)
        {
          n_points[ii] = arr.dimension(ii);
        }
    }

    // replace all occurrences of input_string with output_string in source_string

    Teuchos::RCP<Function > StringFunction::derivative_test(MDArrayString& deriv_spec)
    {
#if 1
      {
        std::string str1="x+y+x";
        std::string eps_string = boost::lexical_cast<std::string>(1.e-10);

        //replace(str1, "x", "(x+1.e-10)");
        Util::replace(str1, "x", "(x+"+eps_string+")");
        std::cout << "str1= " << str1 << std::endl;
        //Util::pause(true,"tmp str1");
      }
#endif      


      std::string fstr = m_func_string;
      // FIXME this is just for a simple test and only works for linear functions
      int outputDim = deriv_spec.dimension(0);
      static std::string s_xyzt[] = {"x", "y", "z", "t"};
      std::string out_str;
      for (int iresult = 0; iresult < deriv_spec.dimension(0); iresult++)
        {
          fstr = m_func_string;
          for (int jderiv = 0; jderiv < deriv_spec.dimension(1); jderiv++)
            {
              //std::cout << "deriv_spec= " << deriv_spec(iresult, jderiv)  << " fstr= " << fstr << std::endl;
              char xyzt = deriv_spec(iresult, jderiv)[0];
              std::replace( fstr.begin(), fstr.end(), xyzt, '1' );
              //std::cout << "xyzt= " << xyzt << " fstr= " << fstr << std::endl;
            }
          for (int jderiv = 0; jderiv < 4; jderiv++)
            {
              char xyzt = s_xyzt[jderiv][0];
              std::replace( fstr.begin(), fstr.end(), xyzt, '0' );
              //std::cout << "xyzt= " << xyzt << " fstr= " << fstr << std::endl;
            }
          if (deriv_spec.dimension(0) > 1)
            {
              out_str += "v["+boost::lexical_cast<std::string>(iresult)+"]= " + fstr+";";
            }
          else
            {
              out_str = fstr;
            }
        }
      std::string fname = getName();
      std::string new_name = "deriv_"+fname;
      //std::cout << "fname= " << fname << " new_name= " << new_name << " out_str= " << out_str << std::endl;
      //Util::pause(true, "tmp:: StringFunction::derivative");
      return Teuchos::rcp(new StringFunction(out_str.c_str(), Name(new_name), Dimensions(3), Dimensions(outputDim) ));
    }



    Teuchos::RCP<Function > StringFunction::derivative_test_fd(MDArrayString& deriv_spec, double eps)
    {
      std::string eps_string = boost::lexical_cast<std::string>(eps);
      std::string fstr = m_func_string;
      std::string fstr_p = m_func_string;
      std::string fstr_m = m_func_string;
      // FIXME this is just for a simple test and only works for linear functions
      int outputDim = deriv_spec.dimension(0);
      static std::string s_xyzt[] = {"x", "y", "z", "t"};
      std::string out_str;
      for (int iresult = 0; iresult < deriv_spec.dimension(0); iresult++)
        {
          fstr = m_func_string;
          fstr_p = m_func_string;
          fstr_m = m_func_string;
          for (int jderiv = 0; jderiv < deriv_spec.dimension(1); jderiv++)
            {
              //std::cout << "deriv_spec= " << deriv_spec(iresult, jderiv)  << " fstr= " << fstr << std::endl;
              std::string rep = "("+deriv_spec(iresult, jderiv)+"+"+eps_string+")";
              Util::replace(fstr_p, deriv_spec(iresult, jderiv), rep);
              rep = "("+deriv_spec(iresult, jderiv)+"-"+eps_string+")";
              Util::replace(fstr_m, deriv_spec(iresult, jderiv), rep);
              fstr = "(("+fstr_p+")-("+fstr_m+"))/(2.0*"+eps_string+")";
              //std::cout << "xyzt= " << xyzt << " fstr= " << fstr << std::endl;
            }
          if (deriv_spec.dimension(0) > 1)
            {
              out_str += "v["+boost::lexical_cast<std::string>(iresult)+"]= " + fstr+";";
            }
          else
            {
              out_str = fstr;
            }
        }
      std::string fname = getName();
      std::string new_name = "deriv_"+fname;
      //std::cout << "fname= " << fname << " new_name= " << new_name << " out_str= " << out_str << std::endl;
      //Util::pause(true, "tmp:: StringFunction::derivative");
      return Teuchos::rcp(new StringFunction(out_str.c_str(), Name(new_name), Dimensions(3), Dimensions(outputDim) ));
    }

    Teuchos::RCP<Function > StringFunction::derivative(MDArrayString& deriv_spec)
    {
      if (m_gradient_string.length() == 0)
        {
          throw std::runtime_error("StringFunction::derivative: must set gradient strings first");
        }
      std::string fname = getName();
      std::string new_name = "deriv_"+fname;
      //std::cout << "fname= " << fname << " new_name= " << new_name << std::endl;
      //Util::pause(true, "tmp:: StringFunction::derivative");
      // FIXME
      return Teuchos::rcp(new StringFunction(m_gradient_string.c_str(), Name(new_name), Dimensions(3), Dimensions(3) ));
    }

    void StringFunction::operator()(MDArray& inp, MDArray& out, double time_value_optional)
    {
      EXCEPTWATCH;
      argsAreValid(inp, out);

      int domain_rank   = m_domain_dimensions.size();
      int codomain_rank = m_codomain_dimensions.size();

      // FIXME move to argsAreValid
      VERIFY_OP(domain_rank, ==, 1, "StringFunction::operator(): must specify domain Dimensions as rank 1 (for now)");
      VERIFY_OP(codomain_rank, ==, 1, "StringFunction::operator(): must specify codomain Dimensions as rank 1 (for now)");

      int inp_rank      = inp.rank();
      int out_rank      = out.rank();
      int inp_offset    = inp_rank - domain_rank;
      int out_offset    = out_rank - codomain_rank;

      // returns 1,1,1,... etc. if that dimension doesn't exist
      static int n_inp_points[3] = {1,1,1};
      static int n_out_points[3] = {1,1,1};
      first_dimensions(inp, inp_offset, n_inp_points);  
      first_dimensions(out, out_offset, n_out_points);  

      // FIXME check n_out_points, n_inp_points consistency

      int nInDim  = m_domain_dimensions[0];
      int nOutDim = m_codomain_dimensions[0];

      static MDArray inp_loc(3);
      static MDArray out_loc(3);
      if (inp_loc.dimension(0) != nInDim)
        inp_loc.resize( nInDim);
      if (out_loc.dimension(0) != nOutDim)
        out_loc.resize( nOutDim);

      enum { maxRank = 3};
      int iDim[maxRank-1] = {0,0};

      for (iDim[0] = 0; iDim[0] < n_inp_points[0]; iDim[0]++)
        {
          for (iDim[1] = 0; iDim[1] < n_inp_points[1]; iDim[1]++)
            {
              for (int id = 0; id < nInDim; id++)
                {
                  switch (inp_rank)
                    {
                    case 3:
                      inp_loc(id) = inp(iDim[0], iDim[1], id);
                      break;
                    case 2:
                      inp_loc(id) = inp(iDim[0], id);
                      break;
                    case 1:
                      inp_loc(id) = inp( id);
                      break;
                    }
                }

              if (inp_rank == 3)
                {
                  switch (nInDim)
                    {
                    case 3:
                      m_x = inp(iDim[0], iDim[1], 0);
                      m_y = inp(iDim[0], iDim[1], 1);
                      m_z = inp(iDim[0], iDim[1], 2);
                      m_t = time_value_optional;
                      break;
                    case 2:
                      m_x = inp(iDim[0], iDim[1], 0);
                      m_y = inp(iDim[0], iDim[1], 1);
                      m_t = time_value_optional;
                      m_z = 0.0;
                      break;
                    case 1:
                      m_x = inp(iDim[0], iDim[1], 0);
                      m_t = time_value_optional;
                      break;
                    case 0:  // does this make sense?
                      m_t = time_value_optional;
                      break;
                    default:
                      {
                        std::ostringstream msg;
                        msg << "StringFunction::operator() wrong number of dimensions = " << inp.dimension(0);
                        throw std::runtime_error(msg.str());
                      }
                      break;
                    }

                }
              else if (inp_rank == 2)
                {
                  switch (nInDim)
                    {
                    case 3:
                      m_x = inp(iDim[0], 0);
                      m_y = inp(iDim[0], 1);
                      m_z = inp(iDim[0], 2);
                      m_t = time_value_optional;
                      break;
                    case 2:
                      m_x = inp(iDim[0], 0);
                      m_y = inp(iDim[0], 1);
                      m_t = time_value_optional;
                      m_z = 0.0;
                      break;
                    case 1:
                      m_x = inp(iDim[0], 0);
                      m_t = time_value_optional;
                      break;
                    case 0:  // does this make sense?
                      m_t = time_value_optional;
                      break;
                    default:
                      {
                        std::ostringstream msg;
                        msg << "StringFunction::operator() wrong number of dimensions = " << inp.dimension(0);
                        throw std::runtime_error(msg.str());
                      }
                      break;
                    }

                }
              else if (inp_rank == 1)
                {
                  switch (nInDim)
                    {
                    case 3:
                      m_x = inp(0);
                      m_y = inp(1);
                      m_z = inp(2);
                      m_t = time_value_optional;
                      break;
                    case 2:
                      m_x = inp(0);
                      m_y = inp(1);
                      m_t = time_value_optional;
                      m_z = 0.0;
                      break;
                    case 1:
                      m_x = inp(0);
                      m_t = time_value_optional;
                      break;
                    case 0:  // does this make sense?
                      m_t = time_value_optional;
                      break;
                    default:
                      {
                        std::ostringstream msg;
                        msg << "StringFunction::operator() wrong number of dimensions = " << inp.dimension(0);
                        throw std::runtime_error(msg.str());
                      }
                      break;
                    }

                }


              if (m_func_to_value.size())
                {
                  evalFunctions(inp_loc, time_value_optional);
                }

              //Save the evaluations slightly differently whether it's a scalar or vector valued function
              if ((int)m_v.size() != nOutDim)
                {
                  m_v.resize(nOutDim);
                }
//               if ((int)m_g.size() != nOutDim)
//                 {
//                   m_g.resize(nOutDim);
//                 }

              if(nOutDim == 1)
                m_v[0] = m_functionExpr.evaluate();
              else
                m_functionExpr.evaluate();

              for (int iOutDim = 0; iOutDim < nOutDim; iOutDim++)
                {
                  if (out_rank == 1)
                    out(iOutDim) = m_v[iOutDim];
                  else if (out_rank == 2)
                    out(iDim[0], iOutDim) = m_v[iOutDim];
                  else if (out_rank == 3)
                    out(iDim[0], iDim[1], iOutDim) = m_v[iOutDim];
                  else
                    {
                      VERIFY_1("StringFunction::operator() bad output rank");
                    }
                }
            }
        }
    }

    ///  Dimensions of parametric_coordinates and input_phy_points are required to be ([P],[D])
    /// output_values: ([P], [DOF])
    void StringFunction::operator()(MDArray& input_phy_points, MDArray& output_values,
                                    const stk::mesh::Entity& element, const MDArray& parametric_coordinates, double time_value_optional) 
    {
      argsAreValid(input_phy_points, output_values);
      argsAreValid(parametric_coordinates, output_values);

      m_element = &element;
      m_have_element = true;
      VERIFY_OP(parametric_coordinates.rank(), ==, 2, "StringFunction::operator() parametric_coordinates rank bad");
      m_parametric_coordinates = MDArray(1, parametric_coordinates.dimension(1));
      int nPoints=parametric_coordinates.dimension(0);
      int spaceDim = parametric_coordinates.dimension(1);
      for (int iPoint = 0; iPoint < nPoints; iPoint++)
        {
          for (int iSpace=0; iSpace < spaceDim; iSpace++)
            m_parametric_coordinates(0, iSpace)=parametric_coordinates(iPoint, iSpace);
          (*this)(input_phy_points, output_values, time_value_optional);
        }

      // reset this else we won't be able to reuse this object correctly
      m_have_element = false;
      m_element = 0;

      // *  Dimensions of parametric_coordinates are required to be ([P],[D])
      //FieldFunction:: void operator()(const stk::mesh::Entity *element, const MDArray& parametric_coordinates, MDArray& out);

    }

    void StringFunction::operator()(MDArray& input_phy_points, MDArray& output_values,
                                    const stk::mesh::Bucket& bucket, const MDArray& parametric_coordinates, double time_value_optional) 
    {
      VERIFY_OP(input_phy_points.rank(), ==, 3, "StringFunction::operator() must pass in input_phy_points(numCells, numPointsPerCell, spaceDim)");
      int nPoints = input_phy_points.dimension(1);
      int nSpaceDim = input_phy_points.dimension(2);
      int nOutDim = output_values.dimension(2);
      MDArray input_phy_points_one(1, nPoints, nSpaceDim);
      MDArray output_values_one   (1, nPoints, nOutDim);
      const unsigned num_elements_in_bucket = bucket.size();
      for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
        {
          stk::mesh::Entity& element = bucket[iElement];
          for (int iPoint = 0; iPoint<nPoints; iPoint++)
            {
              for (int iSpaceDim=0; iSpaceDim < nSpaceDim; iSpaceDim++)
                {
                  input_phy_points_one(0, iPoint, iSpaceDim) = input_phy_points(iElement, iPoint, iSpaceDim);
                }
            }
          (*this)(input_phy_points_one, output_values_one, element, parametric_coordinates, time_value_optional);
          for (int iPoint = 0; iPoint<nPoints; iPoint++)
            {
              for (int iDOF = 0; iDOF < nOutDim; iDOF++)
                {
                  output_values(iElement, iPoint, iDOF) = output_values_one(0, iPoint, iDOF);
                }
            }
        }
    }

    void
    StringFunction::resolve(stk::expreval::VariableMap::iterator & var_it)
    {
      // /* %TRACE% */ Traceback trace__("sierra::Encr::StringFunction::resolve(stk::expreval::VariableMap::iterator & var_it)"); /* %TRACE% */
      //  Fmwk::Parameters * global_params = Fmwk::Domain::singleton()->parameters().get_nested(FUNCTION_ROOT)->get_nested(GLOBAL_FUNCTION_PARAMETERS);
      //  const Fmwk::Parameters * use_funcs = parameters->get_nested(USE_FUNCTIONS);

      std::string name = (*var_it).first;

      bool use_funcs = true;
      //Not a parameter... so see if it matches the normal stuff
#define DOPRINT 0
      if (DOPRINT) std::cout << " resolve: name= " << name << " " ;
      if (!(name).compare("x"))
        {
          if (DOPRINT) std::cout << " bind to x = " << name;
          (*var_it).second->bind(m_x);
        }
      else if (!(name).compare("y"))
        {
          if (DOPRINT) std::cout << " bind to y = " << name;
          (*var_it).second->bind(m_y);
        }
      else if (!(name).compare("z"))
        {
          if (DOPRINT) std::cout << " bind to z = " << name;

          (*var_it).second->bind(m_z);
        }
      else if (!(name).compare("t"))
        {
          if (DOPRINT) std::cout << " bind to t = " << name;

          (*var_it).second->bind(m_t);
        }
      else if (!(name).compare("v"))
        {
          if (DOPRINT) std::cout << "bind to v = " <<  name << std::endl;
          (*var_it).second->bind(m_v[0]);
        }
//       else if (!(name).compare("g"))
//         {
//           if (DOPRINT) std::cout << "bind to g = " << name << std::endl;
//           (*var_it).second->bind(m_g[0]);
//         }
      //   else if (!(name).compare("g"))
      //     (*var_it).second->bind(g[0]);
      /*
        else if (parameters->exists(name))
        {
        const Real &const_param = parameters->get(name)->value();
        Real &param = const_cast<Real &>(const_param);

        (*var_it).second->bind(param);
        }
        else if(global_params)
        {
        if (global_params->exists(name))
        {
        const Real &const_param = global_params->get(name)->value();
        Real &param = const_cast<Real &>(const_param);

        (*var_it).second->bind(param);
        }
        }
        else if(use_funcs)
        {
        if (use_funcs->exists(name))
        {
        const String & func_name = use_funcs->get(name)->value();

        Function * func = FunctionFactory::instance()->getFunction(region, func_name);
        func_to_value[func] = std::vector<double>(func->getDimension(),0.0);

        //Bind that solution in
        (*var_it).second->bind(func_to_value[func][0]);
        }
        }
      */

      else if(use_funcs)
        {
          if (DOPRINT) std::cout << "bind to function= " << name << std::endl;
          Function * func = Function::getNameToFunctionMap()[name];
          if (DOPRINT) std::cout << "bind to function= " << func << std::endl;
          //const String & func_name = use_funcs->get(name)->value();
          int func_rank = func->getCodomainDimensions().size();
          int func_dim = func->getCodomainDimensions()[func_rank-1];
          if (func)
            {
              //m_func_to_value[func] = std::vector<double>(func->getDomainDimensions()[0], 0.0); // FIXME for multi points
              m_func_to_value[func] = std::vector<double>(func_dim, 0.0); // FIXME for multi points

            }
          else
            {
              std::ostringstream msg;
              msg << "StringFunction::resolve: unknown function name = " << name;
              throw new std::runtime_error(msg.str());
            }

          //Bind that solution in
          (*var_it).second->bind(m_func_to_value[func][0]);
        }

      else
        {
          std::ostringstream msg;
          msg << "Unable to resolve symbols in expression: " << name;
          throw new std::runtime_error(msg.str());
        }
      if (DOPRINT) std::cout << " end resolve: "<<std::endl;;

      return;
    }
  }
}
