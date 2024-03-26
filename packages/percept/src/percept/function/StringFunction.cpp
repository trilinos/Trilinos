// Copyright 2002 - 2008, 2010, 2011 National Technology Engineering
// Solutions of Sandia, LLC (NTESS). Under the terms of Contract
// DE-NA0003525 with NTESS, the U.S. Government retains certain rights
// in this software.
//
// Use of this source code is governed by a BSD-style
// license that can be found in the LICENSE file.

#define DOPRINT 0

#if DOPRINT
#define PRINTQ(a) do { std::cout << QUOTE(a) << a << std::endl; } while(0)
#define PRINT(a) do { std::cout << a << std::endl; } while(0)
#else
#define PRINTQ(a) do { } while(0)
#define PRINT(a) do {  } while(0)
#endif

#include <string>
#include <sstream>
#include <stdexcept>
#include <sstream>
#include <percept/function/StringFunction.hpp>
#include <percept/ExceptionWatch.hpp>
#include <stk_mesh/base/Bucket.hpp>

  namespace percept
  {
    StringFunction::StringFunction(const char *function_string,
                                   Name name,
                                   int domain_dimension,
                                   int codomain_dimension,
                                   unsigned integration_order) :
      Function(name.getName().c_str(), Dimensions(domain_dimension), Dimensions(codomain_dimension), integration_order),
      m_func_string(function_string), m_functionExpr(*this), m_gradient_string(""), m_spatialDim(0), current_time(0.0), using_current_time(false)
    {
      init();
    }

    StringFunction::StringFunction(const char *function_string,
                                   Name name,
                                   Dimensions domain_dimensions,
                                   Dimensions codomain_dimensions,
                                   unsigned integration_order) :
      Function(name.getName().c_str(), domain_dimensions, codomain_dimensions, integration_order),
      m_func_string(function_string), m_functionExpr(*this), m_gradient_string("") , m_spatialDim(0), current_time(0.0), using_current_time(false)
    {
      init();
    }

    StringFunction::StringFunction(const StringFunction& s) :
      Function(s.m_name.c_str(), s.m_domain_dimensions, s.m_codomain_dimensions, s.m_integration_order),
      m_func_string(s.m_func_string), m_functionExpr(*this), m_gradient_string(s.m_gradient_string), m_spatialDim(s.m_spatialDim), current_time(0.0), using_current_time(false)
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
      m_element = stk::mesh::Entity();
      m_bucket=0;
      m_parametric_coordinates = MDArray("m_parametric_coordinates",1,3);
      m_have_element = false;
      m_have_bucket = false;
    }

    void StringFunction::set_gradient_strings(std::string gstring[3], int len)
    {
      m_gradient_string = "";
      m_spatialDim = len;
      for (int i = 0; i < len; i++)
        m_gradient_string += "v["+std::to_string(i)+"]= "+gstring[i]+";";
    }

    void
    StringFunction::set_gradient_strings(MDArrayString& gstring)
    {
      if (gstring.rank() != 1) throw std::runtime_error("set_gradient_strings takes a rank 1 matrix (i.e. a vector) of strings (MDArrayString)");
      int len = gstring.extent_int(0);
      m_gradient_string = "";
      m_spatialDim = len;
      for (int i = 0; i < len; i++)
        m_gradient_string += "v["+std::to_string(i)+"]= "+gstring(i)+";";
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
      std::map<Function *, std::vector<double> >::iterator it = m_func_to_value.begin();
      const std::map<Function *, std::vector<double> >::iterator it_end = m_func_to_value.end();

      //Evaluate the function and store the answer in its bound location
      for(; it != it_end; ++it)
        {
          Function& func = *it->first;

          PRINTQ(func.getName());

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

          MDArray f_inp = MDArray("f_inp", numCells, func.getNewDomain().extent_int(0));
          MDArray f_out = MDArray("f_out", numCells, func.getNewCodomain().extent_int(0));

          int nInDim = last_dimension(f_inp);
          // FIXME - do we really need the extra array?
          for (int iDim = 0; iDim < nInDim; iDim++)
            {
              f_inp(0, iDim) = inp( iDim);
            }

          if (m_have_element)
            {
              func(f_inp, f_out, m_element, m_parametric_coordinates, time_value_optional);
              PRINT("func() getName()= " << func.getName() << " f_inp=\n" << f_inp << " f_out=\n " << f_out);
            }
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

    static inline void first_dimensions(MDArray& arr, int arr_offset, int *n_points, int max_rank=3)
    {
      for (int ii = 0; ii < max_rank; ii++)
        {
          n_points[ii] = 1;
        }
      for (int ii = 0; ii < arr_offset; ii++)
        {
          n_points[ii] = arr.extent_int(ii);
        }
    }

    // replace all occurrences of input_string with output_string in source_string

    Teuchos::RCP<Function > StringFunction::derivative_test(MDArrayString& deriv_spec)
    {
      bool debug=true;
      std::string fstr = m_func_string;
      // FIXME this is just for a simple test and only works for linear functions
      int outputDim = deriv_spec.extent_int(0);
      static std::string s_xyzt[] = {"x", "y", "z", "t"};
      std::string out_str;
      //if (debug) std::cout << "deriv_spec= " << deriv_spec << std::endl;
      for (int iresult = 0; iresult < deriv_spec.extent_int(0); iresult++)
        {
          fstr = m_func_string;
          for (int jderiv = 0; jderiv < deriv_spec.extent_int(1); jderiv++)
            {
              if (debug) std::cout << "deriv_spec= " << deriv_spec(iresult, jderiv)  << " fstr= " << fstr << std::endl;
              char xyzt = deriv_spec(iresult, jderiv)[0];
              std::replace( fstr.begin(), fstr.end(), xyzt, '1' );
              if (debug) std::cout << "xyzt= " << xyzt << " fstr= " << fstr << std::endl;
            }
          for (int jderiv = 0; jderiv < 4; jderiv++)
            {
              char xyzt = s_xyzt[jderiv][0];
              std::replace( fstr.begin(), fstr.end(), xyzt, '0' );
              if (debug) std::cout << "xyzt= " << xyzt << " fstr= " << fstr << std::endl;
            }
          if (deriv_spec.extent_int(0) > 1)
            {
              out_str += "v["+std::to_string(iresult)+"]= " + fstr+";";
            }
          else
            {
              out_str = fstr;
            }
        }
      std::string fname = getName();
      std::string new_name = "deriv_"+fname;
      if (debug) std::cout << "fname= " << fname << " new_name= " << new_name << " out_str= " << out_str << std::endl;
      //Util::pause(true, "tmp:: StringFunction::derivative");
      return Teuchos::rcp(new StringFunction(out_str.c_str(), Name(new_name), Dimensions(3), Dimensions(outputDim) ));
    }


    Teuchos::RCP<Function > StringFunction::derivative_test_fd(MDArrayString& deriv_spec, double eps)
    {
      std::string eps_string = std::to_string(eps);
      std::string fstr = m_func_string;
      std::string fstr_p = m_func_string;
      std::string fstr_m = m_func_string;
      // FIXME this is just for a simple test and only works for linear functions
      int outputDim = deriv_spec.extent_int(0);
      //static std::string s_xyzt[] = {"x", "y", "z", "t"};
      std::string out_str;
      for (int iresult = 0; iresult < deriv_spec.extent_int(0); iresult++)
        {
          fstr = m_func_string;
          fstr_p = m_func_string;
          fstr_m = m_func_string;
          for (int jderiv = 0; jderiv < deriv_spec.extent_int(1); jderiv++)
            {
              //std::cout << "deriv_spec= " << deriv_spec(iresult, jderiv)  << " fstr= " << fstr << std::endl;
              std::string rep = "("+deriv_spec(iresult, jderiv)+"+"+eps_string+")";
              Util::replace(fstr_p, deriv_spec(iresult, jderiv), rep);
              rep = "("+deriv_spec(iresult, jderiv)+"-"+eps_string+")";
              Util::replace(fstr_m, deriv_spec(iresult, jderiv), rep);
              fstr = "(("+fstr_p+")-("+fstr_m+"))/(2.0*"+eps_string+")";
              //std::cout << "xyzt= " << xyzt << " fstr= " << fstr << std::endl;
            }
          if (deriv_spec.extent_int(0) > 1)
            {
              out_str += "v["+std::to_string(iresult)+"]= " + fstr+";";
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

      return Teuchos::rcp(new StringFunction(m_gradient_string.c_str(), Name(new_name), Dimensions(3), Dimensions(3) ));
    }

    Teuchos::RCP<Function > StringFunction::gradient(int spatialDim)
    {
      if (m_gradient_string.length() == 0)
        {
          throw std::runtime_error("StringFunction::gradient: must set gradient strings first");
        }
      std::string xyz[] = {"x", "y", "z"};
      MDArrayString mda(spatialDim);
      for (int i = 0; i < spatialDim; i++)
        {
          mda(i) = xyz[i];
        }
      return derivative(mda);

    }

    void StringFunction::operator()(MDArray& inp, MDArray& out, double time_value_optional)
    {
      EXCEPTWATCH;
      argsAreValid(inp, out);
      //PRINT("tmp srk StringFunction::operator() getName()= " << getName() << " inp= " << inp << " out= " << out);

      int domain_rank   = m_domain_dimensions.size();
      int codomain_rank = m_codomain_dimensions.size();

      // FIXME move to argsAreValid
      VERIFY_OP(domain_rank, ==, 1, "StringFunction::operator(): must specify domain Dimensions as rank 1 (for now)");
      VERIFY_OP(codomain_rank, ==, 1, "StringFunction::operator(): must specify codomain Dimensions as rank 1 (for now)");

      int inp_rank      = inp.rank();
      int out_rank      = out.rank();
      int inp_offset    = inp_rank - domain_rank;
      int out_offset    = out_rank - codomain_rank;

      // returns 1,1,1,... etc. if that extent_int doesn't exist
      enum { maxRank = 3};
      VERIFY_OP_ON(inp_rank, <=,  maxRank, "StringFunction::operator() input array rank too large");
      VERIFY_OP_ON(out_rank, <=,  maxRank, "StringFunction::operator() output array rank too large");
      int n_inp_points[maxRank] = {1,1,1};
      int n_out_points[maxRank] = {1,1,1};
      first_dimensions(inp, inp_offset, n_inp_points, maxRank);
      first_dimensions(out, out_offset, n_out_points, maxRank);

      int nInDim  = m_domain_dimensions[0];
      int nOutDim = m_codomain_dimensions[0];

      MDArray inp_loc("inp_loc", nInDim);
      MDArray out_loc("out_loc", nOutDim);

      int iDim[maxRank-1] = {0,0};
      double *xyzt[] = {&m_x, &m_y, &m_z, &m_t};

      PRINT("getName()= " << getName() << " n_inp_points= " << n_inp_points[0] << " " << n_inp_points[1] << " " << n_inp_points[2] << " nInDim= " << nInDim);

      if (using_current_time)
	m_t = current_time;
      else
	m_t = time_value_optional;

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
                      *xyzt[id] = inp(iDim[0], iDim[1], id);
                      break;
                    case 2:
                      inp_loc(id) = inp(iDim[0], id);
                      *xyzt[id] = inp(iDim[0], id);
                      break;
                    case 1:
                      inp_loc(id) = inp( id);
                      *xyzt[id] = inp( id);
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

              PRINT("getName()= " << getName() << " about to evaluate m_functionExpr... iDim[0] = " << iDim[0] << " iDim[1]= " << iDim[1]);
              if(nOutDim == 1)
                m_v[0] = m_functionExpr.evaluate();
              else
                m_functionExpr.evaluate();
              PRINT("getName()= " << getName() << " done to evaluate m_functionExpr... iDim[0] = " << iDim[0] << " iDim[1]= " << iDim[1] << " nOutDim= " << nOutDim);
              for (int iOutDim = 0; iOutDim < nOutDim; iOutDim++)
                {
                  switch (out_rank)
                    {
                    case 1: out(iOutDim)                   = m_v[iOutDim]; break;
                    case 2: out(iDim[0], iOutDim)          = m_v[iOutDim]; break;
                    case 3: out(iDim[0], iDim[1], iOutDim) = m_v[iOutDim]; break;
                    default:
                      VERIFY_1("StringFunction::operator() bad output rank");
                    }
                }
            }
        }
      PRINT("tmp srk StringFunction::operator() getName()= " << getName() << " inp=\n " << inp << " out=\n " << out);
      //PRINT("getName()= " << getName() << " done in operator(), out= " << out);
    }

    ///  Dimensions of parametric_coordinates and input_phy_points are required to be ([C],[P],[D]) or ([P],[D])
    /// output_values: ([C], [P], [DOF]), or ( [P], [DOF])
    void StringFunction::operator()(MDArray& input_phy_points, MDArray& output_values,
                                    const stk::mesh::Entity element, const MDArray& parametric_coordinates, double time_value_optional)
    {
      //PRINT("tmp srk StringFunction::operator(element) getName()= " << getName() << " input_phy_points= " << input_phy_points << " output_values= " << output_values);

      argsAreValid(input_phy_points, output_values);
      argsAreValid(parametric_coordinates, output_values);

      m_element = element;
      m_have_element = true;
      VERIFY_OP(parametric_coordinates.rank(), ==, 2, "StringFunction::operator() parametric_coordinates rank bad");
      m_parametric_coordinates = MDArray("m_parametric_coordinates", 1, parametric_coordinates.extent_int(1));
      int nPoints=parametric_coordinates.extent_int(0);
      int spaceDim = parametric_coordinates.extent_int(1);

      int rank_inp = input_phy_points.rank();
      int rank_out = output_values.rank();
      VERIFY_OP_ON(rank_inp, <=, 3, "bad input rank");
      VERIFY_OP_ON(rank_out, <=, 3, "bad output rank");
      int nOutD = output_values.extent_int(1);
      MDArray input_phy_points_one("input_phy_points_one", 1, spaceDim), output_values_one("output_values_one", 1, nOutD);

      int nCells = 1;
      if (rank_inp == 3)
        {
          nCells = input_phy_points.extent_int(0);
          VERIFY_OP_ON(input_phy_points.extent_int(1), ==, nPoints, "bad inp dim");
          VERIFY_OP_ON(input_phy_points.extent_int(2), ==, spaceDim, "bad inp dim");

          VERIFY_OP_ON(output_values.extent_int(0), ==, nCells, "bad out dim");
          VERIFY_OP_ON(output_values.extent_int(1), ==, nPoints, "bad out dim");
          nOutD = output_values.extent_int(2);
          Kokkos::resize(output_values_one,1,nOutD);
        }
      else
        {
          VERIFY_OP_ON(input_phy_points.extent_int(0), ==, nPoints, "bad inp dim");
          VERIFY_OP_ON(input_phy_points.extent_int(1), ==, spaceDim, "bad inp dim");

          VERIFY_OP_ON(output_values.extent_int(0), ==, nPoints, "bad out dim");
        }

      for (int iCell = 0; iCell < nCells; iCell++)
        {
          for (int iPoint = 0; iPoint < nPoints; iPoint++)
            {
              for (int iSpace=0; iSpace < spaceDim; iSpace++)
                {
                  m_parametric_coordinates(0, iSpace) = parametric_coordinates(iPoint, iSpace);
                  switch(rank_inp)
                    {
                    case 2: input_phy_points_one(0, iSpace) = input_phy_points(iPoint, iSpace); break;
                    case 3: input_phy_points_one(0, iSpace) = input_phy_points(iCell, iPoint, iSpace); break;
                    }
                }
              (*this)(input_phy_points_one, output_values_one, time_value_optional);
              for (int iD = 0; iD < nOutD; iD++)
                {
                  //output_values(iPoint, iD) = output_values_one(0, iD);
                  switch(rank_out)
                    {
                    case 2: output_values(iPoint, iD) = output_values_one(0, iD); break;
                    case 3: output_values(iCell, iPoint, iD) = output_values_one(0, iD); break;
                    }
                }
            }
        }

      // reset this else we won't be able to reuse this object correctly
      m_have_element = false;
      m_element = stk::mesh::Entity();

      // *  Dimensions of parametric_coordinates are required to be ([P],[D])
      //FieldFunction:: void operator()(const stk::mesh::Entity element, const MDArray& parametric_coordinates, MDArray& out);
      PRINT("tmp srk StringFunction::operator(element) getName()= " << getName() << " input_phy_points=\n " << input_phy_points << " output_values=\n " << output_values);

    }

    void StringFunction::operator()(MDArray& input_phy_points, MDArray& output_values,
                                    const stk::mesh::Bucket& bucket, const MDArray& parametric_coordinates, double time_value_optional)
    {
      PRINT("tmp srk StringFunction::operator(bucket) getName()= " << getName() << " input_phy_points= " << input_phy_points << " output_values= " << output_values);

      VERIFY_OP(input_phy_points.rank(), ==, 3, "StringFunction::operator() must pass in input_phy_points(numCells, numPointsPerCell, spaceDim)");
      int nPoints = input_phy_points.extent_int(1);
      int nSpaceDim = input_phy_points.extent_int(2);
      int nOutDim = output_values.extent_int(2);
      MDArray input_phy_points_one("input_phy_points_one", 1, nPoints, nSpaceDim);
      MDArray output_values_one   ("output_values_one", 1, nPoints, nOutDim);
      const unsigned num_elements_in_bucket = bucket.size();
      for (unsigned iElement = 0; iElement < num_elements_in_bucket; iElement++)
        {
          stk::mesh::Entity element = bucket[iElement];
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
      //  Fmwk::Parameters * global_params = sierra::Domain::singleton()->parameters().get_nested(FUNCTION_ROOT)->get_nested(GLOBAL_FUNCTION_PARAMETERS);
      //  const Fmwk::Parameters * use_funcs = parameters->get_nested(USE_FUNCTIONS);

      std::string name = (*var_it).first;

      bool use_funcs = true;
      //Not a parameter... so see if it matches the normal stuff
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
      else if(use_funcs)
        {
          if (DOPRINT) std::cout << "bind to function= " << name << std::endl;
          Function * func = Function::getNameToFunctionMap()[name];
          if (DOPRINT) std::cout << "bind to function= " << func << std::endl;
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
      if (DOPRINT) std::cout << " end resolve: "<<std::endl;

      return;
    }
  }
