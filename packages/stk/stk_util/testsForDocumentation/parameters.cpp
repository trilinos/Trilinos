#include <gtest/gtest.h>
#include <complex>
#include <vector>
#include <string>
#include <mpi.h>
#include <stk_util/util/ParameterList.hpp>

namespace
{
  TEST(StkUtilTestForDocumentation, Parameters)
  {
    //-BEGIN-init
    //+ INITIALIZATION
    std::vector<std::string> exp_name;
    std::vector<stk::util::ParameterType::Type> exp_type;
    
    double pi = 3.14159;
    float e = 2.71828;
    int answer = 42;
    int64_t big_answer = 42000000000001;
    std::string team_name = "STK Transition Team";

    exp_name.push_back("PI"); exp_type.push_back(stk::util::ParameterType::DOUBLE);
    exp_name.push_back("E"); exp_type.push_back(stk::util::ParameterType::FLOAT);
    exp_name.push_back("Answer"); exp_type.push_back(stk::util::ParameterType::INTEGER);
    exp_name.push_back("Answer_64"); exp_type.push_back(stk::util::ParameterType::INT64);
    exp_name.push_back("TeamName"); exp_type.push_back(stk::util::ParameterType::STRING);

    std::vector<double> my_double_vector;
    my_double_vector.push_back(2.78); my_double_vector.push_back(5.30);
    my_double_vector.push_back(6.21);
    exp_name.push_back("some_doubles"); exp_type.push_back(stk::util::ParameterType::DOUBLEVECTOR);
    
    std::vector<float> my_float_vector;
    my_float_vector.push_back(194.0); my_float_vector.push_back(-194.0);
    my_float_vector.push_back(47.0);  my_float_vector.push_back(92.0);
    exp_name.push_back("some_floats"); exp_type.push_back(stk::util::ParameterType::FLOATVECTOR);
    
    std::vector<int> ages;
    ages.push_back(55); ages.push_back(49); ages.push_back(21); ages.push_back(19);
    exp_name.push_back("Ages"); exp_type.push_back(stk::util::ParameterType::INTEGERVECTOR);
    
    std::vector<int64_t> ages_64;
    ages_64.push_back(55); ages_64.push_back(49); ages_64.push_back(21); ages_64.push_back(19);
    exp_name.push_back("Ages_64"); exp_type.push_back(stk::util::ParameterType::INT64VECTOR);
    
    std::vector<std::string> names;
    names.push_back("greg"); names.push_back("chloe"); names.push_back("tuffy");
    names.push_back("liberty"); names.push_back("I have spaces");
    exp_name.push_back("Names"); exp_type.push_back(stk::util::ParameterType::STRINGVECTOR);
    //-END-init
    
    //-BEGIN-define
    //+ Define parameters...
    stk::util::ParameterList params;
    params.set_param("PI", pi);
    params.set_param("E", e);
    params.set_param("Answer", answer);
    params.set_param("Answer_64", big_answer);
    params.set_param("TeamName", team_name);
    params.set_param("some_doubles", my_double_vector);
    params.set_param("some_floats", my_float_vector); 
    params.set_param("Ages", ages); 
    params.set_param("Ages_64", ages_64); 
    params.set_param("Names", names); 
    //-END-define

    //-BEGIN-access
    //+ Write parameters to stdout...
    params.write_parameter_list(std::cout);
      
    //+ Access parameters by name...
    size_t num_param = exp_name.size();
    for (size_t i=0; i < num_param; i++) {
      stk::util::Parameter &param = params.get_param(exp_name[i]);
      EXPECT_EQ(param.type, exp_type[i]);
    }

    //+ Extract some parameter values if know type:
    std::vector<int> pages = params.get_value<std::vector<int> >("Ages");
    for (size_t i=0; i < pages.size(); i++) {
      EXPECT_EQ(pages[i], ages[i]);
    }

    double my_pi = params.get_value<double>("PI");
    EXPECT_EQ(my_pi, pi);

    //+ Change value of an existing parameter
    params.set_value("Answer", 21);

    int new_answer = params.get_value<int>("Answer");
    EXPECT_EQ(new_answer, 21);
    
    //-END-access

    //-BEGIN-error
    //+ If the requested parameter does not exist, 
    //+ an error message is printed to stderr and an invalid
    //+ parameter object is returned
    stk::util::Parameter no_exist = params.get_param("DoesNotExist");
    EXPECT_EQ(no_exist.type, stk::util::ParameterType::INVALID);
    
    //+ In this method of requesting a parameter, no error
    //+ message is printed if the parameter doesn't exist and
    //+ instead the returned iterator is equal to the end of the
    //+ parameter list.
    stk::util::ParameterMapType::iterator it = params.find("DoesNotExist");
    EXPECT_TRUE(it == params.end());
    
    //+ If the parameter types do not match, an error message is
    //+ printed and the value 0 of the requested type is returned.
    int invalid = params.get_value<int>("PI");
    EXPECT_EQ(invalid, 0);
    
    //+ If the parameter types do not match, an error message is
    //+ printed and an empty vector of the requested type is returned.
    std::vector<double> pies = params.get_value<std::vector<double> >("PI");
    EXPECT_EQ(pies.size(), 0u);
    //-END-error

    //-BEGIN-usertype
    //+ Adding a parameter of "unsupported" type...
    stk::util::ParameterList more_params;
    std::complex<double> phase(3.14,2.718);
    more_params.set_param("phase", phase);

    //+ The print system doesn't know about this type, so will print
    //+ a warning message about unrecognized type.
    more_params.write_parameter_list(std::cout);

    //+ However, you can still retrieve the value of the parameter
    //+ if you know what type it is.
    std::complex<double> my_phase = more_params.get_value<std::complex<double> >("phase");
    EXPECT_EQ(my_phase, phase);

    //+ The Parameter class won't help you on determining the type,
    //+ You must know what it is.
    EXPECT_EQ(more_params.get_param("phase").type, stk::util::ParameterType::INVALID);

    //+ If the wrong type is specified, an exception will be thrown...
    EXPECT_THROW(more_params.get_value<std::complex<int> >("phase"), std::exception);
    //-END-usertype
  }
}

// type not supported...
// iterator access valid
