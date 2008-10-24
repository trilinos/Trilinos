
#include <fei_iostream.hpp>
#include <fei_ParameterSet.hpp>
#include <Teuchos_ParameterList.hpp>
#include <fei_Trilinos_Helpers.hpp>

#include <fei_unit_Params.hpp>

#include <vector>
#include <cmath>

#undef fei_file
#define fei_file "fei_unit_Params.cpp"
#include "fei_ErrMacros.hpp"

int test_Params_test1()
{
  FEI_COUT << "testing copying fei::ParameterSet to Teuchos::ParameterList...";

  fei::ParameterSet paramset;
  Teuchos::ParameterList paramlist;

  paramset.add(fei::Param("int1", 1));
  paramset.add(fei::Param("double1",1.0));
  paramset.add(fei::Param("str1", "1.0"));

  Trilinos_Helpers::copy_parameterset(paramset, paramlist);

  try {
    int i1 = paramlist.get<int>("int1");
    double d1 = paramlist.get<double>("double1");
    std::string str1 = paramlist.get<std::string>("str1");

    if (i1 != 1 || std::abs(d1-1.0) >= 1.e-14 || str1 != "1.0") {
      FEI_COUT << "FAILED wrong value" << FEI_ENDL;
      return(-1);
    }
  }
  catch (...) {
    FEI_COUT << "FAILED get" << FEI_ENDL;
    return(-1);
  }

  FEI_COUT << "ok" << FEI_ENDL;
  return(0);
}

int test_Params_test2()
{
  FEI_COUT << "testing copying Teuchos::ParameterList to fei::ParameterSet...";

  fei::ParameterSet paramset;
  Teuchos::ParameterList paramlist;

  paramlist.set("int1", 1);
  paramlist.set("double1",1.0);
  paramlist.set("str1", "1.0");

  Trilinos_Helpers::copy_parameterlist(paramlist, paramset);

  int i1;
  int err = paramset.getIntParamValue("int1", i1);
  double d1;
  err += paramset.getDoubleParamValue("double1", d1);
  std::string str1;
  err += paramset.getStringParamValue("str1", str1);

  if (err != 0) {
    FEI_COUT << "FAILED" << FEI_ENDL;
    return(-1);
  }

  FEI_COUT << "ok" << FEI_ENDL;
  return(0);
}

bool test_Params::run(MPI_Comm /*comm*/)
{
  if (test_Params_test1() != 0) {
    throw std::runtime_error("test_Params_test1 failed.");
  }

  if (test_Params_test2() != 0) {
    throw std::runtime_error("test_Params_test2 failed.");
  }

  return true;
}

