
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_iostream.hpp>
#include <fei_ParameterSet.hpp>
#include <Teuchos_ParameterList.hpp>
#include <fei_Trilinos_Helpers.hpp>

#include <vector>
#include <cmath>

#undef fei_file
#define fei_file "fei_unit_Params.cpp"
#include "fei_ErrMacros.hpp"

TEUCHOS_UNIT_TEST(Params, Params_test1)
{
  fei::ParameterSet paramset;
  Teuchos::ParameterList paramlist;

  paramset.add(fei::Param("int1", 1));
  paramset.add(fei::Param("double1",1.0));
  paramset.add(fei::Param("str1", "1.0"));

  Trilinos_Helpers::copy_parameterset(paramset, paramlist);

  TEUCHOS_TEST_EQUALITY(paramlist.get<int>("int1"), 1, out, success);
  TEUCHOS_TEST_EQUALITY((paramlist.get<double>("double1") -1.0) < 1.e-14, true, out, success);
  TEUCHOS_TEST_EQUALITY(paramlist.get<std::string>("str1"), std::string("1.0"), out, success);
}

TEUCHOS_UNIT_TEST(Params, Params_test2)
{
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

  TEUCHOS_TEST_EQUALITY(err, 0, out, success);
}

