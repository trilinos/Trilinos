
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>

#include <fei_iostream.hpp>
#include <fei_ArrayUtils.hpp>

#include <vector>
#include <cmath>

TEUCHOS_UNIT_TEST(fei_utils, insertion_sort_with_companions)
{
  int len = 5;
  std::vector<int> iarray(len);
  std::vector<double> darray(len);

  iarray[0] = 2;
  iarray[1] = 3;
  iarray[2] = 0;
  iarray[3] = 4;
  iarray[4] = 1;

  darray[0] = 2.0;
  darray[1] = 3.0;
  darray[2] = 0.0;
  darray[3] = 4.0;
  darray[4] = 1.0;

  fei::insertion_sort_with_companions(len, &iarray[0], &darray[0]);

  for(int i=0; i<len; ++i) {
    TEUCHOS_TEST_EQUALITY(iarray[i], i, out, success);

    TEUCHOS_TEST_EQUALITY(std::abs(darray[i] - 1.0*i) < 1.e-49, true, out, success);
  }

  iarray.resize(20);

  len = 4;

  iarray[10] = 91;
  iarray[11] = 2225;
  iarray[12] = 214;
  iarray[13] = 3;

  fei::insertion_sort_with_companions(len, &iarray[10], &darray[0]);
}

