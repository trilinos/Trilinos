
#include <fei_iostream.hpp>
#include <fei_ArrayUtils.hpp>

#include <fei_unit_utils.hpp>

#include <vector>
#include <cmath>

void test_insertion_sort_with_companions()
{
  FEI_COUT << "testing fei::insertion_sort_with_companions...";

  int len = 5;
  std::vector<int> array(len);
  std::vector<double> darray(len);

  array[0] = 2;
  array[1] = 3;
  array[2] = 0;
  array[3] = 4;
  array[4] = 1;

  darray[0] = 2.0;
  darray[1] = 3.0;
  darray[2] = 0.0;
  darray[3] = 4.0;
  darray[4] = 1.0;

  fei::insertion_sort_with_companions(len, &array[0], &darray[0]);

  for(int i=0; i<len; ++i) {
    if (array[i] != i) {
      throw std::runtime_error("insertion_sort test 1 failed.");
    }

    if (std::abs(darray[i] - 1.0*i) > 1.e-49) {
      throw std::runtime_error("insertion_sort test 2 failed.");
    }
  }

  array.resize(20);

  len = 4;

  array[10] = 91;
  array[11] = 2225;
  array[12] = 214;
  array[13] = 3;

  fei::insertion_sort_with_companions(len, &array[10], &darray[0]);

  FEI_COUT << "ok"<<FEI_ENDL;

}

bool test_utils::run(MPI_Comm comm)
{
  test_insertion_sort_with_companions();

  return true;
}

