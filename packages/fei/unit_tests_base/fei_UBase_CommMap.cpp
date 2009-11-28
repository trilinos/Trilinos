
#include <Teuchos_ConfigDefs.hpp>
#include <Teuchos_UnitTestHarness.hpp>
#include <fei_CommMap.hpp>
#include <iostream>

namespace {

TEUCHOS_UNIT_TEST_TEMPLATE_1_DECL(CommMap, test0, T)
{
  typename fei::CommMap<T>::Type comm_map;

  std::vector<T> input_items(4);
  input_items[0] = 2;
  input_items[1] = 3;
  input_items[2] = 1;
  input_items[3] = 0;

  fei::addItemsToCommMap(0, input_items.size(), &input_items[0], comm_map);
  fei::addItemsToCommMap(1, input_items.size(), &input_items[0], comm_map, false);

  std::vector<T>& sorted_items = comm_map[0];
  std::vector<T>& unsorted_items = comm_map[1];

  std::vector<T> expected_sorted_items(4);
  expected_sorted_items[0] = 0;
  expected_sorted_items[1] = 1;
  expected_sorted_items[2] = 2;
  expected_sorted_items[3] = 3;

  bool sorted_correct = sorted_items == expected_sorted_items;
  bool unsorted_correct = unsorted_items == input_items;

  TEUCHOS_TEST_EQUALITY(sorted_correct, true, out, success);
  TEUCHOS_TEST_EQUALITY(unsorted_correct, true, out, success);
}

#define UNIT_TEST_GROUP(TYPE) \
  TEUCHOS_UNIT_TEST_TEMPLATE_1_INSTANT(CommMap,test0,TYPE)

typedef long int longint;
UNIT_TEST_GROUP(int)
UNIT_TEST_GROUP(longint)
UNIT_TEST_GROUP(float)
UNIT_TEST_GROUP(double)

}//namespace <anonymous>

