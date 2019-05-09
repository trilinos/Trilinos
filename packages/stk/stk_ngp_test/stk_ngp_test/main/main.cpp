#include <stk_ngp_test/ngp_test.hpp>

int main(int argc, char** argv) {
  ngp_testing::NgpTestEnvironment testEnv(&argc, argv);
  int returnVal = testEnv.run_all_tests();
  return returnVal;
}
