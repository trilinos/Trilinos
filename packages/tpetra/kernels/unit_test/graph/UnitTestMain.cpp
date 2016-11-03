#include "gtest/gtest.h"


char *input_filename = NULL;
int main(int argc, char **argv) {
  ::testing::InitGoogleTest(&argc, argv);

  for ( int i = 1 ; i < argc ; ++i ) {
    if ( 0 == strcasecmp( argv[i] , "--mtx-file" ) ) {
      input_filename = argv[++i];
    }
  }

  return RUN_ALL_TESTS();
}
