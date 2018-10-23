#ifndef ERRORBDDC_H
#define ERRORBDDC_H

#include <stdio.h>
#include <string.h>

namespace bddc {

#define BDDC_TEST_FOR_EXCEPTION(ierr, x, msg)                           \
    if ((ierr) != 0) {                                                  \
      fprintf(stderr, ">> Error in file %s, line %d, error %d \n",__FILE__,__LINE__,ierr); \
      fprintf(stderr, "   %s\n", msg);                                  \
      throw x(msg);                                                     \
    }

}

#endif // ERRORBDDC_H
