#include <limits>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <sys/time.h>

#include <gtest/gtest.h>
#include <Kokkos_Core.hpp>

#include "StructuredMeshKernel.hpp"

#include <iostream>

TEST(refine, structured_functor)
{
  StructuredMeshKernel kernel;
  kernel.init(1000,4);
  kernel.run_functor();
  kernel.postprocess();
}

TEST(refine, structured_lambda)
{
  StructuredMeshKernel kernel;
  kernel.init(1000,4);
  kernel.run_lambda();
  kernel.postprocess();
}
