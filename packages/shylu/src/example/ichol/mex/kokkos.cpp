#include <Kokkos_Core.hpp>
#include <impl/Kokkos_Timer.hpp>
#include <mex.h>
#include <matrix.h>

using namespace std;

// Testing initialization and finalization of Kokkos

#undef  __FUNCT__ 
#define __FUNCT__ "kokkos"
void mexFunction(int nlhs,
                 mxArray *plhs [],
                 int nrhs,
                 const mxArray *prhs []) {

  Kokkos::Impl::Timer time;

  Kokkos::initialize();
  string name = typeid(Kokkos::DefaultExecutionSpace).name();
  mexPrintf("\n Kokkos is initialized with a default spaceL: %s\n", name.c_str());

  Kokkos::finalize();
  mexPrintf("\n Kokkos is finalized\n");

  plhs[0] = mxCreateDoubleScalar(time.seconds());
}
