
#include <KokkosArray_Host.hpp>
#include <KokkosArray_Cuda.hpp>

#include <TestSpmv.hpp>
#include <TestCG.hpp>

namespace Test {

int test_cuda()
{
  KokkosArray::Cuda::initialize();
  test_spmv_driver<KokkosArray::Cuda>("Cuda");
  test_cgsolve_driver<KokkosArray::Cuda>("Cuda");
  KokkosArray::Cuda::finalize();
  return 0 ;
}

}

