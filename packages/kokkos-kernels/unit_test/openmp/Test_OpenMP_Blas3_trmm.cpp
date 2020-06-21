#include<Test_OpenMP.hpp>
// Remove this ifdef once we have a fall back implementation.
#ifdef KOKKOSKERNELS_ENABLE_TPL_BLAS
#include<Test_Blas3_trmm.hpp>
#endif // KOKKOSKERNELS_ENABLE_TPL_BLAS
