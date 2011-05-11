# Check for CUDA support
INCLUDE(TPLDeclareLibraries)

SET(_CUDA_FAILURE OFF)

# Have CMake find CUDA
IF(NOT _CUDA_FAILURE)
  FIND_PACKAGE(CUDA 3.2 REQUIRED)
  IF (NOT CUDA_FOUND)
    SET(_CUDA_FAILURE ON)
  ENDIF()
ENDIF()

# # Test that CUDA compiler works
# IF(NOT _CUDA_FAILURE)
#   INCLUDE(TrilinosCUDASupport) 
#   SET(SRC "
#     #include <cuda_runtime.h>
#     __global__ void vecAdd(const float* a, const float* b, float* c, int N)
#     {
#         int i = blockDim.x * blockIdx.x + threadIdx.x;
#         if (i < N) c[i] = a[i] + b[i];
#     }
#     __global__ void vecInit(float* x, float val, int N)
#     {
#         int i = blockDim.x * blockIdx.x + threadIdx.x;
#         if (i < N) x[i] = val;
#     }
#     int main() {
#         const int N               = 2048;
#         const int threadsPerBlock = 256;
#         const int blocksPerGrid   = 8;
#         float* a = NULL;
#         float* b = NULL;
#         float* c = NULL;
#         cudaMalloc((void**)&a, N);
#         cudaMalloc((void**)&b, N);
#         cudaMalloc((void**)&c, N);
#         // init
#         vecInit<<<blocksPerGrid, threadsPerBlock>>>(a,1.0f,N);
#         vecInit<<<blocksPerGrid, threadsPerBlock>>>(b,2.0f,N);
#         vecInit<<<blocksPerGrid, threadsPerBlock>>>(c,0.0f,N);
#         // run
#         vecAdd<<<blocksPerGrid, threadsPerBlock>>>(a, b, c, N);
#     }
#   ")
#   CHECK_CUDA_SOURCE_COMPILES(${SRC} _NVCC_SUCCESS)
#   IF(NOT _NVCC_SUCCESS)
#     SET(_CUDA_FAILURE ON)
#   ENDIF()
# ENDIF()

IF(NOT _CUDA_FAILURE)
  # if we haven't met failure
  macro(PACKAGE_ADD_CUDA_LIBRARY cuda_target)
    PACKAGE_ADD_LIBRARY(${cuda_target} ${ARGN} CUDALIBRARY)
  endmacro()
  GLOBAL_SET(TPL_CUDA_LIBRARY_DIRS)
  GLOBAL_SET(TPL_CUDA_INCLUDE_DIRS ${CUDA_TOOLKIT_INCLUDE})
  GLOBAL_SET(TPL_CUDA_LIBRARIES ${CUDA_CUDART_LIBRARY} ${CUDA_cublas_LIBRARY} ${CUDA_cufft_LIBRARY})
ELSE()
  SET(TPL_ENABLE_CUDA OFF PARENT_SCOPE)
  MESSAGE(FATAL_ERROR "\nDid not find acceptable version of CUDA compiler")
ENDIF()
