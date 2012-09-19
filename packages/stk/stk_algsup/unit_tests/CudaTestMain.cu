#include <iostream>
#include <cuda.h>

/*--------------------------------------------------------------------*/

static void cuda_safe_call( cudaError e , const char * name )
{
  if ( cudaSuccess != e ) {
    fprintf(stderr,"%s error: %s\n",name,cudaGetErrorString(e) );
    exit( EXIT_FAILURE );
  }
}

#define CUDA_SAFE_CALL( call )  cuda_safe_call( call , # call )

/*--------------------------------------------------------------------*/

// Kernel that executes on the CUDA device  
__global__ void square_array(float *a, int N)  
{  
  int idx = blockIdx.x * blockDim.x + threadIdx.x;  
  if (idx<N) a[idx] = a[idx] * a[idx];  
}  
  

int main(int argc, char** argv)
{
  int deviceCount;
  CUDA_SAFE_CALL( cudaGetDeviceCount(&deviceCount) );
  printf("deviceCount: %d\n", deviceCount);

  if (deviceCount <= 0) {
    printf("Error, no CUDA devices.\n");
    exit( EXIT_FAILURE );
  }

  float *a_h, *a_d;  // Pointer to host & device arrays  
  const int N = 10;  // Number of elements in arrays  
  size_t size = N * sizeof(float);  

  a_h = (float *)malloc(size);        // Allocate array on host  
  CUDA_SAFE_CALL( cudaMalloc((void **) &a_d, size) );   // Allocate array on device  

  // Initialize host array and copy it to CUDA device  
  for (int i=0; i<N; i++) a_h[i] = (float)i;  

  CUDA_SAFE_CALL( cudaMemcpy(a_d, a_h, size, cudaMemcpyHostToDevice) );

  // Do calculation on device:  
  int block_size = 4;  
  int n_blocks = N/block_size + (N%block_size == 0 ? 0:1);  
  square_array <<< n_blocks, block_size >>> (a_d, N);  

  CUDA_SAFE_CALL( cudaThreadSynchronize() );

  // Retrieve result from device and store it in host array  
  CUDA_SAFE_CALL( cudaMemcpy(a_h, a_d, sizeof(float)*N, cudaMemcpyDeviceToHost) );

  // Print results  
  bool test_passed = true;
  for (int i=0; i<N; i++) {
    if (a_h[i] != i*i) test_passed = false;
    std::cout << "array["<<i<<"]: " << a_h[i] << std::endl;
  }

  if (test_passed) {
    std::cout << "Test passed." << std::endl;
  }
  else {
    std::cout << "Test failed." << std::endl;
  }

  // Cleanup  
  free(a_h); CUDA_SAFE_CALL( cudaFree(a_d) );

  return(0);
}

