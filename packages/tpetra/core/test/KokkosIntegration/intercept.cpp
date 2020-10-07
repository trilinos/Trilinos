/*
 * a cuda intercept for kokkos deep copies, now with counts
 */

#include <string>
#include <map>
#include <cuda.h>
#include <stdio.h>
#include <dlfcn.h>
#include <mpi.h>
#include <Kokkos_Core.hpp>
#include <Kokkos_View.hpp>
#include <Kokkos_Atomic.hpp>

namespace ApiTest {
  enum devCalls {
    cudaDeviceSynchronize,
    cudaMemcpy2DAsync,
    cudaMemcpy3DAsync,
    cudaMemcpyAsync
  };

  std::map<std::string, int> totalCounter;
  const int numCounters = 4;

  __device__ int *cudaCounters;
  int *hostCounters;
}

int pincr(std::map<std::string, int> &store, std::string key, int val) {
  std::map<std::string, int>::iterator it = store.find(key);
  if (it == store.end()) {
    store.insert(std::make_pair(key, 1));
  } else {
    it->second += val;
  }
  return 0;
}

int incr(std::map<std::string, int> &store, std::string key) {
  return pincr(store, key, 1);
}

void printAll(std::map<std::string, int> &store) {
  std::map<std::string, int>::iterator it = store.begin();
  fprintf (stderr, "**************** cuda call analysis ****************\n"); 
  fprintf (stderr, "call\t\ttimes\n");
  while (it != store.end()) {
    fprintf(stderr, "%s:\t%d\n", it->first.c_str(), it->second);
    it++;
  }
  fprintf (stderr, "**************** cuda call analysis ****************\n"); 
}

int MPI_Init(int *argc, char ***argv) {
  int (*o_mpi_init)(int *, char ***);
  o_mpi_init = (int (*)(int *, char ***))dlsym(RTLD_NEXT, "MPI_Init");

  ApiTest::hostCounters = (int *)malloc (sizeof(int)*4);
  for (size_t i = 0; i < 4; i++) {
    ApiTest::hostCounters[i] = 0;
  }

  return o_mpi_init(argc, argv);
}

int MPI_Finalize(void) {
  fprintf(stderr, "MPI_Finalize()\n");
  int (*o_mpi_finalize)(void);
  std::map<std::string, int>::iterator f;
  o_mpi_finalize = (int (*)(void))dlsym(RTLD_NEXT, "MPI_Finalize");

  //copies out the device calls
  pincr(ApiTest::totalCounter, "cudaDeviceSynchronize", 
	ApiTest::hostCounters[ApiTest::devCalls::cudaDeviceSynchronize]);
  pincr(ApiTest::totalCounter, "cudaMemcpy2DAsync",
	ApiTest::hostCounters[ApiTest::devCalls::cudaMemcpy2DAsync]);
  pincr(ApiTest::totalCounter, "cudaMemcpy3DAsync",
	ApiTest::hostCounters[ApiTest::devCalls::cudaMemcpy3DAsync]);
  pincr(ApiTest::totalCounter, "cudaMemcpyAsync",
	ApiTest::hostCounters[ApiTest::devCalls::cudaMemcpyAsync]);

  //baseline: these numbers are magic, happen with kokkos initialization
  //to update, run KokkosDeepcopy with no deep copies performed and comment
  //this section out.
  pincr(ApiTest::totalCounter, "cudaDeviceSynchronize", -5);
  pincr(ApiTest::totalCounter, "cudaMemcpy", -11);
  pincr(ApiTest::totalCounter, "cudaMemcpyToSymbol", -1);
  pincr(ApiTest::totalCounter, "cudaMemcpy2DAsync", -1);
  pincr(ApiTest::totalCounter, "cudaMemcpy3DAsync", -1);
  pincr(ApiTest::totalCounter, "cudaMemcpyAsync", -1);

  free(ApiTest::hostCounters);

  printAll(ApiTest::totalCounter);

  return o_mpi_finalize();
}

namespace Kokkos {
void initialize(int& narg, char* arg[]) {
  void (*o_init)(int&, char **);
  o_init = (void (*)(int&, char **))dlsym(RTLD_NEXT, "_ZN6Kokkos10initializeERiPPc");

  o_init(narg, arg);

  ApiTest::cudaCounters = (int *)Kokkos::kokkos_malloc<Kokkos::CudaSpace>("cuda calls", sizeof(int)*4);
} 

void finalize() {
  fprintf(stderr, "Kokkos::finalize()\n");
  void (*o_finalize)(void);
  o_finalize = (void (*)(void))dlsym(RTLD_NEXT, "_ZN6Kokkos8finalizeEv");

  //copy ApiTest::cudaCounters into ApiTest::hostCounters
  {
    // This scoping is important since it ensures the Views deallocate *before* the real Kokkos::finalize() is called.
    Kokkos::View<int*,Kokkos::HostSpace> h_cudaCounters("host counter temp",ApiTest::numCounters);
    Kokkos::View<int*,Kokkos::CudaSpace> d_cudaCounters(ApiTest::cudaCounters,ApiTest::numCounters);
    Kokkos::deep_copy(h_cudaCounters,d_cudaCounters);

    for(int i=0; i<ApiTest::numCounters; i++)
      ApiTest::hostCounters[i] += h_cudaCounters[i];
    

    /*printf("CMS: Total calls now look like: ");
    for(int i=0; i<ApiTest::numCounters; i++)
      printf("%d ",ApiTest::hostCounters[i]);
    printf("\n");
    */

    Kokkos::kokkos_free<Kokkos::CudaSpace>(ApiTest::cudaCounters);
  }
  o_finalize();
}

};

__host__ __device__ cudaError_t cudaDeviceSynchronize() {
  cudaError_t (*o_cudaDeviceSynchronize)();
  o_cudaDeviceSynchronize = (cudaError_t (*)())dlsym(RTLD_NEXT, "cudaDeviceSynchronize");
#ifdef __CUDA_ARCH__
  Kokkos::atomic_increment(&ApiTest::cudaCounters[ApiTest::devCalls::cudaDeviceSynchronize]);
#else
  incr(ApiTest::totalCounter, "cudaDeviceSynchronize");
#endif
  return o_cudaDeviceSynchronize();
}

//Copies data to the given symbol on the device.
__host__ cudaError_t cudaMemcpy(void* dst, const void* src, size_t count, cudaMemcpyKind kind) {
  cudaError_t (*o_cudaMemcpy)(void*, const void*, size_t, cudaMemcpyKind);  
  o_cudaMemcpy = (cudaError_t (*)(void*, const void*, size_t, cudaMemcpyKind))dlsym(RTLD_NEXT, "cudaMemcpy");
  incr(ApiTest::totalCounter, "cudaMemcpy");
  return o_cudaMemcpy(dst, src, count, kind);
}

//Copies data between host and device.
__host__ cudaError_t cudaMemcpy2D(void* dst, size_t dpitch, const void* src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind) {
  cudaError_t (*o_cudaMemcpy2D)(void*, size_t, const void*, size_t, size_t, size_t, cudaMemcpyKind);
  o_cudaMemcpy2D = (cudaError_t (*)(void*, size_t, const void*, size_t, size_t, size_t, cudaMemcpyKind))dlsym(RTLD_NEXT, "cudaMemcpy2D");
  incr(ApiTest::totalCounter, "cudaMemcpy2D");
  return o_cudaMemcpy2D(dst, dpitch, src, spitch, width, height, kind);
}

//Copies data between host and device.
__host__ cudaError_t cudaMemcpy2DArrayToArray ( cudaArray_t dst, size_t wOffsetDst, size_t hOffsetDst, cudaArray_const_t src, size_t wOffsetSrc, size_t hOffsetSrc, size_t width, size_t height, cudaMemcpyKind kind) {
  cudaError_t (*o_cudaMemcpy2DArrayToArray) (cudaArray_t, size_t, size_t, cudaArray_const_t, size_t, size_t, size_t, size_t, cudaMemcpyKind);
  o_cudaMemcpy2DArrayToArray = (cudaError_t (*)(cudaArray_t, size_t, size_t, cudaArray_const_t, size_t, size_t, size_t, size_t, cudaMemcpyKind))dlsym(RTLD_NEXT, "cudaMemcpy2DArrayToArray");
  incr(ApiTest::totalCounter, "cudaMemcpy2DArrayToArray");
  return o_cudaMemcpy2DArrayToArray(dst, wOffsetDst, hOffsetDst, src, wOffsetSrc, hOffsetSrc, width, height, kind);
}

//Copies data between host and device.
__host__ __device__ cudaError_t cudaMemcpy2DAsync ( void* dst, size_t dpitch, const void* src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind, cudaStream_t stream) {
  cudaError_t (*o_cudaMemcpy2DAsync) (void*, size_t, const void*, size_t, size_t, size_t, cudaMemcpyKind, cudaStream_t);
  o_cudaMemcpy2DAsync = (cudaError_t (*)(void*, size_t, const void*, size_t, size_t, size_t, cudaMemcpyKind, cudaStream_t))dlsym(RTLD_NEXT, "cudaMemcpy2DAsync");
#ifdef __CUDA_ARCH__
  Kokkos::atomic_increment(&ApiTest::cudaCounters[ApiTest::devCalls::cudaMemcpy2DAsync]);
#else
  incr(ApiTest::totalCounter, "cudaMemcpy2DAsync");
#endif
  return o_cudaMemcpy2DAsync(dst, dpitch, src, spitch, width, height, kind, stream);
}

//Copies data between host and device.
__host__ cudaError_t cudaMemcpy2DFromArray ( void* dst, size_t dpitch, cudaArray_const_t src, size_t wOffset, size_t hOffset, size_t width, size_t height, cudaMemcpyKind kind ) {
  cudaError_t (*o_cudaMemcpy2DFromArray) ( void*, size_t, cudaArray_const_t, size_t, size_t, size_t, size_t, cudaMemcpyKind);
  o_cudaMemcpy2DFromArray = (cudaError_t (*)(void*, size_t, cudaArray_const_t, size_t, size_t, size_t, size_t, cudaMemcpyKind))dlsym(RTLD_NEXT, "cudaMemcpy2DFromArray");
  incr(ApiTest::totalCounter, "cudaMemcpy2DFromArray");
  return o_cudaMemcpy2DFromArray(dst, dpitch, src, wOffset, hOffset, width, height, kind);
}

//Copies data between host and device.
__host__ cudaError_t cudaMemcpy2DFromArrayAsync ( void* dst, size_t dpitch, cudaArray_const_t src, size_t wOffset, size_t hOffset, size_t width, size_t height, cudaMemcpyKind kind, cudaStream_t stream) {
  cudaError_t (*o_cudaMemcpy2DFromArrayAsync) ( void*, size_t, cudaArray_const_t, size_t, size_t, size_t, size_t, cudaMemcpyKind, cudaStream_t);
  o_cudaMemcpy2DFromArrayAsync = (cudaError_t (*)(void*, size_t, cudaArray_const_t, size_t, size_t, size_t, size_t, cudaMemcpyKind, cudaStream_t))dlsym(RTLD_NEXT, "cudaMemcpy2DFromArrayAsync");
  incr(ApiTest::totalCounter, "cudaMemcpy2DFromArrayAsync");
  return o_cudaMemcpy2DFromArrayAsync(dst, dpitch, src, wOffset, hOffset, width, height, kind, stream);
}

//Copies data between host and device.
__host__ cudaError_t cudaMemcpy2DToArray ( cudaArray_t dst, size_t wOffset, size_t hOffset, const void* src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind ) {
  cudaError_t (*o_cudaMemcpy2DToArray) ( cudaArray_t, size_t, size_t, const void*, size_t, size_t, size_t, cudaMemcpyKind );
  o_cudaMemcpy2DToArray = (cudaError_t (*)(cudaArray_t, size_t, size_t, const void*, size_t, size_t, size_t, cudaMemcpyKind))dlsym(RTLD_NEXT, "cudaMemcpy2DToArray");
  incr(ApiTest::totalCounter, "cudaMemcpy2DToArray");
  return o_cudaMemcpy2DToArray(dst, wOffset, hOffset, src, spitch, width, height, kind);
}

//Copies data between host and device.
__host__ cudaError_t cudaMemcpy2DToArrayAsync ( cudaArray_t dst, size_t wOffset, size_t hOffset, const void* src, size_t spitch, size_t width, size_t height, cudaMemcpyKind kind, cudaStream_t stream) {
  cudaError_t (*o_cudaMemcpy2DToArrayAsync) ( cudaArray_t, size_t, size_t, const void*, size_t, size_t, size_t, cudaMemcpyKind, cudaStream_t );
  o_cudaMemcpy2DToArrayAsync = (cudaError_t (*)(cudaArray_t, size_t, size_t, const void*, size_t, size_t, size_t, cudaMemcpyKind, cudaStream_t))dlsym(RTLD_NEXT, "cudaMemcpy2DToArrayAsync");
  incr(ApiTest::totalCounter, "cudaMemcpy2DToArrayAsync");
  return o_cudaMemcpy2DToArrayAsync(dst, wOffset, hOffset, src, spitch, width, height, kind, stream);
}

//Copies data between 3D objects.
__host__ cudaError_t cudaMemcpy3D ( const cudaMemcpy3DParms* p ) {
  cudaError_t (*o_cudaMemcpy3D) ( const cudaMemcpy3DParms* );
  o_cudaMemcpy3D = (cudaError_t (*)(const cudaMemcpy3DParms*))dlsym(RTLD_NEXT, "cudaMemcpy3D");
  incr(ApiTest::totalCounter, "cudaMemcpy3D");
  return o_cudaMemcpy3D(p);
}

//Copies data between 3D objects.
__host__ __device__ cudaError_t cudaMemcpy3DAsync ( const cudaMemcpy3DParms* p, cudaStream_t stream ) {
  cudaError_t (*o_cudaMemcpy3DAsync) ( const cudaMemcpy3DParms* , cudaStream_t );
  o_cudaMemcpy3DAsync = (cudaError_t (*)(const cudaMemcpy3DParms* , cudaStream_t))dlsym(RTLD_NEXT, "cudaMemcpy3DAsync");
#ifdef __CUDA_ARCH__
  Kokkos::atomic_increment(&ApiTest::cudaCounters[ApiTest::devCalls::cudaMemcpy3DAsync]);
#else
  incr(ApiTest::totalCounter, "cudaMemcpy3DAsync");
#endif
  return o_cudaMemcpy3DAsync(p, stream);
}

//Copies memory between devices.
__host__ cudaError_t cudaMemcpy3DPeer ( const cudaMemcpy3DPeerParms* p ) {
  cudaError_t (*o_cudaMemcpy3DPeer) ( const cudaMemcpy3DPeerParms* );
  o_cudaMemcpy3DPeer = (cudaError_t (*)(const cudaMemcpy3DPeerParms*))dlsym(RTLD_NEXT, "cudaMemcpy3DPeer");
  incr(ApiTest::totalCounter, "cudaMemcpy3DPeer");
  return o_cudaMemcpy3DPeer(p);
}

//Copies memory between devices asynchronously.
__host__ cudaError_t cudaMemcpy3DPeerAsync ( const cudaMemcpy3DPeerParms* p, cudaStream_t stream) {
  cudaError_t (*o_cudaMemcpy3DPeerAsync) ( const cudaMemcpy3DPeerParms*, cudaStream_t );
  o_cudaMemcpy3DPeerAsync = (cudaError_t (*)(const cudaMemcpy3DPeerParms*, cudaStream_t))dlsym(RTLD_NEXT, "cudaMemcpy3DPeerAsync");
  incr(ApiTest::totalCounter, "cudaMemcpy3DPeerAsync");
  return o_cudaMemcpy3DPeerAsync(p, stream);
}

//Copies data between host and device.
__host__ __device__ cudaError_t cudaMemcpyAsync ( void* dst, const void* src, size_t count, cudaMemcpyKind kind, cudaStream_t stream) {
  cudaError_t (*o_cudaMemcpyAsync) ( void*, const void*, size_t, cudaMemcpyKind, cudaStream_t );
  o_cudaMemcpyAsync = (cudaError_t (*)(void*, const void*, size_t, cudaMemcpyKind, cudaStream_t))dlsym(RTLD_NEXT, "cudaMemcpyAsync");
#ifdef __CUDA_ARCH__
  Kokkos::atomic_increment(&ApiTest::cudaCounters[ApiTest::devCalls::cudaMemcpyAsync]);
#else
  incr(ApiTest::totalCounter, "cudaMemcpyAsync");
#endif
  return o_cudaMemcpyAsync(dst, src, count, kind, stream);
}

//Copies data from the given symbol on the device.
__host__ cudaError_t cudaMemcpyFromSymbol ( void* dst, const void* symbol, size_t count, size_t offset, cudaMemcpyKind kind) {
  cudaError_t (*o_cudaMemcpyFromSymbol) ( void*, const void*, size_t, size_t, cudaMemcpyKind );
  o_cudaMemcpyFromSymbol = (cudaError_t (*)( void*, const void*, size_t, size_t, cudaMemcpyKind ))dlsym(RTLD_NEXT, "cudaMemcpyFromSymbol");
  incr(ApiTest::totalCounter, "cudaMemcpyFromSymbol");
  return o_cudaMemcpyFromSymbol(dst, symbol, count, offset, kind);
}

//Copies data from the given symbol on the device.
__host__ cudaError_t cudaMemcpyFromSymbolAsync ( void* dst, const void* symbol, size_t count, size_t offset, cudaMemcpyKind kind, cudaStream_t stream) {
  cudaError_t (*o_cudaMemcpyFromSymbolAsync) ( void*, const void*, size_t, size_t, cudaMemcpyKind, cudaStream_t );
  o_cudaMemcpyFromSymbolAsync = (cudaError_t (*)( void*, const void*, size_t, size_t, cudaMemcpyKind, cudaStream_t ))dlsym(RTLD_NEXT, "cudaMemcpyFromSymbolAsync");
  incr(ApiTest::totalCounter, "cudaMemcpyFromSymbolAsync");
  return o_cudaMemcpyFromSymbolAsync(dst, symbol, count, offset, kind, stream);
}

//Copies memory between two devices.
__host__ cudaError_t cudaMemcpyPeer ( void* dst, int  dstDevice, const void* src, int  srcDevice, size_t count ) {
  cudaError_t (*o_cudaMemcpyPeer) ( void*, int, const void*, int, size_t );
  o_cudaMemcpyPeer = (cudaError_t (*)( void*, int, const void*, int, size_t ))dlsym(RTLD_NEXT, "cudaMemcpyPeer");
  incr(ApiTest::totalCounter, "cudaMemcpyPeer");
  return o_cudaMemcpyPeer(dst, dstDevice, src, srcDevice, count);
}

//Copies memory between two devices asynchronously.
__host__ cudaError_t cudaMemcpyPeerAsync ( void* dst, int  dstDevice, const void* src, int  srcDevice, size_t count, cudaStream_t stream) {
  cudaError_t (*o_cudaMemcpyPeerAsync) ( void*, int, const void*, int, size_t, cudaStream_t );
  o_cudaMemcpyPeerAsync = (cudaError_t (*)( void*, int, const void*, int, size_t, cudaStream_t ))dlsym(RTLD_NEXT, "cudaMemcpyPeerAsync");
  incr(ApiTest::totalCounter, "cudaMemcpyPeerAsync");
  return o_cudaMemcpyPeerAsync(dst, dstDevice, src, srcDevice, count, stream);
}

//Copies data to the given symbol on the device.
__host__ cudaError_t cudaMemcpyToSymbol ( const void* symbol, const void* src, size_t count, size_t offset, cudaMemcpyKind kind ) {
  cudaError_t (*o_cudaMemcpyToSymbol) ( const void*, const void*, size_t, size_t, cudaMemcpyKind );
  o_cudaMemcpyToSymbol = (cudaError_t (*)( const void*, const void*, size_t, size_t, cudaMemcpyKind ))dlsym(RTLD_NEXT, "cudaMemcpyToSymbol");
  incr(ApiTest::totalCounter, "cudaMemcpyToSymbol");
  return o_cudaMemcpyToSymbol(symbol, src, count, offset, kind);
}

//Copies data to the given symbol on the device.
__host__ cudaError_t cudaMemcpyToSymbolAsync ( const void* symbol, const void* src, size_t count, size_t offset, cudaMemcpyKind kind, cudaStream_t stream ) {
  cudaError_t (*o_cudaMemcpyToSymbolAsync) ( const void*, const void*, size_t, size_t, cudaMemcpyKind, cudaStream_t );
  o_cudaMemcpyToSymbolAsync = (cudaError_t (*)( const void*, const void*, size_t, size_t, cudaMemcpyKind, cudaStream_t ))dlsym(RTLD_NEXT, "cudaMemcpyToSymbolAsync");
  incr(ApiTest::totalCounter, "cudaMemcpyToSymbolAsync");
  return o_cudaMemcpyToSymbolAsync(symbol, src, count, offset, kind, stream);
}
