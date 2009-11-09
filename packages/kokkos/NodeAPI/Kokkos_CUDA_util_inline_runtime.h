#ifndef KOKKOS_CUDA_UTIL_INLINE_FUNCTIONS_RUNTIME_H_
#define KOKKOS_CUDA_UTIL_INLINE_FUNCTIONS_RUNTIME_H_

#include <cstdio>
#include <cstdlib>
#include <cuda_runtime.h>

// We define these calls here, so the user doesn't need to include __FILE__ and __LINE__
// The advantage is the developers gets to use the inline function so they can debug
#define cutilSafeCallNoSync(err)     __cudaSafeCallNoSync(err, __FILE__, __LINE__)
#define cutilSafeCall(err)           __cudaSafeCall      (err, __FILE__, __LINE__)
#define cutilSafeThreadSync()        __cudaSafeThreadSync(__FILE__, __LINE__)
#define cutilCheckMsg(msg)           __cutilCheckMsg     (msg, __FILE__, __LINE__)

inline void __cudaSafeCallNoSync( cudaError err, const char *file, const int line )
{
    if( cudaSuccess != err) {
        fprintf(stderr, "cudaSafeCallNoSync() Runtime API error in file <%s>, line %i : %s.\n",
                file, line, cudaGetErrorString( err) );
        exit(-1);
    }
}

inline void __cudaSafeCall( cudaError err, const char *file, const int line )
{
    if( cudaSuccess != err) {
        fprintf(stderr, "cudaSafeCall() Runtime API error in file <%s>, line %i : %s.\n",
                file, line, cudaGetErrorString( err) );
        exit(-1);
    }
}

inline void __cudaSafeThreadSync( const char *file, const int line )
{
    cudaError err = cudaThreadSynchronize();
    if ( cudaSuccess != err) {
        fprintf(stderr, "cudaThreadSynchronize() Driver API error in file '%s' in line %i : %s.\n",
                file, line, cudaGetErrorString( err) );
        exit(-1);
    }
}

inline void __cutilCheckMsg( const char *errorMessage, const char *file, const int line )
{
    cudaError_t err = cudaGetLastError();
    if( cudaSuccess != err) {
        fprintf(stderr, "cutilCheckMsg() CUTIL CUDA error: %s in file <%s>, line %i : %s.\n",
                errorMessage, file, line, cudaGetErrorString( err) );
        exit(-1);
    }
#ifdef _DEBUG
    err = cudaThreadSynchronize();
    if( cudaSuccess != err) {
        fprintf(stderr, "cutilCheckMsg cudaThreadSynchronize error: %s in file <%s>, line %i : %s.\n",
                errorMessage, file, line, cudaGetErrorString( err) );
        exit(-1);
    }
#endif
}

#endif // _CUTIL_INLINE_FUNCTIONS_RUNTIME_H_
