#ifndef KOKKOS_CUDANODE_IMPL_HPP_
#define KOKKOS_CUDANODE_IMPL_HPP_

#include <Kokkos_CUDANode.hpp>
#include <cuda.h>
#include <cuda_runtime.h>
#include <Kokkos_CUDA_util_inline_runtime.h>
#include <stdlib.h>
#include <stdexcept>
#include <iostream>

// TODO: consider using cudaMallocHost to allocate page-locked host memory
//       this speeds up transfer between device and host, and could be very 
//       useful in the case of Import/Export multivector operations

/////// CUDABuffer template implementations
/////// CUDANode template implementations
template <class T>
CUDABuffer<T> CUDANode::allocBuffer(unsigned int N) {
  T * hptr = (T*)malloc(sizeof(T)*N+1);   // allocate one extra byte for the flags_ char
  T * dptr = NULL;
  // FINISH: if possible, check that there is room; else, boot someone
  cutilSafeCallNoSync( cudaMalloc( (void**)&dptr, sizeof(T)*N ) );
  char *flags = (char *)(hptr+N);
  (*flags) = DIRTY_DEVC;
  return CUDABuffer<T>(N*sizeof(T),flags,dptr,hptr);
}

template <class T>
void CUDANode::freeBuffer(CUDABuffer<T> buf) {
  if (buf.host_ptr_ != NULL) {
    free(buf.host_ptr_); buf.host_ptr_ = NULL;
  }
  if (buf.devc_ptr_ != NULL) {
    cutilSafeCallNoSync( cudaFree(buf.devc_ptr_) ); buf.devc_ptr_ = NULL;
    // FINISH
    // this may return one of the following; check
    // - cudaSuccess 
    // - cudaErrorInvalidDevicePointer 
    // - cudaErrorInitializationError 
  }
}

template <class T>
void CUDANode::copyFromBuffer(unsigned int size, typename buffer<const T>::buffer_t buff, unsigned int offset, T *ptr)
{
  if (buff.length_ < sizeof(T)*(offset+size)) {
    throw std::runtime_error("copyFromBuffer: invalid copy.");
  }
  std::cout << "Copying from buffer" << std::endl;
  cutilSafeCallNoSync( cudaMemcpy( ptr, buff.devc_ptr_+offset, size*sizeof(T), cudaMemcpyDeviceToHost) );
}

template <class T>
void CUDANode::copyToBuffer(unsigned int size, const T *ptr, typename buffer<T>::buffer_t buff, unsigned int offset)
{
  if (buff.length_ < sizeof(T)*(offset+size)) {
    throw std::runtime_error("copyToBuffer: invalid copy.");
  }
  std::cout << "Copying to buffer" << std::endl;
  cutilSafeCallNoSync( cudaMemcpy( buff.devc_ptr_+offset, ptr, size*sizeof(T), cudaMemcpyHostToDevice) );
}

template <class T>
void CUDANode::copyBuffers(unsigned int size, typename buffer<const T>::buffer_t src, unsigned int src_offset, 
                          typename buffer<T>::buffer_t dest, unsigned int dest_offset)
{
  if (sizeof(T)*(dest.dest_offset + size) > dest.length_ || sizeof(T)*(src.src_offset + size) > src.length_) {
    throw std::runtime_error("copyBuffers: invalid copy.");
  }
  std::cout << "Copying between buffers" << std::endl;
  cutilSafeCallNoSync( cudaMemcpy( dest.devc_ptr_+dest_offset, src.devc_ptr_+src_offset, size*sizeof(T), cudaMemcpyDeviceToDevice) );
}

template <class T>
const T * CUDANode::viewBufferConst(unsigned int size, typename buffer<const T>::buffer_t buff, unsigned int offset)
{
  // copy from device to host, if device is dirty.
  if (buff.length_ < sizeof(T)*(offset+size)) {
    throw std::runtime_error("viewBuffer(): invalid view request");
  }
  if ( (*buff.flags_) == DIRTY_DEVC ) {
    copyFromBuffer<char>(buff.length_, reinterpret_cast<CUDABuffer<char> &>(buff), 0, (char*)(buff.host_ptr_));
    (*buff.flags_) = CLEAN;
  }
  // out_views_[(const void *)(buff.host_ptr_+offset)] = buff;
  return buff.host_ptr_+offset;
}

template <class T>
T * CUDANode::viewBuffer(bool writeOnly, unsigned int size, typename buffer<T>::buffer_t buff, unsigned int offset)
{
  // copy from device to host, if device is dirty and there is the potential for reading
  if (buff.length_ < sizeof(T)*(offset+size)) {
    throw std::runtime_error("viewBuffer(): invalid view request");
  }
  if ( writeOnly == false && (*buff.flags_) == DIRTY_DEVC ) {
    copyFromBuffer<char>(buff.length_, reinterpret_cast<CUDABuffer<char> &>(buff), 0, (char *)buff.host_ptr_);
    (*buff.flags_) = CLEAN;
  }
  (*buff.flags_) = DIRTY_HOST;
  // out_views_[(const void *)(buff.host_ptr_+offset)] = buff;
  return buff.host_ptr_+offset;
}

template <class T>
void CUDANode::releaseView(const T * view_ptr)
{
  // FINISH: must delete the view from the map as well
  // std::map<const void *, CUDABuffer<const void> >::iterator it;
  // it = out_views_.find( (const void *) view_ptr );
  // if (it == out_views_.end()) {
  //   throw std::runtime_error("freed invalid buffer view");
  // }
  // FINISH: this is killing performance
  // CUDABuffer<T> &buff = reinterpret_cast<CUDABuffer<T> &>( it->second );
  // if ( (*buff.flags_) == DIRTY_HOST) {
  //   copyToBuffer<char>(buff.length_, (const char *)buff.host_ptr_, reinterpret_cast<CUDABuffer<char> &>(buff), 0);
  //   (*buff.flags_) = CLEAN;
  // }
}

#endif
