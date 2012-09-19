/*------------------------------------------------------------------------*/
/*                 Copyright 2010 Sandia Corporation.                     */
/*  Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive   */
/*  license for use of this work by or on behalf of the U.S. Government.  */
/*  Export of this program may require a license from the                 */
/*  United States Government.                                             */
/*------------------------------------------------------------------------*/

#ifndef stk_algsup_CudaMemoryMgr_hpp
#define stk_algsup_CudaMemoryMgr_hpp

#include <stdio.h>
#include <stdexcept>
#include <map>

#include <stk_algsup/CudaCall.hpp>

namespace stk {

/** Helper class for managing CUDA device memory.
 *
 * Tracks persistent mappings between host(cpu) buffers and CUDA device
 * buffers, intended to be used for mapping bucket/field data pointers on
 * the host to buffers on the device.
 * This allows algorithm/bucket loops to repeatedly reference pointers to
 * field data without performing the device allocation every time.
 *
 * This class also allows for creating and using device buffers that
 * are not mapped to host buffers.
 */
class CudaMemoryMgr {
  public:
    /** Constructor */
    CudaMemoryMgr()
     : host_to_device_map(),
       device_to_host_map()
    {}

    /** Destructor
     * Upon destruction this class de-allocates all device-buffers that
     * it was tracking.
     */
    virtual ~CudaMemoryMgr();

#ifdef STK_HAVE_CUDA

    /** Return a device-pointer corresponding to the given host-ptr and size.
     * The returned device-pointer points to a buffer which has been allocated 
     * on the CUDA device with length buf_size*sizeof(T), but not initialized.
     *
     * If a device-pointer has already been allocated for the given host-pointer
     * (by a previous call to this method) then that (previously-allocated) device-pointer
     * is returned.
     */
    template<class T>
    T* get_buffer(const T* host_ptr, size_t buf_size);

    /** Return a device-pointer allocated with the given size.
     * The returned device-pointer points to a buffer which has been allocated 
     * on the CUDA device with length buf_size*sizeof(T), but not initialized.
     */
    template<class T>
    T* get_buffer(size_t buf_size);

    /** Destroy (free) the specified device-pointer.
     *
     * De-allocates the cuda-device buffer.
     */
    template<class T>
    void destroy_buffer(T*& device_ptr);

    /** Copy the contents of the given host-ptr to the given device-ptr.
     * If the given device-ptr is not known (was not created by a previous
     * call to get_buffer), an exception is thrown.
     */
    template<class T>
    void copy_to_buffer(const T* host_ptr, size_t buf_size, T* device_ptr);

    /** Copy the contents of the given device-ptr to the given host-ptr.
     * If the given device-ptr is not known (was not created by a previous
     * call to get_buffer), an exception is thrown.
     */
    template<class T>
    void copy_from_buffer(T* host_ptr, size_t buf_size, const T* device_ptr);

    static CudaMemoryMgr& get_singleton();

#endif

 private:
  std::map<const void*,void*> host_to_device_map;
  std::map<const void*,const void*> device_to_host_map;
};//class CudaMemoryMgr

#ifdef STK_HAVE_CUDA

//------------------------------------------------------------------------------
template<class T>
inline
T* CudaMemoryMgr::get_buffer(const T* host_ptr, size_t buf_size)
{
  T* device_ptr = NULL;

  std::map<const void*,void*>::iterator iter = host_to_device_map.find(host_ptr);

  if (iter == host_to_device_map.end()) {
    void* void_device_ptr = NULL;
    CUDA_CALL( cudaMalloc( &void_device_ptr, sizeof(T)*buf_size) );
    device_ptr = reinterpret_cast<T*>(void_device_ptr);

    host_to_device_map.insert( std::make_pair(host_ptr, device_ptr) );
    device_to_host_map.insert( std::make_pair(device_ptr, host_ptr) );
  }
  else {
    device_ptr = reinterpret_cast<T*>(iter->second);
  }

  return device_ptr;
}

//------------------------------------------------------------------------------
template<class T>
inline
T* CudaMemoryMgr::get_buffer(size_t buf_size)
{
  T* device_ptr = NULL;

  CUDA_CALL( cudaMalloc( (void**)&device_ptr, sizeof(T)*buf_size) );

  device_to_host_map.insert( std::make_pair(device_ptr, NULL) );

  return device_ptr;
}

//------------------------------------------------------------------------------
template<class T>
inline
void CudaMemoryMgr::destroy_buffer(T*& device_ptr)
{
  std::map<const void*,const void*>::iterator iter = device_to_host_map.find(device_ptr);
  if (iter != device_to_host_map.end()) {
    const void* host_ptr = iter->second;
    if (host_ptr != NULL) {
      std::map<const void*,void*>::iterator iter2 = host_to_device_map.find(host_ptr);
      if (iter2 != host_to_device_map.end()) {
        host_to_device_map.erase(iter2);
      }
    }
    CUDA_CALL( cudaFree(device_ptr) );
    device_ptr = NULL;
    device_to_host_map.erase(iter);
  }
}

//------------------------------------------------------------------------------
template<class T>
inline
void CudaMemoryMgr::copy_to_buffer(const T* host_ptr, size_t buf_size, T* device_ptr)
{
  std::map<const void*,const void*>::iterator iter = device_to_host_map.find(device_ptr);
  if (iter == device_to_host_map.end()) {
    //failed to find device_ptr in device_to_host_map
    throw std::runtime_error("CudaMemoryMgr::copy_to_buffer ERROR, device_ptr not known.");
  }

  CUDA_CALL( cudaMemcpy( device_ptr, host_ptr, sizeof(T)*buf_size, cudaMemcpyHostToDevice) );
}

//------------------------------------------------------------------------------
template<class T>
inline
void CudaMemoryMgr::copy_from_buffer(T* host_ptr, size_t buf_size, const T* device_ptr)
{
  std::map<const void*,const void*>::iterator iter = device_to_host_map.find(device_ptr);
  if (iter == device_to_host_map.end()) {
    //failed to find device_ptr in device_to_host_map
    throw std::runtime_error("CudaMemoryMgr::copy_from_buffer ERROR, device_ptr not known.");
  }

  CUDA_CALL( cudaMemcpy( host_ptr, device_ptr, sizeof(T)*buf_size, cudaMemcpyDeviceToHost) );
}

#endif

}//namespace stk

#endif

