#ifndef cuda_memory_h
#define cuda_memory_h

#include <vector>

float* cuda_alloc_float(size_t n, const float& init_value = 0.0);

void copy_cuda_memory_to_host(size_t n, const float* device_mem, std::vector<float>& host_vec);

#endif
