// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions are
// met:
//
// 1. Redistributions of source code must retain the above copyright
// notice, this list of conditions and the following disclaimer.
//
// 2. Redistributions in binary form must reproduce the above copyright
// notice, this list of conditions and the following disclaimer in the
// documentation and/or other materials provided with the distribution.
//
// 3. Neither the name of the Corporation nor the names of the
// contributors may be used to endorse or promote products derived from
// this software without specific prior written permission.
//
// THIS SOFTWARE IS PROVIDED BY SANDIA CORPORATION "AS IS" AND ANY
// EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL SANDIA CORPORATION OR THE
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Questions? Contact Eric T. Phipps (etphipp@sandia.gov).
// 
// ***********************************************************************
// @HEADER

#include "Stokhos_Sacado_Kokkos.hpp"

// storage options
enum Storage_Method { STATIC, STATIC_FIXED, LOCAL, DYNAMIC, DYNAMIC_STRIDED, DYNAMIC_THREADED };
static const int num_storage_method = 6;
static const Storage_Method storage_method_values[] = { STATIC, STATIC_FIXED, LOCAL, DYNAMIC, DYNAMIC_STRIDED, DYNAMIC_THREADED };
static const char *storage_method_names[] = { "static", "static-fixed", "local", "dynamic", "dynamic-strided", "dynamic-threaded" };

template <int MaxSize, typename node_type>
struct MPVectorTypes {
  // Storage types
  typedef Stokhos::StaticStorage<int,double,MaxSize,node_type> static_storage;
  typedef Stokhos::StaticFixedStorage<int,double,MaxSize,node_type> static_fixed_storage;
  typedef Stokhos::LocalStorage<int,double,MaxSize,node_type> local_storage;
  typedef Stokhos::DynamicStorage<int,double,node_type> dynamic_storage;
  typedef Stokhos::DynamicStridedStorage<int,double,node_type> dynamic_strided_storage;
  typedef Stokhos::DynamicThreadedStorage<int,double,node_type> dynamic_threaded_storage;

  // Vector types
  typedef Sacado::MP::Vector<static_storage, node_type> static_vector;
  typedef Sacado::MP::Vector<static_fixed_storage, node_type> static_fixed_vector;
  typedef Sacado::MP::Vector<local_storage, node_type> local_vector;
  typedef Sacado::MP::Vector<dynamic_storage, node_type> dynamic_vector;
  typedef Sacado::MP::Vector<dynamic_strided_storage, node_type> dynamic_strided_vector;
  typedef Sacado::MP::Vector<dynamic_threaded_storage, node_type> dynamic_threaded_vector;
};

template <int MaxSize, typename node_type> struct MPVectorExample {
  static bool 
  run(Storage_Method storage_method, int n, int sz, int nblocks, int nthreads, 
      bool reset, bool print);
};

// Maximum size of expansion -- currently 2, 4, or 8 for LocalStorage
const int MaxSize = 4;
