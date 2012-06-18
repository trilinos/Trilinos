// @HEADER
// ***********************************************************************
// 
//                           Stokhos Package
//                 Copyright (2009) Sandia Corporation
// 
// Under terms of Contract DE-AC04-94AL85000, there is a non-exclusive
// license for use of this work by or on behalf of the U.S. Government.
// 
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License as
// published by the Free Software Foundation; either version 2.1 of the
// License, or (at your option) any later version.
//  
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//  
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307
// USA
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
