/*
//@HEADER
// ************************************************************************
//
//               KokkosKernels 0.9: Linear Algebra and Graph Kernels
//                 Copyright 2017 Sandia Corporation
//
// Under the terms of Contract DE-AC04-94AL85000 with Sandia Corporation,
// the U.S. Government retains certain rights in this software.
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
// Questions? Contact Siva Rajamanickam (srajama@sandia.gov)
//
// ************************************************************************
//@HEADER
*/
#include <Kokkos_Atomic.hpp>
#include <atomic>
namespace KokkosKernels{

namespace Experimental{


template <typename size_type, typename key_type, typename value_type>
struct HashmapAccumulator{
  size_type hash_key_size;
  size_type max_value_size;
  size_type used_size;

  size_type *hash_begins;
  size_type *hash_nexts;
  key_type *keys;
  value_type *values;
  const int INSERT_SUCCESS;
  const int INSERT_FULL;

  /**
   * Assumption: hash_begins_ are all initialized to -1.
   */
  KOKKOS_INLINE_FUNCTION
  HashmapAccumulator ():
        hash_key_size(),
        max_value_size(),
        used_size(0),
        hash_begins(),
        hash_nexts(),
        keys(),
        values(), INSERT_SUCCESS(0), INSERT_FULL(1){}

  /**
   * Assumption: hash_begins_ are all initialized to -1.
   */
  KOKKOS_INLINE_FUNCTION
  HashmapAccumulator (
      const size_type hash_key_size_,
      const size_type value_size_,
      size_type *hash_begins_,
      size_type *hash_nexts_,
      key_type *keys_,
      value_type *values_):

        hash_key_size(hash_key_size_),
        max_value_size(value_size_),
        used_size(0),
        hash_begins(hash_begins_),
        hash_nexts(hash_nexts_),
        keys(keys_),
        values(values_), INSERT_SUCCESS(0), INSERT_FULL(1){}


  //function to be called from device.
  //Accumulation is OR operation.
  //Insertion is sequential, no race condition for the insertion.
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_mergeOr_TrackHashes (
      size_type hash,
      key_type key,
      value_type value,

      size_type *used_size_,
      const size_type max_value_size_,
      size_type *used_hash_size,
      size_type *used_hashes){

    size_type i = hash_begins[hash];
    for (; i != -1; i = hash_nexts[i]){
      if (keys[i] == key){
        values[i] = values[i] | value;
        return INSERT_SUCCESS;
      }
    }

    if (*used_size_ >= max_value_size_) return INSERT_FULL;
    size_type my_index = (*used_size_)++;

    if (hash_begins[hash] == -1){
      used_hashes[used_hash_size[0]++] = hash;
    }
    hash_nexts[my_index] = hash_begins[hash];

    hash_begins[hash] = my_index;
    keys[my_index] = key;
    values[my_index] = value;
    return INSERT_SUCCESS;
  }

  //function to be called from device.
  //Accumulation is OR operation.
  //TODO: This function is for triangle counting.
  //Assume that there are 2 values for triangle count.
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_mergeOr_TriangleCount_TrackHashes (
      size_type hash,
      key_type key,
      value_type value,
      value_type *values2,
      size_type *used_size_,
      const size_type max_value_size_,
      size_type *used_hash_size,
      size_type *used_hashes){

    size_type i = hash_begins[hash];
    for (; i != -1; i = hash_nexts[i]){
      if (keys[i] == key){
        values2[i] = values2[i] | (values[i] & value);
        values[i] = values[i] | value;
        return INSERT_SUCCESS;
      }
    }

    if (*used_size_ >= max_value_size_) return INSERT_FULL;
    size_type my_index = (*used_size_)++;

    if (hash_begins[hash] == -1){
      used_hashes[used_hash_size[0]++] = hash;
    }
    hash_nexts[my_index] = hash_begins[hash];

    hash_begins[hash] = my_index;
    keys[my_index] = key;
    values[my_index] = value;
    values2[my_index] = 0;
    return INSERT_SUCCESS;
  }


  //this is used in slow triangle counting method.
  //L x Incidence
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_mergeAnd_TriangleCount_TrackHashes (
      size_type hash,
      key_type key,
      value_type value,
      value_type *values2,
      size_type *used_size_,
      const size_type max_value_size_,
      size_type *used_hash_size,
      size_type *used_hashes){
    //this function will only try to do an AND operation with
    //existing keys. If the key is not there, returns INSERT_FULL.
    size_type i = hash_begins[hash];
    for (; i != -1; i = hash_nexts[i]){
      if (keys[i] == key){
        //values2[i] = values2[i] | (values[i] & value);
        values[i] = values[i] & value;
        ++values2[i];
        return INSERT_SUCCESS;
      }
    }
    return INSERT_FULL;
  }

  //this is used in LxL or Incidence^T x L
  KOKKOS_INLINE_FUNCTION
  value_type sequential_insert_into_hash_mergeAnd_TriangleCount_TrackHashes (
      size_type hash,
      key_type key,
      value_type value){
    //this function will only try to do an AND operation with
    //existing keys. If the key is not there, returns INSERT_FULL.
    size_type i = hash_begins[hash];
    for (; i != -1; i = hash_nexts[i]){

      if (keys[i] == key){
        return values[i] & value;
      }
    }
    return 0;
  }

  //this is used in slow triangle counting method.
  //L x Incidence
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_TriangleCount_TrackHashes (
      size_type hash,
      key_type key,
      value_type value,
      value_type *values2,
      size_type *used_size_,
      const size_type max_value_size_,
      size_type *used_hash_size,
      size_type *used_hashes){

    //this function will directly insert, won't check if it exists already.
    if (*used_size_ >= max_value_size_) return INSERT_FULL;
    size_type my_index = (*used_size_)++;

    keys[my_index] = key;
    values[my_index] = value;
    values2[my_index] = 1;

    if (hash_begins[hash] == -1){
      hash_begins[hash] = my_index;
      used_hashes[used_hash_size[0]++] = hash;
    }
    else {
      hash_nexts[my_index] = hash_begins[hash];
      hash_begins[hash] = my_index;

    }

    return INSERT_SUCCESS;
  }

  //this is used in LxL or Incidence^T x L
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_TriangleCount_TrackHashes (
      size_type hash,
      key_type key,
      value_type value,
      size_type *used_size_,
      const size_type max_value_size_,
      size_type *used_hash_size,
      size_type *used_hashes){

    //this function will directly insert, won't check if it exists already.
    if (*used_size_ >= max_value_size_) return INSERT_FULL;
    size_type my_index = (*used_size_)++;

    keys[my_index] = key;
    values[my_index] = value;

    if (hash_begins[hash] == -1){
      hash_begins[hash] = my_index;
      used_hashes[used_hash_size[0]++] = hash;
    }
    else {
      hash_nexts[my_index] = hash_begins[hash];
      hash_begins[hash] = my_index;
    }

    return INSERT_SUCCESS;
  }

  //function to be called from device.
  //Insertion is sequential, no race condition for the insertion.
  //the mergeadd used in the numeric of KKMEM.
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_mergeAdd_TrackHashes (
      size_type hash,
      key_type key,
      value_type value,

      size_type *used_size_,
      const size_type max_value_size_,
      size_type *used_hash_size,
      size_type *used_hashes){

    size_type i = hash_begins[hash];
    for (; i != -1; i = hash_nexts[i]){
      if (keys[i] == key){
        values[i] = values[i] + value;
        return INSERT_SUCCESS;
      }
    }

    //if (*used_size_ >= max_value_size_) return INSERT_FULL;
    size_type my_index = (*used_size_)++;

    if (hash_begins[hash] == -1){
      used_hashes[used_hash_size[0]++] = hash;
    }
    hash_nexts[my_index] = hash_begins[hash];

    hash_begins[hash] = my_index;
    keys[my_index] = key;
    values[my_index] = value;
    return INSERT_SUCCESS;
  }


  //no values. simply adds to the keys. 
  //used in the compression to count the sets. 
  //also used in the symbolic of spgemm if no compression is applied.
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_TrackHashes (
      size_type hash,
      key_type key,

      size_type *used_size_,
      const size_type max_value_size_,
      size_type *used_hash_size,
      size_type *used_hashes){

    size_type i = hash_begins[hash];
    for (; i != -1; i = hash_nexts[i]){
      if (keys[i] == key){
        return INSERT_SUCCESS;
      }
    }

    //if (*used_size_ >= max_value_size_) return INSERT_FULL;
    size_type my_index = (*used_size_)++;

    if (hash_begins[hash] == -1){
      used_hashes[used_hash_size[0]++] = hash;
    }
    hash_nexts[my_index] = hash_begins[hash];

    hash_begins[hash] = my_index;
    keys[my_index] = key;
    return INSERT_SUCCESS;
  }



  //used in the kkmem's numeric phase for second level hashmaps.
  //function to be called from device.
  //Accumulation is Add operation. It is not atomicAdd, as this
  //is for the cases where we know that none of the simultanous
  //insertions will have the same key.
  //Insertion is simulteanous for the vector lanes of a thread.
  //used_size should be a shared pointer among the thread vectors
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_mergeAdd_TrackHashes (
      const team_member_t & teamMember,
      const int vector_size,

      size_type &hash,
      const key_type key,
      const value_type value,
      volatile size_type *used_size_,
      const size_type max_value_size_,
      size_type *used_hash_size,
      size_type *used_hashes
      ){



    if (hash != -1){
      size_type i = hash_begins[hash];

      for (; i != -1; i = hash_nexts[i]){
        if (keys[i] == key){
          values[i] = values[i] + value;
          return INSERT_SUCCESS;
        }
      }
    }
    else {
        return INSERT_SUCCESS;
    }

    size_type my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));


    if (my_write_index >= max_value_size_) {
      return INSERT_FULL;
    }
    else {

      keys[my_write_index] = key;
      values[my_write_index] = value;
#if defined(KOKKOS_ARCH_VOLTA) || defined(KOKKOS_ARCH_VOLTA70) || defined(KOKKOS_ARCH_VOLTA72) 
      //this is an issue on VOLTA because warps do not go in SIMD fashion anymore. 
      //while some thread might insert my_write_index into linked list, another 
      //thread in the warp might be reading keys in above loop.
      //before inserting the new value in liked list -- which is done with atomic exchange below,
      //we make sure that the linked is is complete my assigning the hash_next to current head.
      //the head might be different when we do the atomic exchange.
      //this would cause temporarily skipping a key in the linkedlist until 
      //hash_nexts is updated second time as below.
      //but this is okay for spgemm,
      //because no two keys will be inserted into hashmap at the same time, as rows have unique columns. 
      hash_nexts[my_write_index] = hash_begins[hash];
#endif
      size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
      if (hashbeginning == -1){
        used_hashes[Kokkos::atomic_fetch_add(used_hash_size, size_type(1))] = hash;
      }
      hash_nexts[my_write_index] = hashbeginning;
      return INSERT_SUCCESS;
    }
  }


  //used in kkmem's numeric phase to insert to first level hashmaps.
  //function to be called from device.
  //Accumulation is Add operation. It is not atomicAdd, as this
  //is for the cases where we know that none of the simultanous
  //insertions will have the same key.
  //Insertion is simulteanous for the vector lanes of a thread.
  //used_size should be a shared pointer among the thread vectors
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_mergeAdd (
      const team_member_t & teamMember,
      const int vector_size,
      size_type &hash,
      const key_type key,
      const value_type value,
      volatile size_type *used_size_,
      const size_type max_value_size_ ){
	  if (hash != -1){
		  size_type i = hash_begins[hash];

		  for (; i != -1; i = hash_nexts[i]){
			  if (keys[i] == key){
				  values[i] = values[i] + value;
				  return INSERT_SUCCESS;
			  }
		  }
	  }
	  else {
		  return INSERT_SUCCESS;
	  }


	  if (used_size_[0] >= max_value_size_){
		  return INSERT_FULL;
	  }
	  size_type my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));


	  if (my_write_index >= max_value_size_) {
		  return INSERT_FULL;
	  }
	  else {

		  keys[my_write_index] = key;
		  values[my_write_index] = value;
#if defined(KOKKOS_ARCH_VOLTA) || defined(KOKKOS_ARCH_VOLTA70) || defined(KOKKOS_ARCH_VOLTA72)  
		  //this is an issue on VOLTA because warps do not go in SIMD fashion anymore. 
		  //while some thread might insert my_write_index into linked list, another 
		  //thread in the warp might be reading keys in above loop.
		  //before inserting the new value in liked list -- which is done with atomic exchange below,
		  //we make sure that the linked is is complete my assigning the hash_next to current head.
		  //the head might be different when we do the atomic exchange.
		  //this would cause temporarily skipping a key in the linkedlist until 
		  //hash_nexts is updated second time as below.
		  //but this is okay for spgemm,
		  //because no two keys will be inserted into hashmap at the same time, as rows have unique columns. 
		  hash_nexts[my_write_index] = hash_begins[hash];
#endif

		  size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
		  hash_nexts[my_write_index] = hashbeginning;
		  return INSERT_SUCCESS;
	  }
  }

  //used in symbolic of kkmem if the compression is not applied.
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash (
      const team_member_t & teamMember,
      const int &vector_size,

      const size_type &hash,
      const key_type &key,
      volatile size_type *used_size_,
      const size_type &max_value_size_
      ){


    if (hash != -1){
      size_type i = hash_begins[hash];
      for (; i != -1; i = hash_nexts[i]){
        if (keys[i] == key){
          return INSERT_SUCCESS;
        }
      }
    }
    else {
        return INSERT_SUCCESS;
    }

    size_type my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));

    if (my_write_index >= max_value_size_) {
      return INSERT_FULL;
    }
    else {

      keys[my_write_index] = key;
#if defined(KOKKOS_ARCH_VOLTA) || defined(KOKKOS_ARCH_VOLTA70) || defined(KOKKOS_ARCH_VOLTA72)  
      //this is an issue on VOLTA because warps do not go in SIMD fashion anymore. 
      //while some thread might insert my_write_index into linked list, another 
      //thread in the warp might be reading keys in above loop.
      //before inserting the new value in liked list -- which is done with atomic exchange below,
      //we make sure that the linked is is complete my assigning the hash_next to current head.
      //the head might be different when we do the atomic exchange.
      //this would cause temporarily skipping a key in the linkedlist until 
      //hash_nexts is updated second time as below.
      //but this is okay for spgemm,
      //because no two keys will be inserted into hashmap at the same time, as rows have unique columns. 
      hash_nexts[my_write_index] = hash_begins[hash];
#endif

      size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
      hash_nexts[my_write_index] = hashbeginning;
      return INSERT_SUCCESS;
    }
  }

  //function to be called from device.
  //Accumulation is Add operation. It is not atomicAdd, as this
  //is for the cases where we know that none of the simultanous
  //insertions will have the same key.
  //Insertion is simulteanous for the vector lanes of a thread.
  //used_size should be a shared pointer among the thread vectors
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_mergeOr (
      const team_member_t & teamMember,
      const int &vector_size,

      const size_type &hash,
      const key_type &key,
      const value_type &value,
      volatile size_type *used_size_,
      const size_type &max_value_size_
      ){


    if (hash != -1){
      size_type i = hash_begins[hash];
      for (; i != -1; i = hash_nexts[i]){
        if (keys[i] == key){
          values[i] = (key_type)values[i] | (key_type)value;
          return INSERT_SUCCESS;
        }
      }
    }
    else {
        return INSERT_SUCCESS;
    }

    size_type my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));

    if (my_write_index >= max_value_size_) {
      return INSERT_FULL;
    }
    else {

      keys[my_write_index] = key;
      values[my_write_index] = value;
#if defined(KOKKOS_ARCH_VOLTA) || defined(KOKKOS_ARCH_VOLTA70) || defined(KOKKOS_ARCH_VOLTA72)  
      //this is an issue on VOLTA because warps do not go in SIMD fashion anymore. 
      //while some thread might insert my_write_index into linked list, another 
      //thread in the warp might be reading keys in above loop.
      //before inserting the new value in liked list -- which is done with atomic exchange below,
      //we make sure that the linked is is complete my assigning the hash_next to current head.
      //the head might be different when we do the atomic exchange.
      //this would cause temporarily skipping a key in the linkedlist until 
      //hash_nexts is updated second time as below.
      //but this is okay for spgemm,
      //because no two keys will be inserted into hashmap at the same time, as rows have unique columns. 
      hash_nexts[my_write_index] = hash_begins[hash];
#endif

      size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
      hash_nexts[my_write_index] = hashbeginning;
      return INSERT_SUCCESS;
    }
  }


  //function to be called from device.
  //Accumulation is Add operation. It is not atomicAdd, as this
  //is for the cases where we know that none of the simultanous
  //insertions will have the same key.
  //Insertion is simulteanous for the vector lanes of a thread.
  //used_size should be a shared pointer among the thread vectors
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_mergeOr_TrackHashes (
      const team_member_t & teamMember,
      const int &vector_size,

      const size_type &hash,
      const key_type &key,
      const value_type &value,
      volatile size_type *used_size_,
      const size_type &max_value_size_,
      size_type *used_hash_size,
      size_type *used_hashes){

	  if (hash != -1){
		  size_type i = hash_begins[hash];
		  for (; i != -1; i = hash_nexts[i]){
			  if (keys[i] == key){
				  values[i] = (key_type)values[i] | (key_type)value;
				  return INSERT_SUCCESS;
			  }
		  }
	  }
	  else {
		  return INSERT_SUCCESS;
	  }

	  size_type my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));

	  if (my_write_index >= max_value_size_) {
		  return INSERT_FULL;
	  }
	  else {

		  keys[my_write_index] = key;
		  values[my_write_index] = value;
#if defined(KOKKOS_ARCH_VOLTA) || defined(KOKKOS_ARCH_VOLTA70) || defined(KOKKOS_ARCH_VOLTA72)  
		  //this is an issue on VOLTA because warps do not go in SIMD fashion anymore. 
		  //while some thread might insert my_write_index into linked list, another 
		  //thread in the warp might be reading keys in above loop.
		  //before inserting the new value in liked list -- which is done with atomic exchange below,
		  //we make sure that the linked is is complete my assigning the hash_next to current head.
		  //the head might be different when we do the atomic exchange.
		  //this would cause temporarily skipping a key in the linkedlist until 
		  //hash_nexts is updated second time as below.
		  //but this is okay for spgemm,
		  //because no two keys will be inserted into hashmap at the same time, as rows have unique columns. 
		  hash_nexts[my_write_index] = hash_begins[hash];
#endif


		  size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
		  if (hashbeginning == -1){
			  used_hashes[Kokkos::atomic_fetch_add(used_hash_size, size_type(1))] = hash;
		  }
		  hash_nexts[my_write_index] = hashbeginning;
		  return INSERT_SUCCESS;
	  }
  }

  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_TrackHashes (
      const team_member_t & teamMember,
      const int &vector_size,

      const size_type &hash,
      const key_type &key,
      volatile size_type *used_size_,
      const size_type &max_value_size_,
      size_type *used_hash_size,
      size_type *used_hashes){

	  if (hash != -1){
		  size_type i = hash_begins[hash];
		  for (; i != -1; i = hash_nexts[i]){
			  if (keys[i] == key){
				  //values[i] = (key_type)values[i] | (key_type)value;
				  return INSERT_SUCCESS;
			  }
		  }
	  }
	  else {
		  return INSERT_SUCCESS;
	  }

	  size_type my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));

	  if (my_write_index >= max_value_size_) {
		  return INSERT_FULL;
	  }
	  else {

		  keys[my_write_index] = key;
#if defined(KOKKOS_ARCH_VOLTA) || defined(KOKKOS_ARCH_VOLTA70) || defined(KOKKOS_ARCH_VOLTA72)  
                  //this is an issue on VOLTA because warps do not go in SIMD fashion anymore. 
                  //while some thread might insert my_write_index into linked list, another 
                  //thread in the warp might be reading keys in above loop.
                  //before inserting the new value in liked list -- which is done with atomic exchange below,
                  //we make sure that the linked is is complete my assigning the hash_next to current head.
                  //the head might be different when we do the atomic exchange.
                  //this would cause temporarily skipping a key in the linkedlist until 
                  //hash_nexts is updated second time as below.
                  //but this is okay for spgemm,
                  //because no two keys will be inserted into hashmap at the same time, as rows have unique columns. 
                  hash_nexts[my_write_index] = hash_begins[hash];
#endif

		  size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
		  if (hashbeginning == -1){
			  used_hashes[Kokkos::atomic_fetch_add(used_hash_size, size_type(1))] = hash;
		  }
		  hash_nexts[my_write_index] = hashbeginning;
		  return INSERT_SUCCESS;
	  }
  }

};



}
}

