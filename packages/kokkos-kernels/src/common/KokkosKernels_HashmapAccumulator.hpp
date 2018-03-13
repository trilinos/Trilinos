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

  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_merge (
      size_type hash,
      key_type key,

      size_type *used_size_,
      const size_type max_value_size_){

    size_type i = hash_begins[hash];
    for (; i != -1; i = hash_nexts[i]){
      if (keys[i] == key){
        return INSERT_SUCCESS;
      }
    }

    if (*used_size_ >= max_value_size_) return INSERT_FULL;
    size_type my_index = (*used_size_)++;

    hash_nexts[my_index] = hash_begins[hash];
    hash_begins[hash] = my_index;
    keys[my_index] = key;
    return INSERT_SUCCESS;
  }

  //function to be called from device.
  //Accumulation is OR operation.
  //Insertion is sequential, no race condition for the insertion.
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_mergeOr (
      size_type hash,
      key_type key,
      value_type value,

      size_type *used_size_,
      const size_type max_value_size_){

    size_type i = hash_begins[hash];
    for (; i != -1; i = hash_nexts[i]){
      if (keys[i] == key){
        values[i] = values[i] | value;
        return INSERT_SUCCESS;
      }
    }

    if (*used_size_ >= max_value_size_) return INSERT_FULL;
    size_type my_index = (*used_size_)++;

    hash_nexts[my_index] = hash_begins[hash];
    hash_begins[hash] = my_index;
    keys[my_index] = key;
    values[my_index] = value;
    return INSERT_SUCCESS;
  }

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
  //Accumulation is AND operation.
  //Insertion is sequential, no race condition for the insertion.
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_mergeAnd (
      size_type hash,
      key_type key,
      value_type value){

    size_type i = hash_begins[hash];
    //if it is in the hash perform an and operation.
    for (; i != -1; i = hash_nexts[i]){
      if (keys[i] == key){
        values[i] = values[i] & value;
        return INSERT_SUCCESS;
      }
    }
    //if it is not there, return it is full in the meaning that it is not there.
    return INSERT_FULL;
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
  //Accumulation is OR operation.
  //Insertion is sequential, no race condition for the insertion.
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

  //function to be called from device.
  //Accumulation is OR operation.
  //Insertion is sequential, no race condition for the insertion.
  KOKKOS_INLINE_FUNCTION
  int sequential_sorted_insert_into_hash_mergeAdd_TrackHashes (
      size_type hash,
      key_type key,
      value_type value,

      size_type *used_size_,
      const size_type max_value_size_,
      size_type *used_hash_size,
      size_type *used_hashes){

    size_type i = hash_begins[hash];

    if (i != -1){
        if (keys[i] == key){
          values[i] = values[i] + value;
          return INSERT_SUCCESS;
        }
        else if (keys[i] < key){
        	//if not equal at the same time not larger than the key, than rest is gouing to be even smaller so insert here.
            size_type my_index = (*used_size_)++;
            hash_nexts[my_index] = hash_begins[hash];
            hash_begins[hash] = my_index;
            keys[my_index] = key;
            values[my_index] = value;
        	return INSERT_SUCCESS;
        }
    }
    else {
    	//insert.
        size_type my_index = (*used_size_)++;
        used_hashes[used_hash_size[0]++] = hash;
        hash_nexts[my_index] = hash_begins[hash];
        hash_begins[hash] = my_index;
        keys[my_index] = key;
        values[my_index] = value;
    	return INSERT_SUCCESS;
    }
    size_type prev = i;
    i = hash_nexts[i];
    for (; i != -1; i = hash_nexts[i]){
      if (keys[i] == key){
        values[i] = values[i] + value;
        return INSERT_SUCCESS;
      } else if (keys[i] < key){
    	  size_type my_index = (*used_size_)++;
    	  hash_nexts[my_index] = hash_nexts[prev];
    	  hash_nexts[prev] = my_index;
    	  keys[my_index] = key;
    	  values[my_index] = value;
    	  return INSERT_SUCCESS;
      }
      prev = i;
    }

    size_type my_index = (*used_size_)++;
    hash_nexts[my_index] = hash_nexts[prev];
    hash_nexts[prev] = my_index;
    keys[my_index] = key;
    values[my_index] = value;
    return INSERT_SUCCESS;
  }


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


  //function to be called from device.
  //Accumulation is add operation.
  //Insertion is sequential, no race condition for the insertion.
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_mergeAdd (
      const size_type hash,
      const key_type key,
      const value_type value,

      size_type *used_size_,
      const size_type max_value_size_){


    size_type i = hash_begins[hash];
    for (; i != -1; i = hash_nexts[i]){
      if (keys[i] == key){
        values[i] = values[i] + value;
        return INSERT_SUCCESS;
      }
    }

    if (*used_size_ >= max_value_size_) return INSERT_FULL;
    size_type my_index = (*used_size_)++;

    hash_nexts[my_index] = hash_begins[hash];
    hash_begins[hash] = my_index;
    keys[my_index] = key;
    values[my_index] = value;
    return INSERT_SUCCESS;
  }



  //function to be called from device.
  //Accumulation is Or operation. It is not atomicOr, as this
  //is for the cases where we know that none of the simultanous
  //insertions will have the same key.
  //Insertion is simulteanous for the vector lanes of a thread.
  //used_size should be a shared pointer among the thread vectors
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int vector_insert_into_hash_mergeOr (
      const team_member_t & teamMember,
      const int vector_size,

      size_type &hash,
      const key_type key,
      const value_type value,
      size_type *used_size_,
      const size_type max_value_size_){

    char key_not_found = 1;
    if (hash != -1){
      size_type i = hash_begins[hash];

      for (; i != -1; i = hash_nexts[i]){
        if (keys[i] == key){
          values[i] = values[i] | value;
          key_not_found = 0;
          break;
        }
      }
    }
    else {
      key_not_found = 0;
    }

    int write_pos = 0;
    Kokkos::parallel_scan(
        Kokkos::ThreadVectorRange(teamMember, vector_size),
        [&] (const int threadid, int &update, const bool final) {
      if (final){
        write_pos = update;
      }
      update += key_not_found;
    });

    size_type my_write_index = (*used_size_) + write_pos;
    int num_writes = 0;
    Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &num_writes_) {
          num_writes_ += key_not_found;
        }, num_writes);

    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      (*used_size_) += num_writes;
    });

    if (key_not_found == 0) return INSERT_SUCCESS;


    if (my_write_index >= max_value_size_) {
      return INSERT_FULL;
    }
    else {
      keys[my_write_index] = key;
      values[my_write_index] = value;

      size_type hashbeginning = hash_begins[hash];
      hash_begins[hash] = my_write_index;
      while(hash_begins[hash] != my_write_index){
        hashbeginning = hash_begins[hash];
        hash_begins[hash] = my_write_index;
      }
      hash_nexts[my_write_index] = hashbeginning;
    }
    hash = -1;
    return INSERT_SUCCESS;
  }


  //function to be called from device.
  //Accumulation is Add operation. It is not atomicAdd, as this
  //is for the cases where we know that none of the simultanous
  //insertions will have the same key.
  //Insertion is simulteanous for the vector lanes of a thread.
  //used_size should be a shared pointer among the thread vectors
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int vector_insert_into_hash_mergeAdd (
      const team_member_t & teamMember,
      const int vector_size,

      size_type &hash,
      const key_type key,
      const value_type value,
      size_type *used_size_,
      const size_type max_value_size_){

    char key_not_found = 1;
    if (hash != -1){
      size_type i = hash_begins[hash];

      for (; i != -1; i = hash_nexts[i]){
        if (keys[i] == key){
          values[i] = values[i] + value;
          key_not_found = 0;
          break;
        }
      }
    }
    else {
      key_not_found = 0;
    }

    int write_pos = 0;
    Kokkos::parallel_scan(
        Kokkos::ThreadVectorRange(teamMember, vector_size),
        [&] (const int threadid, int &update, const bool final) {
      if (final){
        write_pos = update;
      }
      update += key_not_found;
    });

    size_type my_write_index = (*used_size_) + write_pos;
    int num_writes = 0;
    Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &num_writes_) {
          if (key_not_found){
            num_writes_ += 1;
          }
        }, num_writes);

    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      (*used_size_) += num_writes;
    });

    if (key_not_found == 0) return INSERT_SUCCESS;


    if (my_write_index >= max_value_size_) {
      return INSERT_FULL;
    }
    else {

      keys[my_write_index] = key;
      values[my_write_index] = value;
      size_type hashbeginning = hash_begins[hash];
      hash_begins[hash] = my_write_index;

      //printf("key:%d write:%d hash:%d hash_begins[hash] :%d\n", key, my_write_index, hash, hash_begins[hash] );
      while(hash_begins[hash] != my_write_index){
        hashbeginning = hash_begins[hash];
        hash_begins[hash] = my_write_index;
      }
      hash_nexts[my_write_index] = hashbeginning;
    }
    hash = -1;
    return INSERT_SUCCESS;
  }

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
      //,bool print = false
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


    /*
    Kokkos::parallel_scan(
        Kokkos::ThreadVectorRange(teamMember, vector_size),
        [&] (const int threadid, int &update, const bool final) {
      if (final){
        write_pos = update;
      }
      update += key_not_found;
    });

    size_type my_write_index = (*used_size_) + write_pos;
    int num_writes = 0;
    Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &num_writes_) {
          if (key_not_found){
            num_writes_ += 1;
          }
        }, num_writes);

    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      (*used_size_) += num_writes;
    });
     */
    size_type my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));


    if (my_write_index >= max_value_size_) {
      return INSERT_FULL;
    }
    else {

      keys[my_write_index] = key;
      values[my_write_index] = value;
      size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
      if (hashbeginning == -1){
        used_hashes[Kokkos::atomic_fetch_add(used_hash_size, size_type(1))] = hash;
      }
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
  int vector_atomic_insert_into_hash_mergeAdd (
      const team_member_t & teamMember,
      const int vector_size,

      size_type &hash,
      const key_type key,
      const value_type value,
      volatile size_type *used_size_,
      const size_type max_value_size_
      //,bool print = false
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


	    /*
	    Kokkos::parallel_scan(
	        Kokkos::ThreadVectorRange(teamMember, vector_size),
	        [&] (const int threadid, int &update, const bool final) {
	      if (final){
	        write_pos = update;
	      }
	      update += key_not_found;
	    });

	    size_type my_write_index = (*used_size_) + write_pos;
	    int num_writes = 0;
	    Kokkos::parallel_reduce(
	            Kokkos::ThreadVectorRange(teamMember, vector_size),
	            [&] (const int threadid, int &num_writes_) {
	          if (key_not_found){
	            num_writes_ += 1;
	          }
	        }, num_writes);

	    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
	      (*used_size_) += num_writes;
	    });
	     */
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
	      size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
	      hash_nexts[my_write_index] = hashbeginning;
	      return INSERT_SUCCESS;
	    }
	  }

  //used_size should be a shared pointer among the thread vectors
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int team_atomic_insert_into_hash_mergeAdd_TrackHashes (
      const team_member_t & teamMember,
      const int vector_size,

      size_type &hash,
      const key_type key,
      const value_type value,
      volatile size_type *used_size_,
      const size_type max_value_size_,
      size_type *used_hash_size,
      size_type *used_hashes
      //,bool print = false
      ){


    char key_not_found = 1;

    if (hash != -1){
      size_type i = hash_begins[hash];
      for (; i != -1; i = hash_nexts[i]){
        if (keys[i] == key){
          key_not_found = 0;
          values[i] = values[i] + value;

          break;
        }
      }
    }
    else {
      key_not_found = 0;
    }

    if ((*used_size_) >= max_value_size_){
    	return key_not_found; // success means 0, full means 1.
#if 0
      if (key_not_found == 0) {
        return INSERT_SUCCESS;
      }
      else {
        return INSERT_FULL;
      }
#endif
    }

    size_type my_write_index = 0;
    if (key_not_found){
    	my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));
    	//my_write_index = used_size_[0]++;
    }
    else {
    	return INSERT_SUCCESS;
    }


    if (my_write_index >= max_value_size_) {
      return INSERT_FULL;
    }
    else {

      keys[my_write_index] = key;
      values[my_write_index] = value;


      size_type hashbeginning = hash_begins[hash];
      hash_nexts[my_write_index] = hashbeginning;

      hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
      hash_nexts[my_write_index] = hashbeginning;

      if (hashbeginning == -1){
        used_hashes[Kokkos::atomic_fetch_add(used_hash_size, 1)] = hash;
      }
      return INSERT_SUCCESS;
    }
  }


  //function to be called from device.
  //Accumulation is Add operation. It is not atomicAdd, as this
  //is for the cases where we know that none of the simultanous
  //insertions will have the same key.
  //Insertion is simulteanous for the threads of a team
  //used_size should be a shared pointer among the thread
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int team_atomic_insert_into_hash_mergeAdd (
      const team_member_t & teamMember,
      const int vector_size,

      size_type &hash,
      const key_type key,
      const value_type value,
      volatile size_type *used_size_,
      const size_type max_value_size_
      //,bool print = false
      ){


    int key_not_found = 1;

    if (hash != -1){
      size_type i = hash_begins[hash];
#if 0
      std::cout << "hash" <<  hash << " hash_begins[hash]:" << hash_begins[hash] << std::endl;
#endif
      for (; i != -1; i = hash_nexts[i]){
#if 0
      std::cout << "hash" << " hash" << " hash_begins[hash]:" << hash_begins[hash] << " keys[i]:" << keys[i] << std::endl;
#endif
        if (keys[i] == key){
          key_not_found = 0;
          values[i] = values[i] + value;
          break;
        }
      }
    }
    else {
      key_not_found = 0;
    }


    if ((*used_size_) >= max_value_size_){
    	return key_not_found; // success means 0, full means 1.
#if 0
      if (key_not_found == 0) {
        return INSERT_SUCCESS;
      }
      else {
        return INSERT_FULL;
      }
#endif
    }

    size_type my_write_index = 0;
    if (key_not_found){
    	my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));
    	//my_write_index = used_size_[0]++;
    }
    else {
    	return INSERT_SUCCESS;
    }

    if (my_write_index >= max_value_size_) {
      return INSERT_FULL;
    }
    else {
      keys[my_write_index] = key;
      values[my_write_index] = value;


      size_type hashbeginning = hash_begins[hash];
      hash_nexts[my_write_index] = hashbeginning;
      //MD Nov. 2017: I need to do this twice. If I skip above below gets errors.
      //This is because until_hash_next are set, some other threads reads hash_begins and it does not find anything useful.
      //I need to set it again with atomic exchange because of race conditions.
      //It is possible that hash_nexts[my_write_index] might point to a further one in the hashmap.
      //But this is okay as no insertions can have same value at the same time.
      hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
      hash_nexts[my_write_index] = hashbeginning;
      return INSERT_SUCCESS;
    }

	//  return 0;
  }

  //function to be called from device.
  //Accumulation is Add operation.
  //It is an atomicAdd, we might have simultanous insertions with the same key.
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int team_atomic_insert_into_hash_mergeAddAtomic (
      const team_member_t & teamMember,
      const int vector_size,

	  size_type &hash,
      const key_type key,
	  const value_type value,
	  volatile size_type *used_size_,
	  const size_type max_value_size_,
	  const size_type & first_bit,
	  const size_type & rest_bit, bool print = false
      ){

#if 0
	  if (print) std::cout << "\tcol:" << key << " val:" << value << " hash:" << hash << " used_size_:" << *used_size_ << std::endl;
#endif
	//in this one we represent the end of link with 0.
    const size_type end_of_link = 0;


    size_type i = end_of_link;

    if (hash != -1){
      i = hash_begins[hash] & rest_bit;

#if 0
      if (print) std::cout << "\t 1 hash_begins[hash]:" << hash_begins[hash] << std::endl;
#endif

      if (i == end_of_link) {
        //we need to insert it to beginning.
      }
      else if (keys[i] == key){
#if 0
        if (print) std::cout << "\t 2 keys[i]:" << keys[i] << std::endl;
#endif
        Kokkos::atomic_add(values + i, value);
        return 0;
      }
      else if (keys[i] > key){
#if 0
    	  if (print) std::cout << "\t 3 keys[i]:" << keys[i] << " key:" << key << std::endl;
#endif
    	  //we need to insert it to beginning.
    	  i = end_of_link;
    	  //revert i.
      } else {
    	  size_type next = hash_nexts[i] & rest_bit;
    	  for (; next != end_of_link; next = hash_nexts[next] & rest_bit){
#if 0
    		  if (print) std::cout << "\t 4 keys[next] :" << keys[next]  << " key:" << key << std::endl;
#endif
    		  if (keys[next] == key){
    			  Kokkos::atomic_add(values + next, value);
    			  return 0;
    		  }
    		  else if (keys[next] > key){
    			  //this means we need to insert between i and next.
    			  //i will be some positive value.
    			  break;
    		  }
                  i = next;
    	  }
    	  //if we come here, it means we need to insert at the end after i.
    	  //i will be some positive value.
      }
    }
    else {
    	return 0;
    }

    //if we are here, it means we could not insert it.
    //we will try to insert it below.


    while(1){
    	if (i == end_of_link){
#if 0
    		if (print) std::cout << "\t 5 i :" << i << std::endl;
#endif

    		//lock hash_begins[hash];
    		volatile size_type initial_hash_begin =  hash_begins[hash] & rest_bit;
    		volatile size_type locked_hash_begin =  initial_hash_begin | first_bit;
    		//if (hash_begins[hash] == initial_hash_begin){ hash_begins[hash] = locked_hash_begin;
    		if (!Kokkos::atomic_compare_exchange_strong(hash_begins + hash, initial_hash_begin, locked_hash_begin)){
    			continue;
    		}
    		//std::atomic_thread_fence(std::memory_order_acquire);
    		//volatile int tmp = 0;
    		//md note somehow commented part
    		if (initial_hash_begin == end_of_link || (keys[initial_hash_begin] > key)){
    			{
    				volatile size_type my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));
#if 0
    				if (print) std::cout << "\t 8 my_write_index :" << my_write_index << " max_value_size_:" << max_value_size_ << std::endl;
#endif
    				if (my_write_index >= max_value_size_){
    					hash_begins[hash] = initial_hash_begin;
    					return 1;
    				}
    				//Kokkos::volatile_store(&(keys[my_write_index]), key);
    				//Kokkos::volatile_store(&(values[my_write_index]), value);
    				//Kokkos::volatile_store(&(hash_nexts[my_write_index]), initial_hash_begin);

    				keys[my_write_index] = key;
    				values[my_write_index] = value;
    				hash_nexts[my_write_index] = initial_hash_begin;

    				//std::atomic_thread_fence(std::memory_order_release);
    				hash_begins[hash] = my_write_index;
    				//Kokkos::volatile_store(&(hash_begins[hash]), my_write_index);

    				/*
    				if( keys[initial_hash_begin] == key)
    					std::cout 	<< "initial_hash_begin:" << initial_hash_begin
									<< " key:" << key
									<< " keys[initial_hash_begin]:" << keys[initial_hash_begin]
							        << std::endl;
    				 */
    				return 0;

    			}

    		}
#if 0
    		else if ((tmp = Kokkos::volatile_load(&(keys[initial_hash_begin]))) > key) {
				if( keys[initial_hash_begin] == key) std::cout << "1 tmp:" << tmp << " initial_hash_begin:" << initial_hash_begin << " key:" << key << " keys[initial_hash_begin]:" << keys[initial_hash_begin] << std::endl;

#if 0
    			if (print) std::cout << "\t 6 initial_hash_begin :" << initial_hash_begin << std::endl;
    			if (print && initial_hash_begin != end_of_link) std::cout << "\t 7 keys[initial_hash_begin] :" << keys[initial_hash_begin] << " key:" << key<< std::endl;
#endif
    			//we need to insert it to beginning.
    			{
    				volatile size_type my_write_index = Kokkos::atomic_fetch_add(used_size_, 1);
#if 0
    				if (print) std::cout << "\t 8 my_write_index :" << my_write_index << " max_value_size_:" << max_value_size_ << std::endl;
#endif
    				if (my_write_index >= max_value_size_){
    					hash_begins[hash] = initial_hash_begin;
    					return 1;
    				}
    				if( keys[initial_hash_begin] == key) std::cout << "2 tmp:" << tmp << " initial_hash_begin:" << initial_hash_begin << " key:" << key << " keys[initial_hash_begin]:" << keys[initial_hash_begin] << std::endl;


    				//Kokkos::volatile_store(&(keys[my_write_index]), key);
    				//Kokkos::volatile_store(&(values[my_write_index]), value);
    				//Kokkos::volatile_store(&(hash_nexts[my_write_index]), initial_hash_begin);

    				keys[my_write_index] = key;
    				values[my_write_index] = value;
    				hash_nexts[my_write_index] = initial_hash_begin;

    				//std::atomic_thread_fence(std::memory_order_release);
    				//Kokkos::volatile_store(&(hash_begins[hash]), my_write_index);

    				hash_begins[hash] = my_write_index;
    				if( keys[initial_hash_begin] == key) std::cout << " 3 tmp:" << tmp << " initial_hash_begin:" << initial_hash_begin << " key:" << key << " keys[initial_hash_begin]:" << keys[initial_hash_begin] << std::endl;
    				return 0;

    			}
    		}
#endif

    		else if (keys[initial_hash_begin] == key){
#if 0
    			if (print) std::cout << "\t 9 keys[initial_hash_begin] :" << keys[initial_hash_begin] << " key:" << key << std::endl;
#endif
    			Kokkos::atomic_add(values + initial_hash_begin, value);
    			hash_begins[hash] = initial_hash_begin;
    			//release lock
    			return 0;
    		}
    		else {
#if 0
    			if (print) std::cout << "\t 9 keys[initial_hash_begin] :" << keys[initial_hash_begin] << " key:" << key << std::endl;
#endif

    			hash_begins[hash] = initial_hash_begin;
    			//release lock
    			i = initial_hash_begin;
    			size_type next = hash_nexts[initial_hash_begin] & rest_bit;
    			for (; next != end_of_link; next = hash_nexts[next] & rest_bit){
#if 0
    				if (print) std::cout << "\t 10 keys[next] :" << keys[next]  << " key:" << key << std::endl;
#endif

    				if (keys[next] == key){
    					Kokkos::atomic_add(values + next, value);
    					return 0;
    				}
    				else if (keys[next] > key){
    					break;
    				}
    				i = next;
    			}
    			//if we come here, it means we need to insert at the end after i.
    			//i will be some positive value.
    		}


    		//if hash_begin is not - 1 anymore,
    		//loop through it
    		//if found atomic
    		//increment the value and spin = 0;
    		//if not found
    		//find new i.
    		//release the lock.
    		//otherwise
    		//allocate new one, insert and spin 0.
    		//release lock
    	}
#if 0
    	else if (1){
    		return 0;
    	}
#endif
    	else {
    		//lock hash_begin[hash];
    		size_type initial_hash_begin = hash_nexts[i] & rest_bit;
    		size_type locked_hash_begin =  initial_hash_begin | first_bit;


    		if (keys[i] == key){
    			Kokkos::atomic_add(values + i, value);
    			return 0;
    		}
    		//std::atomic_thread_fence(std::memory_order_acquire);


    		if (!Kokkos::atomic_compare_exchange_strong(hash_nexts + i, initial_hash_begin, locked_hash_begin)){
#if 0
    			std::cout << "spinning here2 initial_hash_begin:" << initial_hash_begin << " locked_hash_begin:" << locked_hash_begin << std::endl;
#endif
    			continue;
    		}

    		if (initial_hash_begin == end_of_link || keys[initial_hash_begin] > key) {
#if 0
    			if (print) std::cout << "\t 11 initial_hash_begin :" << initial_hash_begin << std::endl;
    			if (print && initial_hash_begin != end_of_link) std::cout << "\t 12 keys[initial_hash_begin] :" << keys[initial_hash_begin] << " key:" << key<< std::endl;
#endif
    			//we need to insert it to beginning.
    			{
    				size_type my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));
#if 0
    				if (print) std::cout << "\t 13 my_write_index :" << my_write_index << " max_value_size_:" << max_value_size_ << std::endl;
#endif
    				if (my_write_index >= max_value_size_){
    					hash_nexts[i] = initial_hash_begin;
    					return 1;
    				}
    				keys[my_write_index] = key;
    				values[my_write_index] = value;
    				hash_nexts[my_write_index] = initial_hash_begin;

    				//std::atomic_thread_fence(std::memory_order_release);
    				hash_nexts[i] = my_write_index;

    				return 0;
    			}
    		}
    		else if (keys[initial_hash_begin] == key){
#if 0
    			if (print) std::cout << "\t 14 keys[initial_hash_begin] :" << keys[initial_hash_begin] << " key:" << key << std::endl;
#endif
    			Kokkos::atomic_add(values + initial_hash_begin, value);
    			hash_nexts[i] = initial_hash_begin;
    			//release lock
    			return 0;
    		}
    		else {
    			hash_nexts[i] = initial_hash_begin;
    			i = initial_hash_begin;
    			//release lock
    			size_type next = hash_nexts[initial_hash_begin] & rest_bit;
    			for (; next != end_of_link; next = hash_nexts[next] & rest_bit){
#if 0
    				if (print) std::cout << "\t 15 keys[next] :" << keys[next]  << " key:" << key << std::endl;
#endif
    				if (keys[next] == key){
    					Kokkos::atomic_add(values + next, value);
    					return 0;
    				}
    				else if (keys[next] > key){
    					break;
    				}
    				i = next;
    			}
    			//if we come here, it means we need to insert at the end after i.
    			//i will be some positive value.
    		}

    		//lock hash_next[i];
    		//if hash_next[hash_next[i]] is now smaller than my value
    		//loop through it
    		//if found atomic
    		//increment the value and spin = 0;
    		//if not found
    		//find new i.
    		//release the lock.
    		//otherwise
    		//allocate new one, insert and spin 0.
    		//release lock
    	}
    }
	//return 0;
  }


  struct BitwiseOrReduction{
    key_type my_key;
    value_type myval;
    const key_type shared_key;
    size_type &hash;
    const bool am_i_owner;
    KOKKOS_INLINE_FUNCTION
    BitwiseOrReduction(
        key_type my_key_,
        value_type myval_,
        const key_type *shared_key_,
        size_type &hash_, bool am_i_owner_):
        my_key(my_key_),
        myval(myval_),
        shared_key(shared_key_),
        hash(hash_),
        am_i_owner(am_i_owner_){}


    KOKKOS_INLINE_FUNCTION
    void operator()(const int &i, value_type &or_reduction) const {
      if (my_key == *shared_key){
        or_reduction = or_reduction | myval;
        if (!am_i_owner){
          hash = -1;
        }
      }
    }

    KOKKOS_INLINE_FUNCTION
    void join (volatile value_type& dst,const volatile value_type& src) const {
      dst = dst | src;
    }

    KOKKOS_INLINE_FUNCTION
    void init (value_type& dst) const
    {
      dst = 0;
    }
  };

  //function to be called from device.
  //If vector insertion might have the same key values,
  //call this function before.
  //this will do the accumulation on the vectors first, so that
  //each key will be uniquely inserted.
  //then call the insertion functions above.
  //for example if vector lanes is inserting,
  //1 1 1 2 2 2 to the hashmap,
  //the values will be merged (or) and
  //threads will insert 1 x x 2 x x with the merged values.
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  void vector_mergeOr (
      const team_member_t & teamMember,
      const int vector_size,
      const int vector_id,
      size_type &hash,
      const key_type &key,
      value_type &value,
      key_type *shared_key
      ){

    //each vectors result is merged
    for (int i = 0; i < vector_size; ++i){
      bool am_i_owner = false;
      if (i == vector_id){
        am_i_owner = true;
        //if its merged previously, or has nothing to do,
        //then sets shared_key to -1 so that we will skip it.
        if (hash == -1 ){
          *shared_key = -1;
        }
        else {
          *shared_key = key;
          //*shared_value = 0;
        }
      }


      if (*shared_key != -1){

        value_type merged_value = 0;
        Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            BitwiseOrReduction(key, value, shared_key, hash, am_i_owner) , merged_value);
        if (am_i_owner) value = merged_value;
      }
    }
  }


  //function to be called from device.
  //If vector insertion might have the same key values,
  //call this function before.
  //this will do the accumulation on the vectors first, so that
  //each key will be uniquely inserted.
  //then call the insertion functions above.
  //for example if vector lanes is inserting,
  //1 1 1 2 2 2 to the hashmap,
  //the values will be merged (sum) and
  //threads will insert 1 x x 2 x x with the merged values.
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  void vector_mergeAdd (
      const team_member_t & teamMember,
      const int vector_size,
      const int vector_id,
      size_type &hash,
      const key_type key,
      value_type &value,
      volatile key_type *shared_key){

    //each vectors result is merged
    for (int i = 0; i < vector_size; ++i){
      bool am_i_owner = false;
      if (i == vector_id){
        am_i_owner = true;
        //if its merged previously, or has nothing to do,
        //then sets shared_key to -1 so that we will skip it.
        if (hash == -1 ){
          *shared_key = -1;
        }
        else {
          *shared_key = key;
        }
      }

      if (*shared_key != -1){
        value_type merged_value = 0;
        Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, value_type &merged_value_) {
          if (key == *shared_key){
            merged_value_ = merged_value_ + value;
            if (!am_i_owner){
              hash = -1;
            }
          }
        }, merged_value);
        if (am_i_owner) value = merged_value;
      }
    }
  }

  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  void vector_mergeAdd_MEM (
      const team_member_t & teamMember,
      const int vector_size,
      //const int vector_id,
      key_type &key,
      value_type &value,
      volatile key_type *result_keys,
      volatile value_type *result_vals){


    size_type new_hash = key % vector_size;

    Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(teamMember, vector_size),
        [&] (const int threadid) {
      result_keys[threadid] = -1;
      key_type r = -1;
      while (true){
        if (result_keys[new_hash] == key){
          //printf("Adding: new_hash:%d, val:%d\n", new_hash, value);
          Kokkos::atomic_add(result_vals + new_hash, value);
          break;
        }
        else if (result_keys[new_hash] == r){
          if (Kokkos::atomic_compare_exchange_strong(result_keys + new_hash, r, key)){
            result_vals[new_hash] = value;
            break;
          }
        }
        else if (++new_hash == vector_size){
          new_hash = 0;
        }
      }
      value = result_vals[threadid];
      key = result_keys[threadid];

    });
  }

  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  void vector_merge_MEM (
      const team_member_t & teamMember,
      const int vector_size,
      key_type &key,
      volatile key_type *result_keys){


    size_type new_hash = key % vector_size;

    Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(teamMember, vector_size),
        [&] (const int threadid) {
      result_keys[threadid] = -1;
      key_type r = -1;
      while (true){
        if (result_keys[new_hash] == key){
          break;
        }
        else if (result_keys[new_hash] == r){
          if (Kokkos::atomic_compare_exchange_strong(result_keys + new_hash, r, key)){
            break;
          }
        }
        else if (++new_hash == vector_size){
          new_hash = 0;
        }
      }
      key = result_keys[threadid];
    });
  }

  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  void vector_mergeOr_MEM (
      const team_member_t & teamMember,
      const int vector_size,
      //const int vector_id,
      key_type &key,
      value_type &value,
      volatile key_type *result_keys,
      volatile value_type *result_vals
      ){


    size_type new_hash = key % vector_size;

    Kokkos::parallel_for(
        Kokkos::ThreadVectorRange(teamMember, vector_size),
        [&] (const int threadid) {
      result_keys[threadid] = -1;
      key_type r = -1;
      while (key!= -1 && true){
        if (result_keys[new_hash] == key){
          //printf("Adding: new_hash:%d, val:%d\n", new_hash, value);
          Kokkos::atomic_fetch_or(result_vals + new_hash, value);

          break;
        }
        else if (result_keys[new_hash] == r){
          if (Kokkos::atomic_compare_exchange_strong(result_keys + new_hash, r, key)){
            result_vals[new_hash] = value;
            break;
          }
        }
        else if (++new_hash == vector_size){
          new_hash = 0;
        }
      }
      value = result_vals[threadid];
      key = result_keys[threadid];

    });
  }

  //function to be called from device.
  //Directly inserts without checking if the value is in the cache already.
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash(
      const team_member_t & teamMember,
      const int &vector_size,

      const size_type &hash,
      const key_type &key,
      const value_type &value,
      volatile size_type *used_size_,
      const size_type &max_value_size_,
      key_type *vals_counts_smem
      ){

    if ((*used_size_) >= max_value_size_){
      return INSERT_FULL;
    }

    char i_will_write;
    hash == -1 ? i_will_write = 0 : i_will_write = 1;


    int write_pos = 0;
    Kokkos::parallel_scan(
        Kokkos::ThreadVectorRange(teamMember, vector_size),
        [&] (const int threadid, int &update, const bool final) {
      if (final){
        write_pos = update;
      }
      update += i_will_write;
    });

    size_type my_write_index = (*used_size_) + write_pos;
    int num_writes = 0;
    Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &num_writes_) {
          if (i_will_write){
            num_writes_ += 1;
          }
        }, num_writes);

    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      (*used_size_) += num_writes;
    });

    if (i_will_write == 0) return INSERT_SUCCESS;


    if (my_write_index >= max_value_size_) {
      return INSERT_FULL;
    }
    else {

      keys[my_write_index] = key;
      values[my_write_index] = value;
      size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
      hash_nexts[my_write_index] = hashbeginning;
      vals_counts_smem[my_write_index] = 1;
      return INSERT_SUCCESS;
    }
  }

  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_TrackHashes(
      const team_member_t & teamMember,
      const int &vector_size,

      const size_type &hash,
      const key_type &key,
      const value_type &value,
      volatile size_type *used_size_,
      const size_type &max_value_size_,
      size_type *used_hash_size,
      size_type *used_hashes,
      key_type *vals_counts_gmem
      ){

    //if ((*used_size_) >= max_value_size_){
    //  return INSERT_FULL;
    //}

    char i_will_write;
    hash == -1 ? i_will_write = 0 : i_will_write = 1;


    int write_pos = 0;
    Kokkos::parallel_scan(
        Kokkos::ThreadVectorRange(teamMember, vector_size),
        [&] (const int threadid, int &update, const bool final) {
      if (final){
        write_pos = update;
      }
      update += i_will_write;
    });

    size_type my_write_index = (*used_size_) + write_pos;
    int num_writes = 0;
    Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &num_writes_) {
          if (i_will_write){
            num_writes_ += 1;
          }
        }, num_writes);

    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      (*used_size_) += num_writes;
    });

    if (i_will_write == 0) return INSERT_SUCCESS;

/*
    if (my_write_index >= max_value_size_) {
      return INSERT_FULL;
    }
    else {
*/
    keys[my_write_index] = key;
    values[my_write_index] = value;
    size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
    if (hashbeginning == -1){
      used_hashes[Kokkos::atomic_fetch_add(used_hash_size, size_type(1))] = hash;
    }
    hash_nexts[my_write_index] = hashbeginning;
    vals_counts_gmem[my_write_index] = 1;
    return INSERT_SUCCESS;
   // }
  }

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
      //values[my_write_index] = value;
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
      size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
      hash_nexts[my_write_index] = hashbeginning;
      return INSERT_SUCCESS;
    }
  }

  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_mergeAND (
      const team_member_t & teamMember,
      const int &vector_size,

      const size_type &hash,
      const key_type &key,
      const value_type &value,
      key_type *vals_counts
      ){

    char key_not_found = 1;
    if (hash != -1){
      size_type i = hash_begins[hash];
      for (; i != -1; i = hash_nexts[i]){
        if (keys[i] == key){
          key_not_found = 0;
          values[i] = (key_type)values[i] & (key_type)value;
          vals_counts[i] += 1;
          break;
        }
      }
    }
    else {
      key_not_found = 0;
    }

    return key_not_found;
    /*
    if (key_not_found == 1) {
      return INSERT_FULL;
    }
    else {
      return INSERT_SUCCESS;
    }
    */
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
		  //values[my_write_index] = value;
		  size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
		  if (hashbeginning == -1){
			  used_hashes[Kokkos::atomic_fetch_add(used_hash_size, size_type(1))] = hash;
		  }
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
  int vector_atomic_insert_into_hash_merge (
      const team_member_t & teamMember,
      const int vector_size,

      size_type &hash,
      const key_type key,
      size_type *used_size_,
      const size_type max_value_size_
      ){

    char key_not_found = 1;
    if (hash != -1){

      size_type i = hash_begins[hash];

      /*
      if (print  && i >= max_value_size || (i < 0 && i != -1))
      printf("BEGING i:%d hash:%d max_value_size:%d\n", i, hash, max_value_size);
      */

      for (; i != -1; i = hash_nexts[i]){
        //if (i >= max_value_size || i < 0)
        //if (print)
        //printf("next of i:%d hash:%d max_value_size:%d\n", i, hash, max_value_size);
        //if (i < 159 || i < 0){
        if (keys[i] == key){
          key_not_found = 0;
          break;
        }

      }
    }
    else {
      key_not_found = 0;
    }

    int write_pos = 0;
    Kokkos::parallel_scan(
        Kokkos::ThreadVectorRange(teamMember, vector_size),
        [&] (const int threadid, int &update, const bool final) {
      if (final){
        write_pos = update;
      }
      update += key_not_found;
    });

    size_type my_write_index = (*used_size_) + write_pos;
    int num_writes = 0;
    Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &num_writes_) {
          if (key_not_found){
            num_writes_ += 1;
          }
        }, num_writes);

    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      (*used_size_) += num_writes;
    });

    if (key_not_found == 0) return INSERT_SUCCESS;


    if (my_write_index >= max_value_size_) {
      return INSERT_FULL;
    }
    else {

      keys[my_write_index] = key;
      size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
      hash_nexts[my_write_index] = hashbeginning;

    }
    hash = -1;
    return INSERT_SUCCESS;
  }
};



}
}

