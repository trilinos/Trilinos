//@HEADER
// ************************************************************************
//
//                        Kokkos v. 4.0
//       Copyright (2022) National Technology & Engineering
//               Solutions of Sandia, LLC (NTESS).
//
// Under the terms of Contract DE-NA0003525 with NTESS,
// the U.S. Government retains certain rights in this software.
//
// Part of Kokkos, under the Apache License v2.0 with LLVM Exceptions.
// See https://kokkos.org/LICENSE for license information.
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception
//
//@HEADER
#ifndef _KOKKOSKERNELS_HASHMAPACCUMULATOR_HPP
#define _KOKKOSKERNELS_HASHMAPACCUMULATOR_HPP
#include <Kokkos_Atomic.hpp>
#include "KokkosKernels_Macros.hpp"
#include <atomic>

namespace KokkosKernels {

namespace Experimental {

/**
 * @brief types of hash operations supported by HashmapAccumulator.
 *
 * /var bitwiseAnd: Performs key & hashOpRHS
 * /var modulo:     Performs key % hashOpRHS
 * /var pow2Modulo: Performs key & (hashOpRHS - 1)
 */
struct HashOpType {
  struct bitwiseAnd {};
  struct modulo {};
  struct pow2Modulo {};
};

template <typename size_type, typename key_type, typename value_type, typename hash_type>
/**
 * \brief HashmapAccumulator class
 * The use of this is described in the paper:
 *   "Performance-portable sparse matrix-matrix multiplication for many-core
 * architectures" ( https://ieeexplore.ieee.org/abstract/document/7965111/ ) in
 * section III.D
 *
 * Public members:
 * \var hash_begins: Holds the beginning indices of the linked lists
 *                   corresponding to hash values [Begins]
 * \var hash_nexts:  Holds the indicies of the next elements
 *                   within the linked list [Nexts]
 * \var keys:        This stores the column indices of the crs matrix [Ids]
 * \var values:      This store the numerical values (matrix elements) [Values]
 *
 * Private members:
 * \var __max_value_size: The length of the two arrays (keys and hash_nexts)
 * \var __hashOpRHS:      The right hand side of the requested hash operation.
 * \var __insert_success: Value to return upon insertion success.
 * \var __insert_full:    Value to return upon insertion failure.
 */
struct HashmapAccumulator {
  // begin public members
  // issue-508, TODO: It's best for used_size to be an internal member of this
  // class but the current use-cases rely on used_size to be a parameter to the
  // below insertion routines. One way to remove used_size as a parameter to the
  // insertion routines is to instantiate multiple HashmapAccumulator objects
  // (one hashmap for each team of threads) instead of using a single
  // HashmapAccumulator object for multiple teams of threads; this entails
  // major refactoring throughout the kokkos-kernels code base.
  // Making used_size a pointer and private member of this
  // class still exposes access to this member outside of the class and is
  // not a good option.
  // size_type used_size;

  // issue-508, TODO: The hash_begins, hash_nexts, keys, values,
  // __insert_success, and __insert_full members should all be private as well.
  // They should be managed solely by this HashmapAccumulator class: initialized
  // in the constructor(s) and only managed by HashmapAccumulator insertion
  // routines. Making these members private requires major refactoring
  // throughout the kokkos-kernels code base. If allocations for these members
  // must really live outside this class, we need new members that break
  // __max_value_size into: hash_begins_len, hash_nexts_len, keys_len, and
  // values_len...!

  size_type *hash_begins;
  size_type *hash_nexts;
  key_type *keys;
  value_type *values;

  /**
   * \brief default constructor HashmapAccumulator
   * Sets used_size to 0, __insert_success to 0, __insert_full to 1, and
   * __hashOpRHS to 0.
   *
   * Assumption: hash_begins_ are all initialized to -1.
   */
  KOKKOS_INLINE_FUNCTION
  HashmapAccumulator() : hash_begins(), hash_nexts(), keys(), values(), __max_value_size(), __hashOpRHS(0) {}

  /**
   * \brief parameterized constructor HashmapAccumulator
   * Sets used_size to 0, __insert_success to 0, and __insert_full to 1.
   *
   * \param max_value_size_: The length of the two arrays (keys and hash_nexts)
   * \param hashOpRHS:       The right hand side of the requested hash
   * operation. \param hash_begins_:    Holds the beginning indices of the
   * linked lists corresponding to hash values [Begins] \param hash_nexts_:
   * Holds the indicies of the next elements within the linked list [Nexts]
   * \param keys_:           This stores the column indices of (??) [Ids]
   * \param values_:         This store the (matrix element?) numerical value of
   * (??) [Values]
   *
   * Assumption: hash_begins_ are all initialized to -1.
   */
  KOKKOS_INLINE_FUNCTION
  HashmapAccumulator(const size_type max_value_size_, const size_type hashOpRHS, size_type *hash_begins_,
                     size_type *hash_nexts_, key_type *keys_, value_type *values_)
      : hash_begins(hash_begins_),
        hash_nexts(hash_nexts_),
        keys(keys_),
        values(values_),
        __max_value_size(max_value_size_),
        __hashOpRHS(hashOpRHS) {
    // Substract 1 and use the bitwiseAnd __compute_hash member.
    if (std::is_same<hash_type, HashOpType::pow2Modulo>::value) {
      __hashOpRHS -= 1;
    }
  }

  // function to be called from device.
  // Accumulation is OR operation.
  // Insertion is sequential, no race condition for the insertion.
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_mergeOr_TrackHashes(key_type key, value_type value, size_type *used_size_,
                                                      size_type *used_hash_size, size_type *used_hashes) {
    size_type hash, i, my_index;

    if (key == -1) return __insert_success;

    hash = __compute_hash(key, __hashOpRHS);
    for (i = hash_begins[hash]; i != -1; i = hash_nexts[i]) {
      if (keys[i] == key) {
        values[i] = values[i] | value;
        return __insert_success;
      }
    }

    if (*used_size_ >= __max_value_size) return __insert_full;
    my_index = (*used_size_)++;

    if (hash_begins[hash] == -1) {
      used_hashes[used_hash_size[0]++] = hash;
    }
    hash_nexts[my_index] = hash_begins[hash];

    hash_begins[hash] = my_index;
    keys[my_index]    = key;
    values[my_index]  = value;
    return __insert_success;
  }

  // function to be called from device.
  // Accumulation is OR operation.
  // TODO: This function is for triangle counting.
  // Assume that there are 2 values for triangle count.
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_mergeOr_TriangleCount_TrackHashes(key_type key, value_type value, value_type *values2,
                                                                    size_type *used_size_, size_type *used_hash_size,
                                                                    size_type *used_hashes) {
    size_type hash, i, my_index;

    if (key == -1) return __insert_success;

    hash = __compute_hash(key, __hashOpRHS);
    for (i = hash_begins[hash]; i != -1; i = hash_nexts[i]) {
      if (keys[i] == key) {
        values2[i] = values2[i] | (values[i] & value);
        values[i]  = values[i] | value;
        return __insert_success;
      }
    }

    if (*used_size_ >= __max_value_size) return __insert_full;
    my_index = (*used_size_)++;

    if (hash_begins[hash] == -1) {
      used_hashes[used_hash_size[0]++] = hash;
    }
    hash_nexts[my_index] = hash_begins[hash];

    hash_begins[hash] = my_index;
    keys[my_index]    = key;
    values[my_index]  = value;
    values2[my_index] = 0;
    return __insert_success;
  }

  // this is used in slow triangle counting method.
  // L x Incidence
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_mergeAnd_TriangleCount_TrackHashes(key_type key, value_type value,
                                                                     value_type *values2, size_type * /*used_size_*/,
                                                                     size_type * /*used_hash_size*/,
                                                                     size_type * /*used_hashes*/) {
    size_type hash, i;

    if (key == -1) return __insert_success;

    // this function will only try to do an AND operation with
    // existing keys. If the key is not there, returns __insert_full.
    hash = __compute_hash(key, __hashOpRHS);
    for (i = hash_begins[hash]; i != -1; i = hash_nexts[i]) {
      if (keys[i] == key) {
        // values2[i] = values2[i] | (values[i] & value);
        values[i] = values[i] & value;
        ++values2[i];
        return __insert_success;
      }
    }
    return __insert_full;
  }

  // this is used in LxL or Incidence^T x L
  KOKKOS_INLINE_FUNCTION
  value_type sequential_insert_into_hash_mergeAnd_TriangleCount_TrackHashes(key_type key, value_type value) {
    size_type hash, i;

    if (key == -1) return __insert_success;

    // this function will only try to do an AND operation with
    // existing keys. If the key is not there, returns __insert_full.
    hash = __compute_hash(key, __hashOpRHS);
    for (i = hash_begins[hash]; i != -1; i = hash_nexts[i]) {
      if (keys[i] == key) {
        return values[i] & value;
      }
    }
    return 0;
  }

  // this is used in slow triangle counting method.
  // L x Incidence
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_TriangleCount_TrackHashes(key_type key, value_type value, value_type *values2,
                                                            size_type *used_size_, size_type *used_hash_size,
                                                            size_type *used_hashes) {
    size_type hash, my_index;

    if (key == -1) return __insert_success;

    // this function will directly insert, won't check if it exists already.
    if (*used_size_ >= __max_value_size) return __insert_full;
    my_index = (*used_size_)++;

    keys[my_index]    = key;
    values[my_index]  = value;
    values2[my_index] = 1;

    hash = __compute_hash(key, __hashOpRHS);
    if (hash_begins[hash] == -1) {
      hash_begins[hash]                = my_index;
      used_hashes[used_hash_size[0]++] = hash;
    } else {
      hash_nexts[my_index] = hash_begins[hash];
      hash_begins[hash]    = my_index;
    }
    return __insert_success;
  }

  // this is used in LxL or Incidence^T x L
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_TriangleCount_TrackHashes(key_type key, value_type value, size_type *used_size_,
                                                            size_type *used_hash_size,
                                                            size_type *used_hashes)  // issue-508, TODO figure out what
                                                                                     // this "used_hashes" is for
  {
    size_type hash, my_index;

    if (key == -1) return __insert_success;

    // this function will directly insert, won't check if it exists already.
    if (*used_size_ >= __max_value_size) return __insert_full;
    my_index = (*used_size_)++;

    keys[my_index]   = key;
    values[my_index] = value;

    hash = __compute_hash(key, __hashOpRHS);
    if (hash_begins[hash] == -1) {
      hash_begins[hash]                = my_index;
      used_hashes[used_hash_size[0]++] = hash;
    } else {
      hash_nexts[my_index] = hash_begins[hash];
      hash_begins[hash]    = my_index;
    }
    return __insert_success;
  }

  // function to be called from device.
  // Insertion is sequential, no race condition for the insertion.
  // the mergeadd used in the numeric of KKMEM.
  KOKKOS_INLINE_FUNCTION
  void sequential_insert_into_hash_mergeAdd_TrackHashes(key_type key, value_type value, size_type *used_size_,
                                                        size_type *used_hash_size, size_type *used_hashes) {
    size_type hash, i, my_index;

    if (key == -1) return;

    // issue-508, TODO: ensure that i < __max_value_size, but
    // need information about length of keys, values, and hash_nexts first!
    hash = __compute_hash(key, __hashOpRHS);
    for (i = hash_begins[hash]; i != -1; i = hash_nexts[i]) {
      if (keys[i] == key) {
        values[i] = values[i] + value;
        return;
      }
    }

    my_index = (*used_size_)++;

    if (hash_begins[hash] == -1) {
      used_hashes[used_hash_size[0]++] = hash;
    }
    hash_nexts[my_index] = hash_begins[hash];

    hash_begins[hash] = my_index;
    keys[my_index]    = key;
    values[my_index]  = value;
  }

  // no values. simply adds to the keys.
  // used in the compression to count the sets.
  // also used in the symbolic of spgemm if no compression is applied.
  KOKKOS_INLINE_FUNCTION
  int sequential_insert_into_hash_TrackHashes(key_type key, size_type *used_size_, size_type *used_hash_size,
                                              size_type *used_hashes) {
    size_type hash, i, my_index;

    if (key == -1) return __insert_success;

    hash = __compute_hash(key, __hashOpRHS);
    for (i = hash_begins[hash]; i != -1; i = hash_nexts[i]) {
      if (keys[i] == key) {
        return __insert_success;
      }
    }

    my_index = (*used_size_)++;

    if (hash_begins[hash] == -1) {
      used_hashes[used_hash_size[0]++] = hash;
    }
    hash_nexts[my_index] = hash_begins[hash];

    hash_begins[hash] = my_index;
    keys[my_index]    = key;
    return __insert_success;
  }

  // used in the kkmem's numeric phase for second level hashmaps.
  // function to be called from device.
  // Accumulation is Add operation. It is not atomicAdd, as this
  // is for the cases where we know that none of the simultanous
  // insertions will have the same key.
  // Insertion is simulteanous for the vector lanes of a thread.
  // used_size should be a shared pointer among the thread vectors
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_mergeAdd_TrackHashes(const key_type key, const value_type value,
                                                          volatile size_type *used_size_, size_type *used_hash_size,
                                                          size_type *used_hashes) {
    size_type hash, i, my_write_index, hashbeginning;

    if (key == -1) return __insert_success;

    hash = __compute_hash(key, __hashOpRHS);
    if (hash != -1) {
      i = hash_begins[hash];

      for (; i != -1; i = hash_nexts[i]) {
        if (keys[i] == key) {
          values[i] = values[i] + value;
          return __insert_success;
        }
      }
    } else {
      return __insert_success;
    }

    my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));

    if (my_write_index >= __max_value_size) {
      return __insert_full;
    } else {
      keys[my_write_index]   = key;
      values[my_write_index] = value;

#ifdef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
      // this is an issue on VOLTA+ because warps do not go in SIMD fashion
      // anymore. while some thread might insert my_write_index into linked
      // list, another thread in the warp might be reading keys in above loop.
      // before inserting the new value in liked list -- which is done with
      // atomic exchange below, we make sure that the linked is is complete my
      // assigning the hash_next to current head. the head might be different
      // when we do the atomic exchange. this would cause temporarily skipping a
      // key in the linkedlist until hash_nexts is updated second time as below.
      // but this is okay for spgemm, because no two keys will be inserted into
      // hashmap at the same time, as rows have unique columns.

      // Neither the compiler nor the execution unit can re-order the line
      // directly below with the next line performing the atomic_exchange as the
      // atomic exchange writes to hash_begins[hash] and this line reads from
      // hash_begins[hash].
      // This line is needed such that threads of execution can still access the
      // old linked list, after hash_begins+hash has been atomically overwritten
      // with my_write_index but before hash_nexts[my_write_index] is
      // overwritten with hashbeginning. If this line was not here, threads may
      // not be able to access the dangling linked list since
      // hash_nexts[my_write_index] would still be -1.
      hash_nexts[my_write_index] = hash_begins[hash];
#endif

      hashbeginning = Kokkos::atomic_exchange(hash_begins + hash, my_write_index);
      if (hashbeginning == -1) {
        used_hashes[Kokkos::atomic_fetch_add(used_hash_size, size_type(1))] = hash;
      }
      hash_nexts[my_write_index] = hashbeginning;
      return __insert_success;
    }
  }

  // just like vector_atomic_insert_into_hash_mergeAdd_TrackHashes
  // except uses atomic addition on updating the value
  // necessary if duplicate key insertions happen simultaneously
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_mergeAtomicAdd_TrackHashes(const key_type key, const value_type value,
                                                                volatile size_type *used_size_,
                                                                size_type *used_hash_size, size_type *used_hashes) {
    size_type hash, i, my_write_index, hashbeginning;

    if (key == -1) return __insert_success;

    hash = __compute_hash(key, __hashOpRHS);
    if (hash != -1) {
      i = hash_begins[hash];

      for (; i != -1; i = hash_nexts[i]) {
        if (keys[i] == key) {
          Kokkos::atomic_add(values + i, value);
          return __insert_success;
        }
      }
    } else {
      return __insert_success;
    }

    my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));

    if (my_write_index >= __max_value_size) {
      return __insert_full;
    } else {
      keys[my_write_index]   = key;
      values[my_write_index] = value;

#ifdef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
      // this is an issue on VOLTA+ and up because warps do not go in SIMD
      // fashion anymore. while some thread might insert my_write_index into
      // linked list, another thread in the warp might be reading keys in above
      // loop. before inserting the new value in liked list -- which is done
      // with atomic exchange below, we make sure that the linked is is complete
      // my assigning the hash_next to current head. the head might be different
      // when we do the atomic exchange. this would cause temporarily skipping a
      // key in the linkedlist until hash_nexts is updated second time as below.
      // but this is okay for spgemm,
      // because no two keys will be inserted into hashmap at the same time, as
      // rows have unique columns.

      // Neither the compiler nor the execution unit can re-order the line
      // directly below with the next line performing the atomic_exchange as the
      // atomic exchange writes to hash_begins[hash] and this line reads from
      // hash_begins[hash].
      // This line is needed such that threads of execution can still access the
      // old linked list, after hash_begins+hash has been atomically overwritten
      // with my_write_index but before hash_nexts[my_write_index] is
      // overwritten with hashbeginning. If this line was not here, threads may
      // not be able to access the dangling linked list since
      // hash_nexts[my_write_index] would still be -1.
      hash_nexts[my_write_index] = hash_begins[hash];
#endif

      hashbeginning = Kokkos::atomic_exchange(hash_begins + hash, my_write_index);
      if (hashbeginning == -1) {
        used_hashes[Kokkos::atomic_fetch_add(used_hash_size, size_type(1))] = hash;
      }
      hash_nexts[my_write_index] = hashbeginning;
      return __insert_success;
    }
  }

  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_mergeAdd_TrackHashes_no_list(const key_type key, const value_type value,
                                                                  size_type *used_hash_size, size_type *used_hashes) {
    size_type hash;

    if (key == -1) return __insert_success;

    hash = __compute_hash(key, __hashOpRHS);
    if (hash != -1) {
      // loop until an empty hash is found and the key insertion succeeds
      // if our capacity is at least some constant multiple of the current used
      // hashes then the expected number of iterations is constant
      int depth = 0;
      // add key to hash to ensure no two keys follow the same paths over hashes
      // add depth to prevent cycles
      for (;; hash = __compute_hash(hash + key + depth++, __hashOpRHS)) {
        if (keys[hash] == key) {
          Kokkos::atomic_add(values + hash, value);
          return __insert_success;
        } else if (keys[hash] == -1) {
          if (Kokkos::atomic_compare_exchange_strong<key_type>(keys + hash, -1, key)) {
            // should only be here if we used a new hash
            used_hashes[Kokkos::atomic_fetch_add(used_hash_size, size_type(1))] = hash;
            Kokkos::atomic_add(values + hash, value);
            return __insert_success;
          }
          // we don't care if we failed if some other thread succeeded with the
          // same key as ours
          if (keys[hash] == key) {
            Kokkos::atomic_add(values + hash, value);
            return __insert_success;
          }
        }
      }
    } else {
      return __insert_success;
    }
  }

  // NOTE: this is an exact copy of vector_atmoic_insert_into_hash_mergeAdd from
  // https://github.com/kokkos/kokkos-kernels/blob/750fe24508a69ed4dba92bb4a9e17a6094b1a083/src/common/KokkosKernels_HashmapAccumulator.hpp#L442-L502
  template <typename team_member_t>
  KOKKOS_INLINE_FUNCTION int vector_atomic_insert_into_hash_mergeAdd_with_team_level_list_length(
      const team_member_t & /* teamMember */, const int /* vector_size */, size_type hash, const key_type key,
      const value_type value, volatile size_type *used_size_, const size_type max_value_size_) {
    // Cannot compute hash here due to impl_speed use-case
    // hash = __compute_hash(key, __hashOpRHS);
    if (key == -1) return __insert_success;

    if (hash != -1) {
      size_type i = hash_begins[hash];
      for (; i != -1; i = hash_nexts[i]) {
        if (keys[i] == key) {
          values[i] = values[i] + value;
          return __insert_success;
        }
      }
    } else {
      return __insert_success;
    }

    // Ensure that threads don't continue incrementing used_size_ if the hashmap
    // is full, used_size_ could overflow and result in undefined behavior.
    if (used_size_[0] >= max_value_size_) {
      return __insert_full;
    }
    size_type my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));

    if (my_write_index >= max_value_size_) {
      return __insert_full;
    } else {
      keys[my_write_index]   = key;
      values[my_write_index] = value;

#ifdef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
      // this is an issue on VOLTA+ and up because warps do not go in SIMD
      // fashion anymore. while some thread might insert my_write_index into
      // linked list, another thread in the warp might be reading keys in above
      // loop. before inserting the new value in liked list -- which is done
      // with atomic exchange below, we make sure that the linked is is complete
      // my assigning the hash_next to current head. the head might be different
      // when we do the atomic exchange. this would cause temporarily skipping a
      // key in the linkedlist until hash_nexts is updated second time as below.
      // but this is okay for spgemm,
      // because no two keys will be inserted into hashmap at the same time, as
      // rows have unique columns.

      // Neither the compiler nor the execution unit can re-order the line
      // directly below with the next line performing the atomic_exchange as the
      // atomic exchange writes to hash_begins[hash] and this line reads from
      // hash_begins[hash].
      // This line is needed such that threads of execution can still access the
      // old linked list, after hash_begins+hash has been atomically overwritten
      // with my_write_index but before hash_nexts[my_write_index] is
      // overwritten with hashbeginning. If this line was not here, threads may
      // not be able to access the dangling linked list since
      // hash_nexts[my_write_index] would still be -1.
      hash_nexts[my_write_index] = hash_begins[hash];
#endif

      // Atomically:
      // hashbeginning = hash_begins[hash]
      // hash_begins[hash] = my_write_index
      // hash_nexts[my_write_index] = hash_begins[hash]
      size_type hashbeginning    = Kokkos::atomic_exchange(hash_begins + hash, my_write_index);
      hash_nexts[my_write_index] = hashbeginning;
      return __insert_success;
    }
  }

  // used in kkmem's numeric phase to insert to first level hashmaps.
  // function to be called from device.
  // Accumulation is Add operation. It is not atomicAdd, as this
  // is for the cases where we know that none of the simultanous
  // insertions will have the same key.
  // Insertion is simulteanous for the vector lanes of a thread.
  // used_size should be a shared pointer among the thread vectors
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_mergeAdd(const key_type key, const value_type value,
                                              volatile size_type *used_size_) {
    if (key == -1) return __insert_success;

    return vector_atomic_insert_into_hash_mergeAdd_with_team_level_list_length(
        nullptr, 0, __compute_hash(key, __hashOpRHS), key, value, used_size_, __max_value_size);
  }

  // used in symbolic of kkmem if the compression is not applied.
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash(const key_type &key, volatile size_type *used_size_) {
    size_type hash, i, my_write_index, hashbeginning;

    if (key == -1) return __insert_success;

    hash = __compute_hash(key, __hashOpRHS);
    for (i = hash_begins[hash]; i != -1; i = hash_nexts[i]) {
      if (keys[i] == key) {
        return __insert_success;
      }
    }

    my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));

    if (my_write_index >= __max_value_size) {
      return __insert_full;
    } else {
      keys[my_write_index] = key;

#ifdef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
      // this is an issue on VOLTA+ and up because warps do not go in SIMD
      // fashion anymore. while some thread might insert my_write_index into
      // linked list, another thread in the warp might be reading keys in above
      // loop. before inserting the new value in liked list -- which is done
      // with atomic exchange below, we make sure that the linked is is complete
      // my assigning the hash_next to current head. the head might be different
      // when we do the atomic exchange. this would cause temporarily skipping a
      // key in the linkedlist until hash_nexts is updated second time as below.
      // but this is okay for spgemm,
      // because no two keys will be inserted into hashmap at the same time, as
      // rows have unique columns.
      hash_nexts[my_write_index] = hash_begins[hash];
#endif

      hashbeginning              = Kokkos::atomic_exchange(hash_begins + hash, my_write_index);
      hash_nexts[my_write_index] = hashbeginning;
      return __insert_success;
    }
  }

  // function to be called from device.
  // Accumulation is Add operation. It is not atomicAdd, as this
  // is for the cases where we know that none of the simultanous
  // insertions will have the same key.
  // Insertion is simulteanous for the vector lanes of a thread.
  // used_size should be a shared pointer among the thread vectors
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_mergeOr(const key_type &key, const value_type &value,
                                             volatile size_type *used_size_) {
    size_type hash, i, my_write_index, hashbeginning;

    if (key == -1) return __insert_success;

    hash = __compute_hash(key, __hashOpRHS);
    for (i = hash_begins[hash]; i != -1; i = hash_nexts[i]) {
      if (keys[i] == key) {
        values[i] = values[i] | value;
        return __insert_success;
      }
    }

    my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));

    if (my_write_index >= __max_value_size) {
      return __insert_full;
    } else {
      keys[my_write_index]   = key;
      values[my_write_index] = value;

#ifdef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
      // this is an issue on VOLTA+ and up because warps do not go in SIMD
      // fashion anymore. while some thread might insert my_write_index into
      // linked list, another thread in the warp might be reading keys in above
      // loop. before inserting the new value in liked list -- which is done
      // with atomic exchange below, we make sure that the linked is is complete
      // my assigning the hash_next to current head. the head might be different
      // when we do the atomic exchange. this would cause temporarily skipping a
      // key in the linkedlist until hash_nexts is updated second time as below.
      // but this is okay for spgemm,
      // because no two keys will be inserted into hashmap at the same time, as
      // rows have unique columns.
      hash_nexts[my_write_index] = hash_begins[hash];
#endif

      hashbeginning              = Kokkos::atomic_exchange(hash_begins + hash, my_write_index);
      hash_nexts[my_write_index] = hashbeginning;
      return __insert_success;
    }
  }

  // function to be called from device.
  // Accumulation is Add operation. It is not atomicAdd, as this
  // is for the cases where we know that none of the simultanous
  // insertions will have the same key.
  // Insertion is simulteanous for the vector lanes of a thread.
  // used_size should be a shared pointer among the thread vectors
  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_mergeOr_TrackHashes(const key_type &key, const value_type &value,
                                                         volatile size_type *used_size_, size_type *used_hash_size,
                                                         size_type *used_hashes) {
    size_type hash, i, my_write_index, hashbeginning;

    if (key == -1) return __insert_success;

    hash = __compute_hash(key, __hashOpRHS);
    for (i = hash_begins[hash]; i != -1; i = hash_nexts[i]) {
      if (keys[i] == key) {
        values[i] = values[i] | value;
        return __insert_success;
      }
    }

    my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));

    if (my_write_index >= __max_value_size) {
      return __insert_full;
    } else {
      keys[my_write_index]   = key;
      values[my_write_index] = value;

#ifdef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
      // this is an issue on VOLTA+ and up because warps do not go in SIMD
      // fashion anymore. while some thread might insert my_write_index into
      // linked list, another thread in the warp might be reading keys in above
      // loop. before inserting the new value in liked list -- which is done
      // with atomic exchange below, we make sure that the linked is is complete
      // my assigning the hash_next to current head. the head might be different
      // when we do the atomic exchange. this would cause temporarily skipping a
      // key in the linkedlist until hash_nexts is updated second time as below.
      // but this is okay for spgemm,
      // because no two keys will be inserted into hashmap at the same time, as
      // rows have unique columns.
      hash_nexts[my_write_index] = hash_begins[hash];
#endif

      hashbeginning = Kokkos::atomic_exchange(hash_begins + hash, my_write_index);
      if (hashbeginning == -1) {
        used_hashes[Kokkos::atomic_fetch_add(used_hash_size, size_type(1))] = hash;
      }
      hash_nexts[my_write_index] = hashbeginning;
      return __insert_success;
    }
  }

  KOKKOS_INLINE_FUNCTION
  int vector_atomic_insert_into_hash_TrackHashes(const key_type &key, volatile size_type *used_size_,
                                                 size_type *used_hash_size, size_type *used_hashes) {
    size_type hash, i, my_write_index, hashbeginning;

    if (key == -1) return __insert_success;

    hash = __compute_hash(key, __hashOpRHS);
    for (i = hash_begins[hash]; i != -1; i = hash_nexts[i]) {
      if (keys[i] == key) {
        return __insert_success;
      }
    }

    my_write_index = Kokkos::atomic_fetch_add(used_size_, size_type(1));

    if (my_write_index >= __max_value_size) {
      return __insert_full;
    } else {
      keys[my_write_index] = key;

#ifdef KOKKOSKERNELS_CUDA_INDEPENDENT_THREADS
      // this is an issue on VOLTA+ and up because warps do not go in SIMD
      // fashion anymore. while some thread might insert my_write_index into
      // linked list, another thread in the warp might be reading keys in above
      // loop. before inserting the new value in liked list -- which is done
      // with atomic exchange below, we make sure that the linked is is complete
      // my assigning the hash_next to current head. the head might be different
      // when we do the atomic exchange. this would cause temporarily skipping a
      // key in the linkedlist until hash_nexts is updated second time as below.
      // but this is okay for spgemm,
      // because no two keys will be inserted into hashmap at the same time, as
      // rows have unique columns.
      hash_nexts[my_write_index] = hash_begins[hash];
#endif

      hashbeginning = Kokkos::atomic_exchange(hash_begins + hash, my_write_index);
      if (hashbeginning == -1) {
        used_hashes[Kokkos::atomic_fetch_add(used_hash_size, size_type(1))] = hash;
      }
      hash_nexts[my_write_index] = hashbeginning;
      return __insert_success;
    }
  }
  // end public members
 private:
  size_type __max_value_size;
  size_type __hashOpRHS;
  static constexpr int __insert_success = 0;
  static constexpr int __insert_full    = 1;

  template <typename U = hash_type, typename std::enable_if<std::is_same<U, HashOpType::bitwiseAnd>::value ||
                                                                std::is_same<U, HashOpType::pow2Modulo>::value,
                                                            std::size_t>::type = 0>
  KOKKOS_INLINE_FUNCTION int __compute_hash(size_type key, size_type bitmask) {
    size_type hash = key & bitmask;
#ifdef HASHMAPACCUMULATOR_ASSERT_ENABLED
    if (hash == -1) Kokkos::abort("__compute_hash: hash = -1");
    if (key == -1) Kokkos::abort("__compute_hash: key = -1");
#endif  // HASHMAPACCUMULATOR_ASSERT_ENABLED
    return hash;
  }

  template <typename U                                                                             = hash_type,
            typename std::enable_if<std::is_same<U, HashOpType::modulo>::value, std::size_t>::type = 0>
  KOKKOS_INLINE_FUNCTION int __compute_hash(size_type key, size_type divisor) {
    size_type hash = key % divisor;
#ifdef HASHMAPACCUMULATOR_ASSERT_ENABLED
    if (hash == -1) Kokkos::abort("__compute_hash: hash = -1");
    if (key == -1) Kokkos::abort("__compute_hash: key = -1");
#endif  // HASHMAPACCUMULATOR_ASSERT_ENABLED
    return hash;
  }
  // private
};  // struct HashmapAccumulator

}  // namespace Experimental
}  // namespace KokkosKernels

#endif  //  _KOKKOSKERNELS_HASHMAPACCUMULATOR_HPP
