#include <Kokkos_Atomic.hpp>
namespace KokkosKernels{

namespace Experimental{

namespace UnorderedHashmap{




template <typename size_type, typename key_type, typename value_type>
struct HashmapAccumulator{
  size_type hash_key_size;
  size_type max_value_size;
  size_type used_size;

  volatile  size_type *hash_begins;
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

      size_type *used_size,
      const size_type max_value_size){

    size_type i = hash_begins[hash];
    for (; i != -1; i = hash_nexts[i]){
      if (keys[i] == key){
        return INSERT_SUCCESS;
      }
    }

    if (*used_size >= max_value_size) return INSERT_FULL;
    size_type my_index = (*used_size)++;

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

      size_type *used_size,
      const size_type max_value_size){

    size_type i = hash_begins[hash];
    for (; i != -1; i = hash_nexts[i]){
      if (keys[i] == key){
        values[i] = values[i] | value;
        return INSERT_SUCCESS;
      }
    }

    if (*used_size >= max_value_size) return INSERT_FULL;
    size_type my_index = (*used_size)++;

    hash_nexts[my_index] = hash_begins[hash];
    hash_begins[hash] = my_index;
    keys[my_index] = key;
    values[my_index] = value;
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

      size_type *used_size,
      const size_type max_value_size){


    size_type i = hash_begins[hash];
    for (; i != -1; i = hash_nexts[i]){
      if (keys[i] == key){
        values[i] = values[i] + value;
        return INSERT_SUCCESS;
      }
    }

    if (*used_size >= max_value_size) return INSERT_FULL;
    size_type my_index = (*used_size)++;

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
      size_type *used_size,
      const size_type max_value_size){

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

    size_type my_write_index = (*used_size) + write_pos;
    int num_writes = 0;
    Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &num_writes) {
          num_writes += key_not_found;
        }, num_writes);

    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      (*used_size) += num_writes;
    });

    if (key_not_found == 0) return INSERT_SUCCESS;


    if (my_write_index >= max_value_size) {
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
      size_type *used_size,
      const size_type max_value_size){

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

    size_type my_write_index = (*used_size) + write_pos;
    int num_writes = 0;
    Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &num_writes) {
          if (key_not_found){
            num_writes += 1;
          }
        }, num_writes);

    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      (*used_size) += num_writes;
    });

    if (key_not_found == 0) return INSERT_SUCCESS;


    if (my_write_index >= max_value_size) {
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
      volatile size_type *used_size,
      const size_type max_value_size,
      size_type *used_hash_size,
      size_type *used_hashes
      //,bool print = false
      ){


    char key_not_found = 1;

    if (hash != -1){
      size_type i = hash_begins[hash];
      //int count = 0;
      for (; i != -1; i = hash_nexts[i]){
        /*
        if (print && hash == 31 && key == 2436 && (value - 0.197059 < 0.00001)){
          printf("inside --- hash:%d key:%d val:%lf i:%d hash_nexts[i]:%d keys[i]:%d count:%d\n", hash, key, value, i, hash_nexts[i], keys[i], count);
        }

        count++;


        if (print && count > 1000) {
          break;
        }
        */

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


    if ((*used_size) >= max_value_size){
      if (key_not_found == 0) {
        return INSERT_SUCCESS;
      }
      else {
        return INSERT_FULL;
      }
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

    size_type my_write_index = (*used_size) + write_pos;
    int num_writes = 0;
    Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &num_writes) {
          if (key_not_found){
            num_writes += 1;
          }
        }, num_writes);

    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      (*used_size) += num_writes;
    });

    if (key_not_found == 0) return INSERT_SUCCESS;


    if (my_write_index >= max_value_size) {
      return INSERT_FULL;
    }
    else {

      keys[my_write_index] = key;
      values[my_write_index] = value;
      size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
      if (hashbeginning == -1){
        used_hashes[Kokkos::atomic_fetch_add(used_hash_size, 1)] = hash;
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
      volatile size_type *used_size,
      const size_type max_value_size
      //,bool print = false
      ){


    char key_not_found = 1;

    if (hash != -1){
      size_type i = hash_begins[hash];
      //int count = 0;
      for (; i != -1; i = hash_nexts[i]){
        /*
        if (print && hash == 31 && key == 2436 && (value - 0.197059 < 0.00001)){
          printf("inside --- hash:%d key:%d val:%lf i:%d hash_nexts[i]:%d keys[i]:%d count:%d\n", hash, key, value, i, hash_nexts[i], keys[i], count);
        }

        count++;


        if (print && count > 1000) {
          break;
        }
        */

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


    if ((*used_size) >= max_value_size){
      if (key_not_found == 0) {
        return INSERT_SUCCESS;
      }
      else {
        return INSERT_FULL;
      }
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

    size_type my_write_index = (*used_size) + write_pos;
    int num_writes = 0;
    Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &num_writes) {
          if (key_not_found){
            num_writes += 1;
          }
        }, num_writes);

    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      (*used_size) += num_writes;
    });

    if (key_not_found == 0) return INSERT_SUCCESS;


    if (my_write_index >= max_value_size) {
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
            [&] (const int threadid, value_type &merged_value) {
          if (key == *shared_key){
            merged_value = merged_value + value;
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
      volatile size_type *used_size,
      const size_type &max_value_size
      ){

    char key_not_found = 1;

    if (hash != -1){
      size_type i = hash_begins[hash];
      for (; i != -1; i = hash_nexts[i]){
        if (keys[i] == key){
          key_not_found = 0;
          values[i] = (key_type)values[i] | (key_type)value;
          break;
        }
      }
    }
    else {
      key_not_found = 0;
    }


    if ((*used_size) >= max_value_size){
      if (key_not_found == 0) {
        return INSERT_SUCCESS;
      }
      else {
        return INSERT_FULL;
      }
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

    size_type my_write_index = (*used_size) + write_pos;
    int num_writes = 0;
    Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &num_writes) {
          if (key_not_found){
            num_writes += 1;
          }
        }, num_writes);

    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      (*used_size) += num_writes;
    });

    if (key_not_found == 0) return INSERT_SUCCESS;


    if (my_write_index >= max_value_size) {
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
      volatile size_type *used_size,
      const size_type &max_value_size,
      size_type *used_hash_size,
      size_type *used_hashes
      ){

    char key_not_found = 1;

    if (hash != -1){
      size_type i = hash_begins[hash];
      for (; i != -1; i = hash_nexts[i]){
        if (keys[i] == key){
          key_not_found = 0;
          values[i] = (key_type)values[i] | (key_type)value;
          break;
        }
      }
    }
    else {
      key_not_found = 0;
    }


    if ((*used_size) >= max_value_size){
      if (key_not_found == 0) {
        return INSERT_SUCCESS;
      }
      else {
        return INSERT_FULL;
      }
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

    size_type my_write_index = (*used_size) + write_pos;
    int num_writes = 0;
    Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &num_writes) {
          if (key_not_found){
            num_writes += 1;
          }
        }, num_writes);

    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      (*used_size) += num_writes;
    });

    if (key_not_found == 0) return INSERT_SUCCESS;


    if (my_write_index >= max_value_size) {
      return INSERT_FULL;
    }
    else {

      keys[my_write_index] = key;
      values[my_write_index] = value;
      size_type hashbeginning = Kokkos::atomic_exchange(hash_begins+hash, my_write_index);
      if (hashbeginning == -1){
        //Kokkos::atomic_fetch_add(used_hash_size, 1);
        used_hashes[Kokkos::atomic_fetch_add(used_hash_size, 1)] = hash;
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
      size_type *used_size,
      const size_type max_value_size
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

    size_type my_write_index = (*used_size) + write_pos;
    int num_writes = 0;
    Kokkos::parallel_reduce(
            Kokkos::ThreadVectorRange(teamMember, vector_size),
            [&] (const int threadid, int &num_writes) {
          if (key_not_found){
            num_writes += 1;
          }
        }, num_writes);

    Kokkos::single(Kokkos::PerThread(teamMember),[=] () {
      (*used_size) += num_writes;
    });

    if (key_not_found == 0) return INSERT_SUCCESS;


    if (my_write_index >= max_value_size) {
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
}
