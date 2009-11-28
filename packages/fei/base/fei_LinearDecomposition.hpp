#ifndef _fei_LinearDecomposition_hpp_
#define _fei_LinearDecomposition_hpp_

#include <fei_macros.hpp>
#include <fei_CommUtils.hpp>

namespace fei {

template<typename GlobalIDType>
class LinearDecomposition {
public:
  LinearDecomposition(int localProc, int numProcs,
                      GlobalIDType lowest_global_id,
                      GlobalIDType highest_global_id);
  ~LinearDecomposition() {}

  GlobalIDType first_local_id() const {return first_local;}
  GlobalIDType last_local_id() const {return last_local;}

  GlobalIDType first_global_id() const {return first_global;}
  GlobalIDType last_global_id() const {return last_global;}

  /** Return the mpi rank of the proc that holds 'id' in the linear-decomposition.
    Returns -1 if id < first_global_id() or id > last_global_id().
   */
  int which_proc(GlobalIDType id) const;

private:
  GlobalIDType first_global;
  GlobalIDType last_global;
  GlobalIDType first_local;
  GlobalIDType last_local;
  std::vector<GlobalIDType> proc_offsets;
};//class LinearDecomposition

template<typename GlobalIDType>
LinearDecomposition<GlobalIDType>::LinearDecomposition(int localProc, int numProcs,
                                     GlobalIDType lowest_global_id,
                                     GlobalIDType highest_global_id)
 : first_local(0),
   last_local(0),
   proc_offsets()
{
  GlobalIDType num_global = highest_global_id - lowest_global_id + 1;
  GlobalIDType num_local = num_global/numProcs;
  GlobalIDType remainder = num_global%numProcs;

  //First have each entry in proc_offsets contain the number of local ids:
  proc_offsets.assign(numProcs, num_local);
  for(GlobalIDType i=0; i<remainder; ++i) {
    ++proc_offsets[i];
  }

  //Now convert proc_offsets so that proc_offsets[i] is the i-th proc's
  //offset into the global space of ids:
  GlobalIDType offset = 0;
  for(size_t i=0; i<proc_offsets.size(); ++i) {
    GlobalIDType tmp = proc_offsets[i];
    proc_offsets[i] = offset;
    offset += tmp;
  }

  first_global = lowest_global_id;
  last_global = highest_global_id;
  first_local = lowest_global_id  + proc_offsets[localProc];
  last_local  = highest_global_id + proc_offsets[localProc] + num_local - 1;
}

template<typename GlobalIDType>
int LinearDecomposition<GlobalIDType>::which_proc(GlobalIDType id) const
{
  if (id < first_global || id > last_global) return -1;

  for(size_t i=1; i<proc_offsets.size(); ++i) {
    if (first_global+proc_offsets[i] > id) return i-1;
  }

  int last_proc = proc_offsets.size() - 1;
  return last_proc;
}

}//namespace fei

#endif

