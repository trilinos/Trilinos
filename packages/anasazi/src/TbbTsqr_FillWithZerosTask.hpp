#ifndef __TSQR_TBB_FillWithZerosTask_hpp
#define __TSQR_TBB_FillWithZerosTask_hpp

#include <tbb/task.h>
#include <TSQR/TBB/TbbTsqr_Partitioner.hpp>
#include <TSQR/Tsqr_SequentialTsqr.hpp>

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

namespace TSQR {
  namespace TBB {
    
    template< class LocalOrdinal, class Scalar >
    class FillWithZerosTask : public tbb::task {
    private:
      typedef MatView< LocalOrdinal, Scalar > mat_view;
      typedef std::pair< mat_view, mat_view > split_type;

    public:
      FillWithZerosTask (const size_t P_first, 
			 const size_t P_last,
			 MatView< LocalOrdinal, Scalar > C,
			 const SequentialTsqr< LocalOrdinal, Scalar >& seq,
			 const bool contiguous_cache_blocks = false)
	: P_first_ (P_first), 
	  P_last_ (P_last), 
	  C_ (C),
	  seq_ (seq), 
	  contiguous_cache_blocks_ (contiguous_cache_blocks)
      {}

      tbb::task* 
      execute () 
      {
	if (P_first_ > P_last_ || C_.empty())
	  return NULL;
	else if (P_first_ == P_last_)
	  {
	    // Fill my partition with zeros.
	    seq_.fill_with_zeros (C_.nrows(), C_.ncols(), C_.get(), 
				  C_.lda(), contiguous_cache_blocks_);
	    return NULL;
	  }
	else
	  {
	    // "c": continuation task
	    tbb::empty_task& c = *new( allocate_continuation() ) tbb::empty_task;

	    // Recurse on two intervals: [P_first, P_mid] and [P_mid+1, P_last]
	    const size_t P_mid = (P_first_ + P_last_) / 2;
	    split_type C_split = 
	      partitioner_.split (C_, P_first_, P_mid, P_last_, 
				  contiguous_cache_blocks_);
	    FillWithZerosTask& topTask = *new( c.allocate_child() )
	      FillWithZerosTask (P_first_, P_mid, C_split.first, seq_, 
				 contiguous_cache_blocks_);
	    FillWithZerosTask& botTask = *new( c.allocate_child() )
	      FillWithZerosTask (P_mid+1, P_last_, C_split.second, seq_, 
				 contiguous_cache_blocks_);

	    // Set reference count of parent (in this case, the
	    // continuation task) to 2 (since 2 children -- no
	    // additional task since no waiting).
	    c.set_ref_count (2);
	    c.spawn (botTask);
	    return &topTask; // scheduler bypass optimization
	  }
      }

    private:
      size_t P_first_, P_last_;
      mat_view C_;
      SequentialTsqr< LocalOrdinal, Scalar > seq_;
      Partitioner< LocalOrdinal, Scalar > partitioner_;
      bool contiguous_cache_blocks_;
    };
  } // namespace TBB
} // namespace TSQR


#endif // __TSQR_TBB_FillWithZerosTask_hpp
