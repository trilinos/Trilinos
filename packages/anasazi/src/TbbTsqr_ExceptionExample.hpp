

void example()
{
  try {
    tbb::parallel_for (blocked_range<size_t>(0, ncores, 1), 
		       ParallelTask (...),
		       tbb::simple_partitioner());
  } 
  // TBB can't guarantee on all systems that an exception thrown in
  // another thread will have its type correctly propagated to this
  // thread.  If it can't, then it captures the exception as a
  // tbb:captured_exception, and propagates it to here.  It may be
  // able to propagate the exception, though, so be prepared for that.
  catch (tbb::captured_exception& ex) {
    throw std::runtime_error (ex.what()); // just catch and rethrow
  }
}
