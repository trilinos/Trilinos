#include "Kokkos_TBBNode.hpp"

tbb::task_scheduler_init TBBNode::tsi_(tbb::task_scheduler_init::deferred);
