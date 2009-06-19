#include "Kokkos_TBBNode.hpp"

tbb::task_scheduler_init Kokkos::TBBNode::tsi_(tbb::task_scheduler_init::deferred);
