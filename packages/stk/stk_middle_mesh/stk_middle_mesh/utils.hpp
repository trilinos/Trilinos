#ifndef TRIANGULATOR_UTILS_H
#define TRIANGULATOR_UTILS_H

#include "mpi.h"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

int comm_rank(MPI_Comm comm);

int comm_size(MPI_Comm comm);

MPI_Comm comm_dup(MPI_Comm);

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
#endif
