#include "utils.hpp"

namespace stk {
namespace middle_mesh {
namespace utils {
namespace impl {

int comm_rank(MPI_Comm comm)
{
  int myrank;
  MPI_Comm_rank(comm, &myrank);
  return myrank;
}

int comm_size(MPI_Comm comm)
{
  int commSize;
  MPI_Comm_size(comm, &commSize);
  return commSize;
}

MPI_Comm comm_dup(MPI_Comm comm)
{
  MPI_Comm newComm;
  MPI_Comm_dup(comm, &newComm);

  return newComm;
}

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
