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

int get_rank_on_other_comm(MPI_Comm comm1, MPI_Comm comm2, int rankOnComm1)
{
  MPI_Group comm1Group, comm2Group;
  MPI_Comm_group(comm1, &comm1Group);
  MPI_Comm_group(comm2, &comm2Group);

  int rankOnComm2;
  MPI_Group_translate_ranks(comm1Group, 1, &rankOnComm1, comm2Group, &rankOnComm2);

  return rankOnComm2;
}

} // namespace impl

} // namespace utils
} // namespace middle_mesh
} // namespace stk
