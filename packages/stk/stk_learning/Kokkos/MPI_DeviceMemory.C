#include <stk_util/parallel/Parallel.hpp>
#include "mtk_kokkos.h"

void mpi_device(int M, MPI_Comm comm) {
  int numProcs = stk::parallel_machine_size(comm);
  int localProc = stk::parallel_machine_rank(comm);

  if (numProcs != 2) {
    if (localProc == 0) {
      std::cout<<"This test only runs on 2 mpi procs." <<std::endl;
    }
    return;
  }

  typedef Kokkos::View<double*>   ViewVectorType;
  ViewVectorType devx("devx", M);

  //Host mirror
  ViewVectorType::HostMirror x =  Kokkos::create_mirror_view(devx);

  // Initialize x vector on host
  for (int i = 0; i < M; ++i) {
    x( i ) = 11.0 + 5*localProc;
  }

  //Deep copy host view to device view
  Kokkos::deep_copy(devx, x);

  if (localProc == 0) {
    EXPECT_NEAR(11.0, x(0), 1.e-6);
    EXPECT_NEAR(11.0, x(M-1), 1.e-6);
  }
  else {
    EXPECT_NEAR(16.0, x(0), 1.e-6);
    EXPECT_NEAR(16.0, x(M-1), 1.e-6);
  }

  MPI_Status status;

  int mpi_tag = 32767;

  if (localProc == 0) {
    MPI_Send(devx.data(), devx.size(), MPI_DOUBLE, 1, mpi_tag, comm);
  }
  else {
    MPI_Recv(devx.data(), devx.size(), MPI_DOUBLE, 0, mpi_tag, comm, &status);
  }

  //Deep copy device view to host
  Kokkos::deep_copy(x, devx);

  EXPECT_NEAR(11.0, x(0), 1.e-6);
  EXPECT_NEAR(11.0, x(M-1), 1.e-6);
}

TEST(MPI_DeviceMemory, sendrecv) {
  int M = 20;
  mpi_device(M, MPI_COMM_WORLD);
}

