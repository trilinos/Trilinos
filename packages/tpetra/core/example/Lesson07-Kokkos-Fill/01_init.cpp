#include <Teuchos_DefaultMpiComm.hpp> // or mpi.h
// This is the only header file you need to include for the "core"
// part of Kokkos.  That includes Kokkos::View, Kokkos::parallel_*,
// and atomic updates.
#include <Kokkos_Core.hpp>

int main (int argc, char* argv[]) {
  // Initialize MPI and Kokkos.  Both take command line arguments.
  // Kokkos should be initialized after MPI so that processes are
  // correctly bound to cores.
  MPI_Init (&argc, &argv);
  Kokkos::initialize (argc, argv);

  //
  // TODO define the problem and discretization, implement the
  // discretization, put the discretization into Tpetra data
  // structures, solve the linear system, and correct the solution for
  // the nonhomogenous Dirichlet boundary conditions.
  //

  // Shut down Kokkos and MPI, in reverse order from their initialization.
  Kokkos::finalize ();
  MPI_Finalize ();
}
