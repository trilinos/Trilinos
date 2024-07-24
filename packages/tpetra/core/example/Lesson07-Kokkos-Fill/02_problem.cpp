// @HEADER
// *****************************************************************************
//          Tpetra: Templated Linear Algebra Services Package
//
// Copyright 2008 NTESS and the Tpetra contributors.
// SPDX-License-Identifier: BSD-3-Clause
// *****************************************************************************
// @HEADER

#include <Tpetra_Core.hpp>
// This is the only header file you need to include for the "core"
// part of Kokkos.  That includes Kokkos::View, Kokkos::parallel_*,
// and atomic updates.
#include <Kokkos_Core.hpp>

int main (int argc, char* argv[]) {
  // Once we fill into Tpetra data structures, we will need to respect
  // Tpetra's choices of local (LO) and global (GO) indices.  For now,
  // we define these to int.
  typedef int LO;
  typedef int GO;

  Tpetra::ScopeGuard tpetraScope (&argc, &argv);
  {
    auto comm = Tpetra::getDefaultComm ();
    const int numProcs = comm->getSize ();
    int numLclElements = 10000;
    int numGblElements = numProcs * numLclElements;

    // We happened to choose a discretization with the same number of
    // nodes and degrees of freedom as elements.  That need not always
    // be the case.
    int numLclNodes = numLclElements;

    // Describe the physical problem (heat equation with nonhomogeneous
    // Dirichlet boundary conditions) and its discretization.
    Kokkos::View<double*> temperature ("temperature", numLclNodes);
    const double diffusionCoeff = 1.0;
    const double x_left = 0.0; // position of the left boundary
    const double T_left = 0.0; // temperature at the left boundary
    const double x_right = 1.0; // position of the right boundary
    const double T_right = 1.0; // temperature at the right boundary
    const double dx = (x_right - x_left) / numGblElements;
    Kokkos::View<double*> forcingTerm ("forcing term", numLclNodes);

    // Set the forcing term.  We picked it so that we can know the exact
    // solution of the heat equation, namely
    //
    // u(x) = -4.0 * (x - 0.5) * (x - 0.5) + 1.0 + (T_left - T_right)x.
    Kokkos::parallel_for (numLclNodes, KOKKOS_LAMBDA (const LO node) {
	const double x = x_left + node * dx;

	forcingTerm(node) = -8.0;
	// We multiply dx*dx into the forcing term, so the matrix's
	// entries don't need to know it.
	forcingTerm(node) *= dx * dx;
      });

    //
    // TODO implement the discretization, put the discretization into
    // Tpetra data structures, solve the linear system, and correct the
    // solution for the nonhomogenous Dirichlet boundary conditions.
    //
  }
  return 0;
}
